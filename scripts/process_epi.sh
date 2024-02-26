#!/bin/bash

usage="$(basename "$0") [-h] [-t T1_FILE] [-e] DIFF_FILE PA_FILE -- Axon diameter estimation from EPI data

where:
    DIFF_FILE    EPI diffusion images
    PA_FILE      EPI with inverted PE
    -h           show this help text
    -t T1_FILE   Do white matter masking with T1 data
    -e           Use maximum likelihood estimator instead of denoising & Rician bias correction"

# option
mle_flag=false
while getopts "ht:me" opt; do
    case $opt in
        h) echo "$usage"
        exit
        ;;
        t) T1_FILE=${OPTARG} ;;
        e) mle_flag=true ;;
        \?) printf "illegal option: -%s\n" "$OPTARG" >&2
        echo "$usage" >&2
        exit 1
        ;;
    esac
done
shift "$((OPTIND-1))"

# input
if [ "$#" -lt 1 ]; then
    echo "AP file missing"
    exit 1
elif [ "$#" -lt 2 ]; then
    echo "PA file missing"
    exit 1
else
    IN_FILE="$1"
    PA_FILE="$2"
fi

IN_FILE_PREFIX=${IN_FILE%%.*}
IN_FILE_PATH=$(dirname $IN_FILE)
SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

if $mle_flag; then
    # Estimate noise for spherical harmonic decomposition
    dwiextract -force -fslgrad "${IN_FILE_PREFIX}.bvec" "${IN_FILE_PREFIX}.bval" -shells 0,6000 ${IN_FILE} "${IN_FILE_PREFIX}_mag_lowb.nii.gz"
    dwidenoise -force -noise "${IN_FILE_PREFIX}_noise_map.nii.gz" "${IN_FILE_PREFIX}_mag_lowb.nii.gz"  "${IN_FILE_PREFIX}_mag_lowb_denoised.nii.gz"
    cp ${IN_FILE} "${IN_FILE_PREFIX}_denoised_mag.nii.gz"
else
    dwidenoise -force -noise "${IN_FILE_PREFIX}_noise_map.nii.gz" ${IN_FILE} "${IN_FILE_PREFIX}_denoised_mag.nii.gz"
fi

# Convert nii to mif
mrconvert -force "${IN_FILE_PREFIX}_denoised_mag.nii.gz" -fslgrad "${IN_FILE_PREFIX}.bvec" "${IN_FILE_PREFIX}.bval" "${IN_FILE_PREFIX}_denoised_mag.mif"

# Gibbs-Ringing removal
mrdegibbs -force "${IN_FILE_PREFIX}_denoised_mag.mif" "${IN_FILE_PREFIX}_denoised_mag_gr.mif"

# Brain mask for eddy
dwiextract -force "${IN_FILE_PREFIX}_denoised_mag_gr.mif" - -bzero | mrmath -force - mean "${IN_FILE_PREFIX}_denoised_mag_meanb0.mif" -axis 3
mrconvert -force "${IN_FILE_PREFIX}_denoised_mag_meanb0.mif" "${IN_FILE_PREFIX}_denoised_mag_meanb0.nii.gz"
bet "${IN_FILE_PREFIX}_denoised_mag_meanb0.nii.gz" "${IN_FILE_PREFIX}_denoised_mag_meanb0_bet.nii.gz" -f 0.4 -m
gunzip -f "${IN_FILE_PREFIX}_denoised_mag_meanb0_bet_mask.nii.gz" # unzip as otherwise the mask gets corrupted by mrconvert in dwifslpreproc

# this file contains a list of the simultaneously acquired slices in acquisition order
slspec="$SCRIPTPATH/slspec_epi.txt"

# input for topup
dwiextract -force "${IN_FILE}" - -bzero -fslgrad "${IN_FILE_PREFIX}.bvec" "${IN_FILE_PREFIX}.bval" | mrmath -force - mean "${IN_FILE_PREFIX}_AP_b0.mif" -axis 3
mrmath -force "${PA_FILE}" mean "${IN_FILE_PREFIX}_PA_b0.mif" -axis 3
mrcat -force "${IN_FILE_PREFIX}_AP_b0.mif" "${IN_FILE_PREFIX}_PA_b0.mif" "${IN_FILE_PREFIX}_b0.mif" -axis 3

# let mrtrix take care of providing eddy and topup input and output
# readout_time is copied from the json file for our data
# mporder is recommended to be somewhere between N/4 and N/2, where N is the number of excitations
# "dont_sep_offs_move" is used as it leads to reduced shifts between the two shells
SCRATCH_DIR="$IN_FILE_PATH/dwifslpreproc_tmp"
dwifslpreproc -force "${IN_FILE_PREFIX}_denoised_mag_gr.mif" "${IN_FILE_PREFIX}_moco.mif" -nocleanup -scratch "$SCRATCH_DIR" -rpe_pair -se_epi "${IN_FILE_PREFIX}_b0.mif" -pe_dir ap -readout_time 0.0227833 -eddy_mask "${IN_FILE_PREFIX}_denoised_mag_meanb0_bet_mask.nii" -eddy_slspec $slspec -eddyqc_all "$IN_FILE_PATH/eddy_params" -eddy_options " --dfields --flm=cubic --repol --dont_sep_offs_move --fwhm=10,0,0,0,0 --data_is_shelled --mporder=13 --ol_type=both " -topup_options " --nthr=16"
mkdir -p "$IN_FILE_PATH/eddy_params/dfields"
cp "$SCRATCH_DIR/dwi_post_eddy.eddy_displacement_fields"* "$IN_FILE_PATH/eddy_params/dfields"
rm -rf $SCRATCH_DIR

# Convert mrtrix output to nii and bvec/bval
mrconvert -force "${IN_FILE_PREFIX}_moco.mif" -export_grad_fsl "${IN_FILE_PREFIX}_moco_unwarped_bet.bvec" "${IN_FILE_PREFIX}_moco_unwarped_bet.bval" "${IN_FILE_PREFIX}_moco.nii.gz"

# Gradient nonlinearity correction
${SCRIPTPATH}/GradientDistortionUnwarp.sh --workingdir="$IN_FILE_PATH/unwarp_wd" --in="${IN_FILE_PREFIX}_moco" --out="${IN_FILE_PREFIX}_moco_unwarped" --coeffs="${SCRIPTPATH}/../connectom_coeff.grad" --owarp="${IN_FILE_PREFIX}_owarp"

# Calculate and apply new brain mask after eddy and nonlinearity correction
mrconvert -force "${IN_FILE_PREFIX}_moco_unwarped.nii.gz" -fslgrad "${IN_FILE_PREFIX}_moco_unwarped_bet.bvec" "${IN_FILE_PREFIX}_moco_unwarped_bet.bval" "${IN_FILE_PREFIX}_moco_unwarped.mif"
dwiextract -force "${IN_FILE_PREFIX}_moco_unwarped.mif" - -bzero | mrmath -force - mean "${IN_FILE_PREFIX}_moco_unwarped_meanb0.mif" -axis 3
mrconvert -force "${IN_FILE_PREFIX}_moco_unwarped_meanb0.mif" "${IN_FILE_PREFIX}_moco_unwarped_meanb0.nii.gz"
bet "${IN_FILE_PREFIX}_moco_unwarped_meanb0.nii.gz" "${IN_FILE_PREFIX}_moco_unwarped_meanb0_bet.nii.gz" -f 0.4 -m
fslmaths ${IN_FILE_PREFIX}_noise_map.nii.gz -bin -mul "${IN_FILE_PREFIX}_moco_unwarped_meanb0_bet_mask.nii.gz" "${IN_FILE_PREFIX}_moco_unwarped_meanb0_bet_mask.nii.gz" -odt int
fslmaths "${IN_FILE_PREFIX}_moco_unwarped.nii.gz" -mul "${IN_FILE_PREFIX}_moco_unwarped_meanb0_bet_mask.nii.gz" "${IN_FILE_PREFIX}_moco_unwarped_bet.nii.gz"

# Run an initial spherical harmonic decomposition to do an additional registration of the b=30000 shell 
# as sometimes eddy seems to cause a shift between the two shells
amp2sh -force -lmax 6 -shells 0,6000 -rician "${IN_FILE_PREFIX}_noise_map.nii.gz" -fslgrad "${IN_FILE_PREFIX}_moco_unwarped_bet.bvec" "${IN_FILE_PREFIX}_moco_unwarped_bet.bval" "${IN_FILE_PREFIX}_moco_unwarped_bet.nii.gz" "${IN_FILE_PATH}/tmp_b6000.nii.gz"
amp2sh -force -lmax 6 -shells 0,30450 -rician "${IN_FILE_PREFIX}_noise_map.nii.gz" -fslgrad "${IN_FILE_PREFIX}_moco_unwarped_bet.bvec" "${IN_FILE_PREFIX}_moco_unwarped_bet.bval" "${IN_FILE_PREFIX}_moco_unwarped_bet.nii.gz" "${IN_FILE_PATH}/tmp_b30000.nii.gz"
fslsplit "${IN_FILE_PATH}/tmp_b6000.nii.gz" "${IN_FILE_PATH}/tmp_b6000_split"
fslsplit "${IN_FILE_PATH}/tmp_b30000.nii.gz" "${IN_FILE_PATH}/tmp_b30000_split"
fslmaths "${IN_FILE_PATH}/tmp_b6000_split0000.nii.gz" -nan -thr 0 "${IN_FILE_PATH}/tmp_b6000.nii.gz"
fslmaths "${IN_FILE_PATH}/tmp_b30000_split0000.nii.gz" -nan -thr 0 "${IN_FILE_PATH}/tmp_b30000.nii.gz"
flirt -ref "${IN_FILE_PATH}/tmp_b6000.nii.gz" -in "${IN_FILE_PATH}/tmp_b30000.nii.gz" -schedule ${FSLDIR}/etc/flirtsch/ytransonly.sch -out "${IN_FILE_PATH}/tmp_b30000.nii.gz" -omat "${IN_FILE_PATH}/flirt_extra_alignment.mat"
/bin/rm "${IN_FILE_PATH}/tmp_b"*

# Concatenate and apply all warps thus far
python "${SCRIPTPATH}/concat_warps.py" "$IN_FILE_PATH/eddy_params" "$IN_FILE_PATH/unwarp_wd" "${IN_FILE_PREFIX}_moco_unwarped.nii.gz" -f "${IN_FILE_PATH}/flirt_extra_alignment.mat" -b "${IN_FILE_PREFIX}_moco_unwarped_bet.bval"
fslmaths "${IN_FILE_PREFIX}_moco_unwarped.nii.gz" -mul "${IN_FILE_PREFIX}_moco_unwarped_meanb0_bet_mask.nii.gz" "${IN_FILE_PREFIX}_moco_unwarped_bet.nii.gz"

if $mle_flag; then
    # Spherical harmonic decomposition with MLE
    matlab -nodisplay -r "addpath ${SCRIPTPATH}/../AxonRadiusMapping/;fitSH('${IN_FILE_PREFIX}_moco_unwarped_bet.nii.gz', '${IN_FILE_PREFIX}_moco_unwarped_meanb0_bet_mask.nii.gz', '${IN_FILE_PREFIX}_noise_map.nii.gz', '${IN_FILE_PREFIX}_moco_unwarped_bet.bval', '${IN_FILE_PREFIX}_moco_unwarped_bet.bvec', '${IN_FILE_PREFIX}');exit"
    gzip -k -f "${IN_FILE_PREFIX}_sh_b6000_powderavg.nii"
    gzip -k -f "${IN_FILE_PREFIX}_sh_b30000_powderavg.nii"
else
    # Spherical harmonic decomposition
    # Rician bias correction needs up to commit aea92a8 from https://github.com/lukeje/mrtrix3
    amp2sh -force -lmax 6 -shells 0,6000 -rician "${IN_FILE_PREFIX}_noise_map.nii.gz" -fslgrad "${IN_FILE_PREFIX}_moco_unwarped_bet.bvec" "${IN_FILE_PREFIX}_moco_unwarped_bet.bval" "${IN_FILE_PREFIX}_moco_unwarped_bet.nii.gz" "${IN_FILE_PREFIX}_sh_b6000.nii.gz"
    amp2sh -force -lmax 6 -shells 0,30450 -rician "${IN_FILE_PREFIX}_noise_map.nii.gz" -fslgrad "${IN_FILE_PREFIX}_moco_unwarped_bet.bvec" "${IN_FILE_PREFIX}_moco_unwarped_bet.bval" "${IN_FILE_PREFIX}_moco_unwarped_bet.nii.gz" "${IN_FILE_PREFIX}_sh_b30000.nii.gz"
    # Extract 0th order coefficients
    fslsplit "${IN_FILE_PREFIX}_sh_b6000.nii.gz" "${IN_FILE_PREFIX}_sh_b6000_split"
    fslsplit "${IN_FILE_PREFIX}_sh_b30000.nii.gz" "${IN_FILE_PREFIX}_sh_b30000_split"
    /bin/mv "${IN_FILE_PREFIX}_sh_b6000_split0000.nii.gz" "${IN_FILE_PREFIX}_sh_b6000.nii.gz"
    /bin/mv "${IN_FILE_PREFIX}_sh_b30000_split0000.nii.gz" "${IN_FILE_PREFIX}_sh_b30000.nii.gz"
    /bin/rm "${IN_FILE_PATH}/"*"split"*
    # Divide by sqrt(4pi) to get powder average
    fslmaths "${IN_FILE_PREFIX}_sh_b6000.nii.gz" -div 3.5449077018110318 "${IN_FILE_PREFIX}_sh_b6000_powderavg.nii.gz"
    fslmaths "${IN_FILE_PREFIX}_sh_b30000.nii.gz" -div 3.5449077018110318 "${IN_FILE_PREFIX}_sh_b30000_powderavg.nii.gz"
fi

# Normalize to mean b0
fslmaths "${IN_FILE_PREFIX}_sh_b6000_powderavg.nii.gz" -div "${IN_FILE_PREFIX}_moco_unwarped_meanb0.nii.gz" "${IN_FILE_PREFIX}_sh_b6000_powderavg.nii.gz"
fslmaths "${IN_FILE_PREFIX}_sh_b30000_powderavg.nii.gz" -div "${IN_FILE_PREFIX}_moco_unwarped_meanb0.nii.gz" "${IN_FILE_PREFIX}_sh_b30000_powderavg.nii.gz"

# Calculate axon diameters
matlab -nodisplay -r "addpath ${SCRIPTPATH}/../AxonRadiusMapping/;calcAxonMaps('${IN_FILE_PREFIX}_sh_b6000_powderavg.nii.gz', '${IN_FILE_PREFIX}_sh_b30000_powderavg.nii.gz', '${IN_FILE_PREFIX}_moco_unwarped_bet.bval', '${IN_FILE_PREFIX}_moco_unwarped_bet.bvec', '${IN_FILE_PATH}/grad_dev.nii.gz');exit"

# Do quantification if T1 file provided
if test -f ${T1_FILE}; then
    ${SCRIPTPATH}/quantification.sh $IN_FILE $T1_FILE
fi

# Do some DTI fitting to check FA
mkdir -p ${IN_FILE_PATH}/dtifit
dtifit -k "${IN_FILE_PREFIX}_moco_unwarped_bet.nii.gz" -o ${IN_FILE_PATH}/dtifit/FSL_dtifit -m "${IN_FILE_PREFIX}_moco_unwarped_meanb0_bet_mask.nii.gz" -r "${IN_FILE_PREFIX}_moco_unwarped_bet.bvec" -b "${IN_FILE_PREFIX}_moco_unwarped_bet.bval"
