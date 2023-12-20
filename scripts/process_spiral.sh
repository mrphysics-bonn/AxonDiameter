#!/bin/bash

usage="$(basename "$0") [-h] [-t T1_FILE] [-m] [-e] DIFF_FILE -- Axon diameter estimation from spiral data

where:
    DIFF_FILE    Spiral diffusion images
    -h           show this help text
    -t T1_FILE   Do white matter masking with T1 data
    -e           Use maximum likelihood estimator instead of denoising & Rician bias correction
    -m           Do denoising on magnitude instead of complex data (only if -e is not selected)"

# option
mag_flag=false
mle_flag=false
while getopts "ht:me" opt; do
    case $opt in
        h) echo "$usage"
        exit
        ;;
        t) T1_FILE=${OPTARG} ;;
        d) mag_flag=true ;;
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
    echo "Input file missing"
    exit 1
else
    IN_FILE="$1"
fi

IN_FILE_PREFIX=${IN_FILE%%.*}
IN_FILE_PATH=$(dirname $IN_FILE)
SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

if $mle_flag; then
    # Convert to magnitude
    python "${SCRIPTPATH}/nifti2mag.py" ${IN_FILE} "${IN_FILE_PREFIX}_mag.nii.gz"
    # Estimate noise for SH decomposition
    dwiextract -force -fslgrad "${IN_FILE_PREFIX}.bvec" "${IN_FILE_PREFIX}.bval" -shells 0,6000 "${IN_FILE_PREFIX}_mag.nii.gz" "${IN_FILE_PREFIX}_mag_lowb.nii.gz"
    dwidenoise -force -noise "${IN_FILE_PREFIX}_noise_map.nii.gz" "${IN_FILE_PREFIX}_mag_lowb.nii.gz"  "${IN_FILE_PREFIX}_mag_lowb_denoised.nii.gz"
    cp "${IN_FILE_PREFIX}_mag.nii.gz" "${IN_FILE_PREFIX}_denoised_mag.nii.gz"
elif $mag_flag; then
    # Convert to magnitude
    python "${SCRIPTPATH}/nifti2mag.py" ${IN_FILE} "${IN_FILE_PREFIX}_mag.nii.gz"
    # Denoising
    dwidenoise -force -noise "${IN_FILE_PREFIX}_noise_map.nii.gz" "${IN_FILE_PREFIX}_mag.nii.gz" "${IN_FILE_PREFIX}_denoised_mag.nii.gz"
else
    # Denoising
    dwidenoise -force -noise "${IN_FILE_PREFIX}_noise_map.nii.gz" ${IN_FILE} "${IN_FILE_PREFIX}_denoised.nii.gz"
    # Convert to magnitude
    python "${SCRIPTPATH}/nifti2mag.py" "${IN_FILE_PREFIX}_denoised.nii.gz" "${IN_FILE_PREFIX}_denoised_mag.nii.gz"
fi

# convert nii to mif
mrconvert -force "${IN_FILE_PREFIX}_denoised_mag.nii.gz" -fslgrad "${IN_FILE_PREFIX}.bvec" "${IN_FILE_PREFIX}.bval" "${IN_FILE_PREFIX}_denoised_mag.mif"

# Gibbs-Ringing removal
mrdegibbs -force "${IN_FILE_PREFIX}_denoised_mag.mif" "${IN_FILE_PREFIX}_denoised_mag_gr.mif"

# Brain mask for eddy
dwiextract -force "${IN_FILE_PREFIX}_denoised_mag_gr.mif" - -bzero | mrmath -force - mean "${IN_FILE_PREFIX}_denoised_mag_meanb0.mif" -axis 3
mrconvert -force "${IN_FILE_PREFIX}_denoised_mag_meanb0.mif" "${IN_FILE_PREFIX}_denoised_mag_meanb0.nii.gz"
bet "${IN_FILE_PREFIX}_denoised_mag_meanb0.nii.gz" "${IN_FILE_PREFIX}_denoised_mag_meanb0_bet.nii.gz" -f 0.4 -m
gunzip -f "${IN_FILE_PREFIX}_denoised_mag_meanb0_bet_mask.nii.gz" # unzip as otherwise the mask gets corrupted by mrconvert in dwifslpreproc

# # this file contains a list of the simultaneously acquired slices in acquisition order
slspec="$SCRIPTPATH/slspec_spiral.txt"

# Motion correction with eddy
# let mrtrix take care of providing eddy input and output
# first level model (flm) is set to movement to avoid eddy current correction, as spiral data has already been corrected
# mporder is recommended to be somewhere between N/4 and N/2, where N is the number of excitations
dwifslpreproc -force "${IN_FILE_PREFIX}_denoised_mag_gr.mif" "${IN_FILE_PREFIX}_moco.mif" -rpe_none -pe_dir ap -eddy_mask "${IN_FILE_PREFIX}_denoised_mag_meanb0_bet_mask.nii" -eddy_slspec $slspec -eddyqc_all "$IN_FILE_PATH/eddy_params" -eddy_options " --flm=movement --repol --data_is_shelled --mporder=13 --ol_type=both "

# Convert mrtrix output to nii and bvec/bval
mrconvert -force "${IN_FILE_PREFIX}_moco.mif" -export_grad_fsl "${IN_FILE_PREFIX}_moco_unwarped_bet.bvec" "${IN_FILE_PREFIX}_moco_unwarped_bet.bval" "${IN_FILE_PREFIX}_moco.nii.gz"

# Gradient nonlinearity correction
${SCRIPTPATH}/GradientDistortionUnwarp.sh --workingdir="$IN_FILE_PATH/unwarp_wd" --in="${IN_FILE_PREFIX}_moco" --out="${IN_FILE_PREFIX}_moco_unwarped" --coeffs="${SCRIPTPATH}/../connectom_coeff.grad" --owarp="${IN_FILE_PREFIX}_owarp"

# Calculate and apply new brain mask after eddy and nonlinearity correction
mrconvert -force "${IN_FILE_PREFIX}_moco_unwarped.nii.gz" -fslgrad "${IN_FILE_PREFIX}_moco_unwarped_bet.bvec" "${IN_FILE_PREFIX}_moco_unwarped_bet.bval" "${IN_FILE_PREFIX}_moco_unwarped.mif"
dwiextract -force "${IN_FILE_PREFIX}_moco_unwarped.mif" - -bzero | mrmath -force - mean "${IN_FILE_PREFIX}_moco_unwarped_meanb0.mif" -axis 3
mrconvert -force "${IN_FILE_PREFIX}_moco_unwarped_meanb0.mif" "${IN_FILE_PREFIX}_moco_unwarped_meanb0.nii.gz"
bet "${IN_FILE_PREFIX}_moco_unwarped_meanb0.nii.gz" "${IN_FILE_PREFIX}_moco_unwarped_meanb0_bet.nii.gz" -f 0.4 -m
fslmaths "${IN_FILE_PREFIX}_moco_unwarped.nii.gz" -mul "${IN_FILE_PREFIX}_moco_unwarped_meanb0_bet_mask.nii.gz" "${IN_FILE_PREFIX}_moco_unwarped_bet.nii.gz"

# Spherical harmonic decomposition
if $mle_flag; then
    matlab -nodisplay -r "addpath ${SCRIPTPATH}/../AxonRadiusMapping/;fitSH('${IN_FILE_PREFIX}_moco_unwarped_bet.nii.gz', '${IN_FILE_PREFIX}_moco_unwarped_meanb0_bet_mask.nii.gz', '${IN_FILE_PREFIX}_noise_map.nii.gz', '${IN_FILE_PREFIX}_moco_unwarped_bet.bval', '${IN_FILE_PREFIX}_moco_unwarped_bet.bvec', '${IN_FILE_PREFIX}');exit"
    gzip -k -f "${IN_FILE_PREFIX}_sh_b6000_powderavg.nii"
    gzip -k -f "${IN_FILE_PREFIX}_sh_b30000_powderavg.nii"
else
    if $mag_flag; then
        amp2sh -force -lmax 6 -shells 0,6000 -normalise -rician "${IN_FILE_PREFIX}_noise_map.nii.gz" -fslgrad "${IN_FILE_PREFIX}_moco_unwarped_bet.bvec" "${IN_FILE_PREFIX}_moco_unwarped_bet.bval" "${IN_FILE_PREFIX}_moco_unwarped_bet.nii.gz" "${IN_FILE_PREFIX}_sh_b6000.nii.gz"
        amp2sh -force -lmax 6 -shells 0,30450 -normalise -rician "${IN_FILE_PREFIX}_noise_map.nii.gz" -fslgrad "${IN_FILE_PREFIX}_moco_unwarped_bet.bvec" "${IN_FILE_PREFIX}_moco_unwarped_bet.bval" "${IN_FILE_PREFIX}_moco_unwarped_bet.nii.gz" "${IN_FILE_PREFIX}_sh_b30000.nii.gz"
    else
        amp2sh -force -lmax 6 -shells 0,6000 -normalise -fslgrad "${IN_FILE_PREFIX}_moco_unwarped_bet.bvec" "${IN_FILE_PREFIX}_moco_unwarped_bet.bval" "${IN_FILE_PREFIX}_moco_unwarped_bet.nii.gz" "${IN_FILE_PREFIX}_sh_b6000.nii.gz"
        amp2sh -force -lmax 6 -shells 0,30450 -normalise -fslgrad "${IN_FILE_PREFIX}_moco_unwarped_bet.bvec" "${IN_FILE_PREFIX}_moco_unwarped_bet.bval" "${IN_FILE_PREFIX}_moco_unwarped_bet.nii.gz" "${IN_FILE_PREFIX}_sh_b30000.nii.gz"
    fi
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

# Calculate axon diameters
matlab -nodisplay -r "addpath ${SCRIPTPATH}/../AxonRadiusMapping/;calcAxonMaps('${IN_FILE_PREFIX}_sh_b6000_powderavg.nii.gz', '${IN_FILE_PREFIX}_sh_b30000_powderavg.nii.gz', '${IN_FILE_PREFIX}_moco_unwarped_bet.bval', '${IN_FILE_PREFIX}_moco_unwarped_bet.bvec', '${IN_FILE_PATH}/grad_dev.nii.gz');exit"

# Do quantification if T1 file provided
if test -f ${T1_FILE}; then
    ${SCRIPTPATH}/quantification.sh $IN_FILE $T1_FILE
fi

# Do some DTI fitting to check FA
mkdir ${IN_FILE_PATH}/dtifit
dtifit -k "${IN_FILE_PREFIX}_moco_unwarped_bet.nii.gz" -o ${IN_FILE_PATH}/dtifit/FSL_dtifit -m "${IN_FILE_PREFIX}_moco_unwarped_meanb0_bet_mask.nii.gz" -r "${IN_FILE_PREFIX}_moco_unwarped_bet.bvec" -b "${IN_FILE_PREFIX}_moco_unwarped_bet.bval"
