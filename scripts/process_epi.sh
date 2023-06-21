#!/bin/bash

# input
if [ "$#" -lt 1 ]; then
    echo "AP file missing"
    exit 1
elif [ "$#" -lt 2 ]; then
    echo "PA file missing"
    exit 1
elif [ "$#" -lt 3 ]; then
    echo "No T1 file specified, dont do white matter masking."
    IN_FILE="$1"
    PA_FILE="$2"
else
    IN_FILE="$1"
    PA_FILE="$2"
    T1_FILE="$3"
fi

IN_FILE_PREFIX=${IN_FILE%%.*}
IN_FILE_PATH=$(dirname $IN_FILE)
SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

# Denoising
dwidenoise -force -noise "${IN_FILE_PREFIX}_noise_map.nii.gz" ${IN_FILE} "${IN_FILE_PREFIX}_denoised_mag.nii.gz"

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
slspec="$SCRIPTPATH/example_slspec.txt"

# input for topup
dwiextract -force "${IN_FILE}" - -bzero -fslgrad "${IN_FILE_PREFIX}.bvec" "${IN_FILE_PREFIX}.bval" | mrmath -force - mean "${IN_FILE_PREFIX}_AP_b0.mif" -axis 3
mrmath -force "${PA_FILE}" mean "${IN_FILE_PREFIX}_PA_b0.mif" -axis 3
mrcat -force "${IN_FILE_PREFIX}_AP_b0.mif" "${IN_FILE_PREFIX}_PA_b0.mif" "${IN_FILE_PREFIX}_b0.mif" -axis 3

# let mrtrix take care of providing eddy and topup input and output
# readout_time is copied from the json file for our data
# mporder is recommended to be somewhere between N/4 and N/2, where N is the number of excitations
dwifslpreproc -force "${IN_FILE_PREFIX}_denoised_mag_gr.mif" "${IN_FILE_PREFIX}_moco.mif" -rpe_pair -se_epi "${IN_FILE_PREFIX}_b0.mif" -pe_dir ap -readout_time 0.0227833 -eddy_mask "${IN_FILE_PREFIX}_denoised_mag_meanb0_bet_mask.nii" -eddy_slspec $slspec -eddyqc_all "$IN_FILE_PATH/eddy_params" -eddy_options " --flm=cubic --repol --dont_sep_offs_move --fwhm=10,0,0,0,0 --data_is_shelled --mporder=13 --ol_type=both "

# Convert mrtrix output to nii and bvec/bval
mrconvert -force "${IN_FILE_PREFIX}_moco.mif" -export_grad_fsl "${IN_FILE_PREFIX}_moco_unwarped_bet.bvec" "${IN_FILE_PREFIX}_moco_unwarped_bet.bval" "${IN_FILE_PREFIX}_moco.nii.gz"

# Gradient nonlinearity correction
${SCRIPTPATH}/GradientDistortionUnwarp.sh --workingdir="$IN_FILE_PATH/unwarp_wd" --in="${IN_FILE_PREFIX}_moco" --out="${IN_FILE_PREFIX}_moco_unwarped" --coeffs="${SCRIPTPATH}/../connectom_coeff.grad" --owarp="${IN_FILE_PREFIX}_owarp"

# Spherical harmonic decomposition
# Rician bias correction needs up to commit aea92a8 from https://github.com/lukeje/mrtrix3
# Calculate one dataset without normalization for relative noise estimation
amp2sh -force -lmax 6 -shells 0,6000 -normalise -rician "${IN_FILE_PREFIX}_noise_map.nii.gz" -fslgrad "${IN_FILE_PREFIX}_moco_unwarped_bet.bvec" "${IN_FILE_PREFIX}_moco_unwarped_bet.bval" "${IN_FILE_PREFIX}_moco_unwarped.nii.gz" "${IN_FILE_PREFIX}_sh_b6000.nii.gz"
amp2sh -force -lmax 6 -shells 0,30450 -normalise -rician "${IN_FILE_PREFIX}_noise_map.nii.gz" -fslgrad "${IN_FILE_PREFIX}_moco_unwarped_bet.bvec" "${IN_FILE_PREFIX}_moco_unwarped_bet.bval" "${IN_FILE_PREFIX}_moco_unwarped.nii.gz" "${IN_FILE_PREFIX}_sh_b30000.nii.gz"

# Register decomposed shells as sometimes eddy causes a shift between the two shells - remove nans and negatives first
fslmaths "${IN_FILE_PREFIX}_sh_b6000.nii.gz" -nan -thr 0 -uthr 1.5 "${IN_FILE_PREFIX}_sh_b6000.nii.gz"
fslmaths "${IN_FILE_PREFIX}_sh_b30000.nii.gz" -nan -thr 0 -uthr 1.5 "${IN_FILE_PREFIX}_sh_b30000.nii.gz"
flirt -ref "${IN_FILE_PREFIX}_sh_b6000.nii.gz" -in "${IN_FILE_PREFIX}_sh_b30000.nii.gz" -schedule ${FSLDIR}/etc/flirtsch/ytransonly.sch -out "${IN_FILE_PREFIX}_sh_b30000.nii.gz"

# Brain masking with mean b0
mrconvert -force "${IN_FILE_PREFIX}_moco_unwarped.nii.gz" -fslgrad "${IN_FILE_PREFIX}_moco_unwarped_bet.bvec" "${IN_FILE_PREFIX}_moco_unwarped_bet.bval" "${IN_FILE_PREFIX}_moco_unwarped.mif"
dwiextract -force "${IN_FILE_PREFIX}_moco_unwarped.mif" - -bzero | mrmath -force - mean "${IN_FILE_PREFIX}_moco_unwarped_meanb0.mif" -axis 3
mrconvert -force "${IN_FILE_PREFIX}_moco_unwarped_meanb0.mif" "${IN_FILE_PREFIX}_moco_unwarped_meanb0.nii.gz"
bet "${IN_FILE_PREFIX}_moco_unwarped_meanb0.nii.gz" "${IN_FILE_PREFIX}_moco_unwarped_meanb0_bet.nii.gz" -f 0.3 -m

# Apply brain mask
fslmaths "${IN_FILE_PREFIX}_moco_unwarped.nii.gz" -mul "${IN_FILE_PREFIX}_moco_unwarped_meanb0_bet_mask.nii.gz" "${IN_FILE_PREFIX}_moco_unwarped_bet.nii.gz"
fslmaths "${IN_FILE_PREFIX}_sh_b6000.nii.gz" -mul "${IN_FILE_PREFIX}_moco_unwarped_meanb0_bet_mask.nii.gz" "${IN_FILE_PREFIX}_sh_b6000.nii.gz"
fslmaths "${IN_FILE_PREFIX}_sh_b30000.nii.gz" -mul "${IN_FILE_PREFIX}_moco_unwarped_meanb0_bet_mask.nii.gz" "${IN_FILE_PREFIX}_sh_b30000.nii.gz"

# Extract 0th order coefficients
fslsplit "${IN_FILE_PREFIX}_sh_b6000.nii.gz" "${IN_FILE_PREFIX}_sh_b6000_split"
fslsplit "${IN_FILE_PREFIX}_sh_b30000.nii.gz" "${IN_FILE_PREFIX}_sh_b30000_split"
/bin/mv "${IN_FILE_PREFIX}_sh_b6000_split0000.nii.gz" "${IN_FILE_PREFIX}_sh_b6000.nii.gz"
/bin/mv "${IN_FILE_PREFIX}_sh_b30000_split0000.nii.gz" "${IN_FILE_PREFIX}_sh_b30000.nii.gz"
/bin/rm "${IN_FILE_PATH}/"*"split"*

# Divide by sqrt(4pi) to get powder average
fslmaths "${IN_FILE_PREFIX}_sh_b6000.nii.gz" -div 3.5449077018110318 "${IN_FILE_PREFIX}_sh_b6000_powderavg.nii.gz"
fslmaths "${IN_FILE_PREFIX}_sh_b30000.nii.gz" -div 3.5449077018110318 "${IN_FILE_PREFIX}_sh_b30000_powderavg.nii.gz"

# Calculate axon diameters
matlab -nodisplay -r "addpath ${SCRIPTPATH}/../AxonRadiusMapping/;calcAxonMaps('${IN_FILE_PREFIX}_sh_b6000_powderavg.nii.gz', '${IN_FILE_PREFIX}_sh_b30000_powderavg.nii.gz', '${IN_FILE_PREFIX}_moco_unwarped_bet.bval', '${IN_FILE_PREFIX}_moco_unwarped_bet.bvec', '${IN_FILE_PATH}/grad_dev.nii.gz');exit"

# Calculate relative SNR maps, white matter masks & do tractography
if test -f ${T1_FILE}; then
    ${SCRIPTPATH}/wm_axons.sh "${IN_FILE_PREFIX}_moco_unwarped_bet.nii.gz" "${IN_FILE_PATH}/AxonRadiusMap.nii" ${T1_FILE}
fi
${SCRIPTPATH}/relative_snr.sh "${IN_FILE_PREFIX}_moco_unwarped_bet.nii.gz" "${IN_FILE_PREFIX}_sh_b6000_powderavg.nii.gz" "${IN_FILE_PREFIX}_sh_b30000_powderavg.nii.gz" "${IN_FILE_PREFIX}_noise_map.nii.gz" ${T1_FILE}
# ${SCRIPTPATH}/tractography.sh "${IN_FILE_PREFIX}_moco_unwarped_bet.nii.gz"

# improved tractography with MrTrix & along-fibre quantification (only CST left atm)
${SCRIPTPATH}/along_tract/along_tract_CST.sh "${IN_FILE_PREFIX}_moco_unwarped_bet.nii.gz" "${T1_FILE%%.*}_bet.nii.gz" "${IN_FILE_PREFIX}_sh_b6000_powderavg.nii.gz" "${IN_FILE_PREFIX}_sh_b30000_powderavg.nii.gz" "${IN_FILE_PREFIX}_moco_unwarped_meanb0_bet_mask.nii.gz"

# recalculate axon diameters with along fibre quantified powder-averaged shells
PATH_TRACT=${IN_FILE_PATH}/along_tract/CST_L
matlab -nodisplay -r "addpath ${SCRIPTPATH}/../AxonRadiusMapping/;calcAxonAlongTracts('$PATH_TRACT/CST_stats_sh6000.txt', '$PATH_TRACT/CST_stats_sh30000.txt', '$PATH_TRACT/CST_grads_6000.txt', '$PATH_TRACT/CST_grads_30000.txt', 'CST_stats_Axon_afterAFQonShells.txt');exit"
matlab -nodisplay -r "addpath ${SCRIPTPATH}/../AxonRadiusMapping/;calcAxonAlongTracts('$PATH_TRACT/CST_stats_between_regions_sh6000.txt', '$PATH_TRACT/CST_stats_between_regions_sh30000.txt', '$PATH_TRACT/CST_grads_between_regions_6000.txt', '$PATH_TRACT/CST_grads_between_regions_30000.txt', 'CST_stats_between_regions_Axon_afterAFQonShells.txt');exit"
