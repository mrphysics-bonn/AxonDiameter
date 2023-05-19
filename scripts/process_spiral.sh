#!/bin/bash

# option
mag_flag=false
while getopts "m" mag; do
    mag_flag=true
    echo "Denoising after conversion to magnitude images."
done
shift "$((OPTIND-1))"

# input
if [ "$#" -lt 1 ]; then
    echo "Input file missing"
    exit 1
elif [ "$#" -lt 2 ]; then
    echo "No T1 file specified, dont do white matter masking."
    IN_FILE="$1"
else
    IN_FILE="$1"
    T1_FILE="$2"
fi

IN_FILE_PREFIX=${IN_FILE%%.*}
IN_FILE_PATH=$(dirname $IN_FILE)
SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

if $mag_flag; then
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
slspec="$SCRIPTPATH/example_slspec.txt"

# Motion correction with eddy
# let mrtrix take care of providing eddy input and output
# first level model (flm) is set to movement to avoid eddy current correction, as spiral data has already been corrected
# mporder is recommended to be somewhere between N/4 and N/2, where N is the number of excitations
dwifslpreproc -force "${IN_FILE_PREFIX}_denoised_mag_gr.mif" "${IN_FILE_PREFIX}_moco.mif" -rpe_none -pe_dir ap -eddy_mask "${IN_FILE_PREFIX}_denoised_mag_meanb0_bet_mask.nii" -eddy_slspec $slspec -eddyqc_all "$IN_FILE_PATH/eddy_params" -eddy_options " --flm=movement --repol --data_is_shelled --mporder=13 --ol_type=both "

# Convert mrtrix output to nii and bvec/bval
mrconvert -force "${IN_FILE_PREFIX}_moco.mif" -export_grad_fsl "${IN_FILE_PREFIX}_moco_unwarped_bet.bvec" "${IN_FILE_PREFIX}_moco_unwarped_bet.bval" "${IN_FILE_PREFIX}_moco.nii.gz"

# Gradient nonlinearity correction
${SCRIPTPATH}/GradientDistortionUnwarp.sh --workingdir="$IN_FILE_PATH/unwarp_wd" --in="${IN_FILE_PREFIX}_moco" --out="${IN_FILE_PREFIX}_moco_unwarped" --coeffs="${SCRIPTPATH}/../connectom_coeff.grad" --owarp="${IN_FILE_PREFIX}_owarp"

# Spherical harmonic decomposition
# Use Rician bias correction only for magnitude data - needs up to commit aea92a8 from https://github.com/lukeje/mrtrix3
# Calculate one dataset without normalization for relative noise estimation
if $mag_flag; then
    amp2sh -force -lmax 6 -shells 0,6000 -normalise -rician "${IN_FILE_PREFIX}_noise_map.nii.gz" -fslgrad "${IN_FILE_PREFIX}_moco_unwarped_bet.bvec" "${IN_FILE_PREFIX}_moco_unwarped_bet.bval" "${IN_FILE_PREFIX}_moco_unwarped.nii.gz" "${IN_FILE_PREFIX}_sh_b6000.nii.gz"
    amp2sh -force -lmax 6 -shells 0,30450 -normalise -rician "${IN_FILE_PREFIX}_noise_map.nii.gz" -fslgrad "${IN_FILE_PREFIX}_moco_unwarped_bet.bvec" "${IN_FILE_PREFIX}_moco_unwarped_bet.bval" "${IN_FILE_PREFIX}_moco_unwarped.nii.gz" "${IN_FILE_PREFIX}_sh_b30000.nii.gz"
else
    amp2sh -force -lmax 6 -shells 0,6000 -normalise -fslgrad "${IN_FILE_PREFIX}_moco_unwarped_bet.bvec" "${IN_FILE_PREFIX}_moco_unwarped_bet.bval" "${IN_FILE_PREFIX}_moco_unwarped.nii.gz" "${IN_FILE_PREFIX}_sh_b6000.nii.gz"
    amp2sh -force -lmax 6 -shells 0,30450 -normalise -fslgrad "${IN_FILE_PREFIX}_moco_unwarped_bet.bvec" "${IN_FILE_PREFIX}_moco_unwarped_bet.bval" "${IN_FILE_PREFIX}_moco_unwarped.nii.gz" "${IN_FILE_PREFIX}_sh_b30000.nii.gz"
fi

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
${SCRIPTPATH}/tractography.sh "${IN_FILE_PREFIX}_moco_unwarped_bet.nii.gz"
