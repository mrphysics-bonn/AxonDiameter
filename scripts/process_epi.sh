#!/bin/bash

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

# Denoising
dwidenoise -noise "${IN_FILE_PREFIX}_noise_map.nii.gz" ${IN_FILE} "${IN_FILE_PREFIX}_denoised_mag.nii.gz"

# Convert nii to mif
mrconvert "${IN_FILE_PREFIX}_denoised_mag.nii.gz" -fslgrad "${IN_FILE_PREFIX}.bvec" "${IN_FILE_PREFIX}.bval" "${IN_FILE_PREFIX}_denoised_mag.mif"

# Low-SNR Rician bias correction (from DESIGNER pipeline)
mrconvert ${IN_FILE} -fslgrad "${IN_FILE_PREFIX}.bvec" "${IN_FILE_PREFIX}.bval" "${IN_FILE_PREFIX}_ricinput.mif"
dwiextract -shell 30450 "${IN_FILE_PREFIX}_ricinput.mif" dwilowb.mif
dwidenoise -noise lowbnoisemap.mif dwilowb.mif dwitmp.mif
mrcalc "${IN_FILE_PREFIX}_denoised_mag.mif" 2 -pow lowbnoisemap.mif 2 -pow -sub -abs -sqrt - | mrcalc - -finite - 0 -if "${IN_FILE_PREFIX}_denoised_mag_rician.mif"
/bin/rm dwilowb.mif lowbnoisemap.mif dwitmp.mif dwirc.mif "${IN_FILE_PREFIX}_ricinput.mif"

# Gibbs-Ringing removal
mrdegibbs "${IN_FILE_PREFIX}_denoised_mag_rician.mif" "${IN_FILE_PREFIX}_denoised_mag_gr.mif"

# this file contains a list of the simultaneously acquired slices in acquisition order
slspec="$SCRIPTPATH/example_slspec.txt"

# input for topup
dwiextract "${IN_FILE}" - -bzero -fslgrad "${IN_FILE_PREFIX}.bvec" "${IN_FILE_PREFIX}.bval" | mrmath - mean "${IN_FILE_PREFIX}_AP_b0.mif" -axis 3
mrmath "${PA_FILE}" mean "${IN_FILE_PREFIX}_PA_b0.mif" -axis 3
mrcat "${IN_FILE_PREFIX}_AP_b0.mif" "${IN_FILE_PREFIX}_PA_b0.mif" "${IN_FILE_PREFIX}_b0.mif" -axis 3

# let mrtrix take care of providing eddy and topup input and output
# readout_time is copied from the json file for our data
# mporder is recommended to be somewhere between N/4 and N/2, where N is the number of excitations
dwifslpreproc "${IN_FILE_PREFIX}_denoised_mag_gr.mif" "${IN_FILE_PREFIX}_moco.mif" -rpe_pair -se_epi "${IN_FILE_PREFIX}_b0.mif" -pe_dir ap -readout_time 0.0227833 -eddy_slspec $slspec -eddyqc_all $IN_FILE_PATH -eddy_options " --flm=cubic --slm=linear --repol --data_is_shelled --mporder=13 --ol_type=both "

# Convert mrtrix output to nii and bvec/bval
mrconvert "${IN_FILE_PREFIX}_moco.mif" -export_grad_fsl "${IN_FILE_PREFIX}_moco.bvec" "${IN_FILE_PREFIX}_moco.bval" "${IN_FILE_PREFIX}_moco.nii.gz"

# Gradient nonlinearity correction
${SCRIPTPATH}/GradientDistortionUnwarp.sh --workingdir="${SCRIPTPATH}/../data/unwarp_wd" --in="${IN_FILE_PREFIX}_moco" --out="${IN_FILE_PREFIX}_moco_unwarped" --coeffs="${SCRIPTPATH}/../connectom_coeff.grad" --owarp="${IN_FILE_PREFIX}_owarp"

# Brain masking
fslsplit "${IN_FILE_PREFIX}_moco_unwarped.nii.gz" "${IN_FILE_PREFIX}_splitted_vol"
bet "${IN_FILE_PREFIX}_splitted_vol0000.nii.gz" "${IN_FILE_PREFIX}_splitted_vol0000.nii.gz" -f 0.3 -m
fslmaths "${IN_FILE_PREFIX}_moco_unwarped.nii.gz" -mul "${IN_FILE_PREFIX}_splitted_vol0000_mask.nii.gz" "${IN_FILE_PREFIX}_moco_unwarped.nii.gz"
/bin/rm "${IN_FILE_PREFIX}_splitted_vol"*

# Spherical harmonic decomposition - disable Rician noise correction as its currently not correct
amp2sh -lmax 6 -shells 0,6000 -normalise -fslgrad "${IN_FILE_PREFIX}_moco.bvec" "${IN_FILE_PREFIX}_moco.bval" "${IN_FILE_PREFIX}_moco_unwarped.nii.gz" "${IN_FILE_PREFIX}_sh_b6000.nii.gz"
amp2sh -lmax 6 -shells 0,30450 -normalise -fslgrad "${IN_FILE_PREFIX}_moco.bvec" "${IN_FILE_PREFIX}_moco.bval" "${IN_FILE_PREFIX}_moco_unwarped.nii.gz" "${IN_FILE_PREFIX}_sh_b30000.nii.gz"

# Divide by sqrt(4pi) to get powder average
fslmaths "${IN_FILE_PREFIX}_sh_b6000.nii.gz" -div 3.5449077018110318 "${IN_FILE_PREFIX}_sh_b6000_powderavg.nii.gz"
fslmaths "${IN_FILE_PREFIX}_sh_b30000.nii.gz" -div 3.5449077018110318 "${IN_FILE_PREFIX}_sh_b30000_powderavg.nii.gz"

# Calculate axon diameters
matlab -nodisplay -r "addpath ${SCRIPTPATH}/../AxonRadiusMapping/;calcAxonMaps('${IN_FILE_PREFIX}_sh_b6000_powderavg.nii.gz', '${IN_FILE_PREFIX}_sh_b30000_powderavg.nii.gz', '${IN_FILE_PREFIX}.bval', '${IN_FILE_PREFIX}.bvec', '${IN_FILE_PATH}/grad_dev.nii.gz');exit"
