#!/bin/bash

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

# Denoising
dwidenoise -noise "${IN_FILE_PREFIX}_noise_map.nii.gz" ${IN_FILE} "${IN_FILE_PREFIX}_denoised.nii.gz"

# Convert to magnitude
python "${SCRIPTPATH}/nifti2mag.py" "${IN_FILE_PREFIX}_denoised.nii.gz" "${IN_FILE_PREFIX}_denoised_mag.nii.gz"

# Motion correction with eddy - convert nii to mif; note must also provide bvecs and bvals
mrconvert "${IN_FILE_PREFIX}_denoised_mag.nii.gz" -fslgrad "${IN_FILE_PREFIX}.bvec" "${IN_FILE_PREFIX}.bval" "${IN_FILE_PREFIX}_denoised_mag.mif"

# this file contains a list of the simultaneously acquired slices in acquisition order
slspec="$SCRIPTPATH/example_slspec.txt"

# let mrtrix take care of providing eddy input and output
# readout_time may need to be longer to satisfy eddy's sanity checking (but not so long that it starts expecting large shifts)
# first level model (flm) is set to linear to avoid overfitting as spiral data should already have been corrected for eddy currents
# mporder is recommended to be somewhere between N/4 and N/2, where N is the number of excitations
dwifslpreproc "${IN_FILE_PREFIX}_denoised_mag.mif" "${IN_FILE_PREFIX}_denoised_mag_moco.mif" -rpe_none -pe_dir ap -readout_time 0.01 -eddy_slspec $slspec -eddyqc_all $IN_FILE_PATH -eddy_options " --flm=linear --repol --data_is_shelled --mporder=13 --ol_type=both "

# Convert mrtrix output to nii and bvec/bval
mrconvert "${IN_FILE_PREFIX}_denoised_mag_moco.mif" -export_grad_fsl "${IN_FILE_PREFIX}_denoised_mag_moco.bvec" "${IN_FILE_PREFIX}_denoised_mag_moco.bval" "${IN_FILE_PREFIX}_denoised_mag_moco.nii.gz"

# Gradient nonlinearity correction
${SCRIPTPATH}/GradientDistortionUnwarp.sh --workingdir="${SCRIPTPATH}/../data/unwarp_wd" --in="${IN_FILE_PREFIX}_denoised_mag_moco" --out="${IN_FILE_PREFIX}_denoised_mag_moco_unwarped" --coeffs="${SCRIPTPATH}/../connectom_coeff.grad" --owarp="${IN_FILE_PREFIX}_owarp"

# Brain masking
fslsplit "${IN_FILE_PREFIX}_denoised_mag_moco_unwarped.nii.gz" "${IN_FILE_PREFIX}_splitted_vol"
bet "${IN_FILE_PREFIX}_splitted_vol0000.nii.gz" "${IN_FILE_PREFIX}_splitted_vol0000.nii.gz" -f 0.3 -m
fslmaths "${IN_FILE_PREFIX}_denoised_mag_moco_unwarped.nii.gz" -mul "${IN_FILE_PREFIX}_splitted_vol0000_mask.nii.gz" "${IN_FILE_PREFIX}_denoised_mag_moco_unwarped.nii.gz"
/bin/rm "${IN_FILE_PREFIX}_splitted_vol"*

# Spherical harmonic decomposition
amp2sh -lmax 6 -shells 0,6000 -normalise -fslgrad "${IN_FILE_PREFIX}_denoised_mag_moco.bvec" "${IN_FILE_PREFIX}_denoised_mag_moco.bval" -rician "${IN_FILE_PREFIX}_noise_map.nii.gz" "${IN_FILE_PREFIX}_denoised_mag_moco_unwarped.nii.gz" "${IN_FILE_PREFIX}_sh_b6000.nii.gz"
amp2sh -lmax 6 -shells 0,30450 -normalise -fslgrad "${IN_FILE_PREFIX}_denoised_mag_moco.bvec" "${IN_FILE_PREFIX}_denoised_mag_moco.bval" -rician "${IN_FILE_PREFIX}_noise_map.nii.gz" "${IN_FILE_PREFIX}_denoised_mag_moco_unwarped.nii.gz" "${IN_FILE_PREFIX}_sh_b30000.nii.gz"

# Divide by sqrt(4pi) to get powder average
fslmaths "${IN_FILE_PREFIX}_sh_b6000.nii.gz" -div 3.5449077018110318 "${IN_FILE_PREFIX}_sh_b6000_powderavg.nii.gz"
fslmaths "${IN_FILE_PREFIX}_sh_b30000.nii.gz" -div 3.5449077018110318 "${IN_FILE_PREFIX}_sh_b30000_powderavg.nii.gz"
