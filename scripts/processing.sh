#!/bin/bash

# input
if [ "$#" -lt 1 ]; then
    echo "Input file missing"
    exit 1
else
    IN_FILE="$1"
fi

IN_FILE_pf=${IN_FILE%%.*}
SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

# Denoising
dwidenoise -noise "${IN_FILE_pf}_noise_map.nii.gz" ${IN_FILE} "${IN_FILE_pf}_denoised.nii.gz"

# Convert to magnitude
python "${SCRIPTPATH}/nifti2mag.py" "${IN_FILE_pf}_denoised.nii.gz" "${IN_FILE_pf}_denoised_mag.nii.gz"

# Motion Correction - results are pretty awful, so currently not used 
# mcflirt -in "${IN_FILE_pf}_denoised_mag.nii.gz" -o "${IN_FILE_pf}_denoised_mag_moco.nii.gz" -cost corratio -refvol 0 -spline_final -plots
# python "${SCRIPTPATH}/rot_bvec.py" "${IN_FILE_pf}.bvec" "${IN_FILE_pf}_denoised_mag_moco.nii.gz.par" # b-value rotation

# Gradient nonlinearity correction
${SCRIPTPATH}/GradientDistortionUnwarp.sh --workingdir="${SCRIPTPATH}/../data/unwarp_wd" --in="${IN_FILE_pf}_denoised_mag" --out="${IN_FILE_pf}_unwarped" --coeffs="${SCRIPTPATH}/../connectom_coeff.grad" --owarp="${IN_FILE_pf}_owarp"