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

# MoCo - not really working well
# mcflirt -in "${IN_FILE_pf}_denoised_mag.nii.gz" -o "${IN_FILE_pf}_denoised_mag_moco.nii.gz" -cost corratio -refvol 0 -spline_final -plots
# python "${SCRIPTPATH}/rot_bvec.py" "${IN_FILE_pf}.bvec" "${IN_FILE_pf}_denoised_mag_moco.nii.gz.par"