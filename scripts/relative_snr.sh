#!/bin/bash

# input
if [ "$#" -lt 1 ]; then
    echo "Diffusion file missing"
    exit 1
elif [ "$#" -lt 2 ]; then
    echo "Powder average shell 1 missing"
    exit 1
elif [ "$#" -lt 3 ]; then
    echo "Powder average shell 2 missing"
    exit 1
else
    DIF_FILE="$1"
    PA_FILE1="$2"
    PA_FILE2="$3"
fi

DIF_FILE_PREFIX=${DIF_FILE%%.*}
PA_FILE1_PREFIX=${PA_FILE1%%.*}
PA_FILE2_PREFIX=${PA_FILE2%%.*}
DIF_FILE_PATH=$(dirname $DIF_FILE)

# Calculate noise map
dwidenoise -noise "${DIF_FILE_PREFIX}_noise_map.nii.gz" ${DIF_FILE} "${DIF_FILE_PREFIX}_denoised.nii.gz"

# Mean b0 relative SNR map
mrconvert $DIF_FILE -force -fslgrad "${DIF_FILE_PREFIX}.bvec" "${DIF_FILE_PREFIX}.bval" "${DIF_FILE_PREFIX}.mif"
dwiextract -force "${DIF_FILE_PREFIX}.mif" - -bzero | mrmath -force - mean "${DIF_FILE_PREFIX}_meanb0.mif" -axis 3
mrconvert -force "${DIF_FILE_PREFIX}_meanb0.mif" "${DIF_FILE_PREFIX}_meanb0.nii.gz"
fslmaths "${DIF_FILE_PREFIX}_meanb0.nii.gz" -div "${DIF_FILE_PREFIX}_noise_map.nii.gz" "${DIF_FILE_PREFIX}_meanb0_relativeSNR.nii.gz"

# Powder-averages relative SNR map
fslmaths ${PA_FILE1} -div "${DIF_FILE_PREFIX}_noise_map.nii.gz" "${PA_FILE1_PREFIX}_relativeSNR.nii.gz"
fslmaths ${PA_FILE2} -div "${DIF_FILE_PREFIX}_noise_map.nii.gz" "${PA_FILE2_PREFIX}_relativeSNR.nii.gz"