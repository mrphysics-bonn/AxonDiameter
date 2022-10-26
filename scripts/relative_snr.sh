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
elif [ "$#" -lt 4 ]; then
    echo "No T1 file specified, dont register SNR maps"
    DIF_FILE="$1"
    PA_FILE1="$2"
    PA_FILE2="$3"
else
    DIF_FILE="$1"
    PA_FILE1="$2"
    PA_FILE2="$3"
    T1_FILE="$4"
fi

DIF_FILE_PREFIX=${DIF_FILE%%.*}
PA_FILE1_PREFIX=${PA_FILE1%%.*}
PA_FILE2_PREFIX=${PA_FILE2%%.*}
DIF_FILE_PATH=$(dirname $DIF_FILE)

# Calculate noise map
mrconvert -force $DIF_FILE -fslgrad "${DIF_FILE_PREFIX}.bvec" "${DIF_FILE_PREFIX}.bval" "${DIF_FILE_PREFIX}.mif"
dwiextract -force -shells 0,6000 "${DIF_FILE_PREFIX}.mif" "${DIF_FILE_PREFIX}_lowb.mif"
dwidenoise -force -noise "${DIF_FILE_PREFIX}_noise_map.nii.gz" "${DIF_FILE_PREFIX}_lowb.mif" "${DIF_FILE_PREFIX}_lowb_denoised.mif"

# Mean b0 relative SNR map
dwiextract -force "${DIF_FILE_PREFIX}.mif" - -bzero | mrmath -force - mean "${DIF_FILE_PREFIX}_meanb0.mif" -axis 3
mrconvert -force "${DIF_FILE_PREFIX}_meanb0.mif" "${DIF_FILE_PREFIX}_meanb0.nii.gz"
fslmaths "${DIF_FILE_PREFIX}_meanb0.nii.gz" -div "${DIF_FILE_PREFIX}_noise_map.nii.gz" "${DIF_FILE_PREFIX}_meanb0_relativeSNR.nii.gz"

# Powder-averages relative SNR map (scale by mean b0 as powder-average is also scaled to that)
fslmaths ${PA_FILE1} -div "${DIF_FILE_PREFIX}_noise_map.nii.gz" "${PA_FILE1_PREFIX}_relativeSNR.nii.gz"
fslmaths "${PA_FILE1_PREFIX}_relativeSNR.nii.gz" -mul "${DIF_FILE_PREFIX}_meanb0.nii.gz" "${PA_FILE1_PREFIX}_relativeSNR.nii.gz"
fslmaths ${PA_FILE2} -div "${DIF_FILE_PREFIX}_noise_map.nii.gz" "${PA_FILE2_PREFIX}_relativeSNR.nii.gz"
fslmaths "${PA_FILE2_PREFIX}_relativeSNR.nii.gz" -mul "${DIF_FILE_PREFIX}_meanb0.nii.gz" "${PA_FILE2_PREFIX}_relativeSNR.nii.gz"

if test -f "${T1_FILE_PREFIX}_bet.nii.gz" && test -f "${DIF_FILE_PREFIX}_meanb0_reg.mat"; then
    flirt -in "${DIF_FILE_PREFIX}_meanb0_relativeSNR.nii.gz" -ref "${T1_FILE_PREFIX}_bet.nii.gz" -out "${DIF_FILE_PREFIX}_meanb0_relativeSNR_reg.nii.gz" -applyxfm -init "${DIF_FILE_PREFIX}_meanb0_reg.mat"
    flirt -in "${PA_FILE1_PREFIX}_relativeSNR.nii.gz" -ref "${T1_FILE_PREFIX}_bet.nii.gz" -out "${PA_FILE1_PREFIX}_relativeSNR_reg.nii.gz" -applyxfm -init "${DIF_FILE_PREFIX}_meanb0_reg.mat"
    flirt -in "${PA_FILE2_PREFIX}_relativeSNR.nii.gz" -ref "${T1_FILE_PREFIX}_bet.nii.gz" -out "${PA_FILE2_PREFIX}_relativeSNR_reg.nii.gz" -applyxfm -init "${DIF_FILE_PREFIX}_meanb0_reg.mat"
fi
