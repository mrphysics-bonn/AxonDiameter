#!/bin/bash

# input
if [ "$#" -lt 1 ]; then
    echo "Diffusion file missing"
    exit 1
elif [ "$#" -lt 2 ]; then
    echo "SH shell 1 missing"
    exit 1
elif [ "$#" -lt 3 ]; then
    echo "SH shell 2 missing"
    exit 1
elif [ "$#" -lt 4 ]; then
    echo "Noise map missing"
    exit 1
elif [ "$#" -lt 5 ]; then
    echo "No T1 file specified, dont register SNR maps"
    DIF_FILE="$1"
    SH_FILE1="$2"
    SH_FILE2="$3"
    NOISE_FILE="$4"
else
    DIF_FILE="$1"
    SH_FILE1="$2"
    SH_FILE2="$3"
    NOISE_FILE="$4"
    T1_FILE="$5"
fi

DIF_FILE_PREFIX=${DIF_FILE%%.*}
SH_FILE1_PREFIX=${SH_FILE1%%.*}
SH_FILE2_PREFIX=${SH_FILE2%%.*}
T1_FILE_PREFIX=${T1_FILE%%.*}

# Mean b0 relative SNR map
mrconvert -force $DIF_FILE -fslgrad "${DIF_FILE_PREFIX}.bvec" "${DIF_FILE_PREFIX}.bval" "${DIF_FILE_PREFIX}.mif"
dwiextract -force "${DIF_FILE_PREFIX}.mif" - -bzero | mrmath -force - mean "${DIF_FILE_PREFIX}_meanb0.mif" -axis 3
mrconvert -force "${DIF_FILE_PREFIX}_meanb0.mif" "${DIF_FILE_PREFIX}_meanb0.nii.gz"
fslmaths "${DIF_FILE_PREFIX}_meanb0.nii.gz" -div $NOISE_FILE "${DIF_FILE_PREFIX}_meanb0_relativeSNR.nii.gz"

# Spherical-harmonic decomposition relative SNR map
fslmaths ${SH_FILE1} -div $NOISE_FILE "${SH_FILE1_PREFIX}_relativeSNR.nii.gz"
fslmaths ${SH_FILE2} -div $NOISE_FILE "${SH_FILE2_PREFIX}_relativeSNR.nii.gz"
# Scale with mean b0 as SH were normalized with that
fslmaths "${SH_FILE1_PREFIX}_relativeSNR.nii.gz" -mul "${DIF_FILE_PREFIX}_meanb0.nii.gz" "${SH_FILE1_PREFIX}_relativeSNR.nii.gz"
fslmaths "${SH_FILE2_PREFIX}_relativeSNR.nii.gz" -mul "${DIF_FILE_PREFIX}_meanb0.nii.gz" "${SH_FILE2_PREFIX}_relativeSNR.nii.gz"

if test -f $T1_FILE && test -f "${DIF_FILE_PREFIX}_meanb0_reg.mat"; then
    # Register
    flirt -in "${DIF_FILE_PREFIX}_meanb0_relativeSNR.nii.gz" -ref $T1_FILE -out "${DIF_FILE_PREFIX}_meanb0_relativeSNR_reg.nii.gz" -applyxfm -init "${DIF_FILE_PREFIX}_meanb0_reg.mat"
    flirt -in "${SH_FILE1_PREFIX}_relativeSNR.nii.gz" -ref $T1_FILE -out "${SH_FILE1_PREFIX}_relativeSNR_reg.nii.gz" -applyxfm -init "${DIF_FILE_PREFIX}_meanb0_reg.mat"
    flirt -in "${SH_FILE2_PREFIX}_relativeSNR.nii.gz" -ref $T1_FILE -out "${SH_FILE2_PREFIX}_relativeSNR_reg.nii.gz" -applyxfm -init "${DIF_FILE_PREFIX}_meanb0_reg.mat"

    # Apply mask
    fslmaths "${DIF_FILE_PREFIX}_meanb0_relativeSNR_reg.nii.gz" -mul "${T1_FILE_PREFIX}_mask.nii.gz" "${DIF_FILE_PREFIX}_meanb0_relativeSNR_reg.nii.gz"
    fslmaths "${SH_FILE1_PREFIX}_relativeSNR_reg.nii.gz" -mul "${T1_FILE_PREFIX}_mask.nii.gz" "${SH_FILE1_PREFIX}_relativeSNR_reg.nii.gz"
    fslmaths "${SH_FILE2_PREFIX}_relativeSNR_reg.nii.gz" -mul "${T1_FILE_PREFIX}_mask.nii.gz" "${SH_FILE2_PREFIX}_relativeSNR_reg.nii.gz"

    # Mask only white matter
    fslmaths "${DIF_FILE_PREFIX}_meanb0_relativeSNR_reg.nii.gz" -mul "${T1_FILE_PREFIX}_seg.nii.gz" "${DIF_FILE_PREFIX}_meanb0_relativeSNR_reg_wm.nii.gz"
    fslmaths "${SH_FILE1_PREFIX}_relativeSNR_reg.nii.gz" -mul "${T1_FILE_PREFIX}_seg.nii.gz" "${SH_FILE1_PREFIX}_relativeSNR_reg_wm.nii.gz"
    fslmaths "${SH_FILE2_PREFIX}_relativeSNR_reg.nii.gz" -mul "${T1_FILE_PREFIX}_seg.nii.gz" "${SH_FILE2_PREFIX}_relativeSNR_reg_wm.nii.gz"
fi
