#!/bin/bash

# input
if [ "$#" -lt 1 ]; then
    echo "Diffusion file missing"
    exit 1
elif [ "$#" -lt 2 ]; then
    echo "Axon radius map missing"
    exit 1
elif [ "$#" -lt 3 ]; then
    echo "T1 file missing"
    exit 1
else
    DIF_FILE="$1"
    AXON_FILE="$2"
    T1_FILE="$3"
fi

DIF_FILE_PREFIX=${DIF_FILE%%.*}
AXON_FILE_PREFIX=${AXON_FILE%%.*}
T1_FILE_PREFIX=${T1_FILE%%.*}
DIF_FILE_PATH=$(dirname $DIF_FILE)
SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

# Register mean b0 to MPRAGE using the epi_reg script
mrconvert -force $DIF_FILE -fslgrad "${DIF_FILE_PREFIX}.bvec" "${DIF_FILE_PREFIX}.bval" "${DIF_FILE_PREFIX}.mif"
dwiextract -force "${DIF_FILE_PREFIX}.mif" - -bzero | mrmath -force - mean "${DIF_FILE_PREFIX}_meanb0.mif" -axis 3
mrconvert -force "${DIF_FILE_PREFIX}_meanb0.mif" "${DIF_FILE_PREFIX}_meanb0.nii.gz"
epi_reg --epi="${DIF_FILE_PREFIX}_meanb0.nii.gz" --t1=$T1_FILE --t1brain="${T1_FILE_PREFIX}_bet.nii.gz" --wmseg="${T1_FILE_PREFIX}_bet_seg.nii.gz" --out="${DIF_FILE_PREFIX}_meanb0_reg.nii.gz"

# # Apply registration on axon radius maps
flirt -in $AXON_FILE -ref "${T1_FILE_PREFIX}_bet.nii.gz" -out "${AXON_FILE_PREFIX}_reg.nii.gz" -applyxfm -init "${DIF_FILE_PREFIX}_meanb0_reg.mat"

# # Apply white matter mask on axon radius maps
fslmaths "${AXON_FILE_PREFIX}_reg.nii.gz" -mul "${T1_FILE_PREFIX}_bet_seg.nii.gz" "${AXON_FILE_PREFIX}_wm.nii.gz"
