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

if test -f "${T1_FILE_PREFIX}_bet_seg.nii.gz"; then
    echo "White matter mask already exists. Skip calculation."
else
    # Do BET masking on T1 data
    bet $T1_FILE "${T1_FILE_PREFIX}_bet.nii.gz" -f 0.4 -m -B

    # Create white matter mask from T1 data
    fast -g "${T1_FILE_PREFIX}_bet.nii.gz"
fi

# Mask diffusion data
fslsplit $DIF_FILE "${DIF_FILE_PREFIX}_splitted_vol"
bet "${DIF_FILE_PREFIX}_splitted_vol0000.nii.gz" "${DIF_FILE_PREFIX}_splitted_vol0000.nii.gz" -f 0.4 -m
/bin/cp "${DIF_FILE_PREFIX}_splitted_vol0000_mask.nii.gz" "${DIF_FILE_PREFIX}_mask.nii.gz"
/bin/rm "${DIF_FILE_PREFIX}_splitted_vol"*
fslmaths $DIF_FILE -mul "${DIF_FILE_PREFIX}_mask.nii.gz" "${DIF_FILE_PREFIX}_bet.nii.gz"

# Register mean b0 to MPRAGE using the epi_reg script (currently not used)
mrconvert "${DIF_FILE_PREFIX}_bet.nii.gz" -force -fslgrad "${DIF_FILE_PREFIX}.bvec" "${DIF_FILE_PREFIX}.bval" "${DIF_FILE_PREFIX}.mif"
dwiextract -force "${DIF_FILE_PREFIX}.mif" - -bzero | mrmath -force - mean "${DIF_FILE_PREFIX}_meanb0.mif" -axis 3
mrconvert -force "${DIF_FILE_PREFIX}_meanb0.mif" "${DIF_FILE_PREFIX}_meanb0.nii.gz"
epi_reg --epi="${DIF_FILE_PREFIX}_meanb0.nii.gz" --t1=$T1_FILE --t1brain="${T1_FILE_PREFIX}_bet.nii.gz" --wmseg="${T1_FILE_PREFIX}_bet_seg_2.nii.gz" --out="${DIF_FILE_PREFIX}_meanb0_reg.nii.gz"

# Or register FA to MPRAGE (6 parameter model) - currently not used
# calc_FA -i "${DIF_FILE_PREFIX}_bet.nii.gz" -o "${DIF_FILE_PATH}/FA.nii.gz" --bvals "${DIF_FILE_PREFIX}.bval" --bvecs "${DIF_FILE_PREFIX}.bvec" --brain_mask "${DIF_FILE_PREFIX}_mask.nii.gz"
# flirt -ref "${T1_FILE_PREFIX}_bet.nii.gz" -in "${DIF_FILE_PATH}/FA.nii.gz" -out "${DIF_FILE_PATH}/FA.nii.gz" -omat "${DIF_FILE_PREFIX}_meanb0_reg.mat" -dof 6 -cost mutualinfo -searchcost mutualinfo
# flirt -ref "${T1_FILE_PREFIX}_bet.nii.gz" -in "${DIF_FILE_PREFIX}_meanb0.nii.gz" -out "${DIF_FILE_PREFIX}_meanb0_reg.nii.gz" -applyxfm -init "${DIF_FILE_PREFIX}_meanb0_reg.mat" -dof 6

# # Apply registration on axon radius maps
flirt -in $AXON_FILE -ref "${T1_FILE_PREFIX}_bet.nii.gz" -out "${AXON_FILE_PREFIX}_reg.nii.gz" -applyxfm -init "${DIF_FILE_PREFIX}_meanb0_reg.mat"

# # Apply white matter mask on axon radius maps
fslmaths "${AXON_FILE_PREFIX}_reg.nii.gz" -mul "${T1_FILE_PREFIX}_bet_seg_2.nii.gz" "${AXON_FILE_PREFIX}_wm.nii.gz"
