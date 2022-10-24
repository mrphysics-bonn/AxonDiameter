#!/bin/bash

# input
if [ "$#" -lt 1 ]; then
    echo "Input diffusion file missing"
    exit 1
else
    IN_FILE="$1"
fi

IN_FILE_PREFIX=${IN_FILE%%.*}
IN_FILE_PATH=$(dirname $IN_FILE)
SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

# Brain masking
fslsplit $IN_FILE "${IN_FILE_PREFIX}_splitted_vol"
bet "${IN_FILE_PREFIX}_splitted_vol0000.nii.gz" "${IN_FILE_PREFIX}_splitted_vol0000.nii.gz" -f 0.3 -m
/bin/cp "${IN_FILE_PREFIX}_splitted_vol0000_mask.nii.gz" "${IN_FILE_PREFIX}_mask.nii.gz"
/bin/rm "${IN_FILE_PREFIX}_splitted_vol"*

# Register diffusion images to MNI space
calc_FA -i $IN_FILE -o "${IN_FILE_PATH}/FA.nii.gz" --bvals "${IN_FILE_PREFIX}.bval" --bvecs "${IN_FILE_PREFIX}.bvec"  --brain_mask "${IN_FILE_PREFIX}_mask.nii.gz"
flirt -ref "${SCRIPTPATH}/../TractSeg/tractseg/resources/MNI_FA_template.nii.gz" -in "${IN_FILE_PATH}/FA.nii.gz" -out "${IN_FILE_PATH}/FA_MNI.nii.gz" -omat "${IN_FILE_PATH}/FA_MNI.mat" -dof 6 -cost mutualinfo -searchcost mutualinfo
flirt -ref "${SCRIPTPATH}/../TractSeg/tractseg/resources/MNI_FA_template.nii.gz" -in $IN_FILE -out "${IN_FILE_PREFIX}_MNI.nii.gz" -applyxfm -init "${IN_FILE_PATH}/FA_MNI.mat" -dof 6
/bin/cp "${IN_FILE_PREFIX}.bval" "${IN_FILE_PREFIX}_MNI.bval"
rotate_bvecs -i "${IN_FILE_PREFIX}.bvec" -t "${IN_FILE_PATH}/FA_MNI.mat" -o "${IN_FILE_PREFIX}_MNI.bvec"

# Do tractography with TractSeg
TractSeg -i "${IN_FILE_PREFIX}_MNI.nii.gz" -o "$IN_FILE_PATH/tracts" --bvals "${IN_FILE_PREFIX}_MNI.bval" --bvecs "${IN_FILE_PREFIX}_MNI.bvec" --raw_diffusion_input
