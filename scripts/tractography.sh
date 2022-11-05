#!/bin/bash

# option
preproc_flag=false
endings_flag=false
while getopts "pe" opt; do
    case $opt in
        p) preproc_flag=true ;;
        e) endings_flag=true ;;
    esac
done
shift "$((OPTIND-1))"

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

if $preproc_flag; then
    echo "Images are registered to MNI space before segmentation."

    # Register diffusion images to MNI space with FA map - not recommended, does not seem to work well
    fslsplit $IN_FILE "${IN_FILE_PREFIX}_splitted_vol"
    bet "${IN_FILE_PREFIX}_splitted_vol0000.nii.gz" "${IN_FILE_PREFIX}_splitted_vol0000.nii.gz" -f 0.3 -m
    /bin/cp "${IN_FILE_PREFIX}_splitted_vol0000_mask.nii.gz" "${IN_FILE_PREFIX}_mask.nii.gz"
    /bin/rm "${IN_FILE_PREFIX}_splitted_vol"*
    calc_FA -i $IN_FILE -o "${IN_FILE_PATH}/FA.nii.gz" --bvals "${IN_FILE_PREFIX}.bval" --bvecs "${IN_FILE_PREFIX}.bvec"  --brain_mask "${IN_FILE_PREFIX}_mask.nii.gz"
    flirt -ref "${SCRIPTPATH}/../TractSeg/tractseg/resources/MNI_FA_template.nii.gz" -in "${IN_FILE_PATH}/FA.nii.gz" -out "${IN_FILE_PATH}/FA_MNI.nii.gz" -omat "${IN_FILE_PATH}/FA_MNI.mat" -dof 6 -cost mutualinfo -searchcost mutualinfo
    flirt -ref "${SCRIPTPATH}/../TractSeg/tractseg/resources/MNI_FA_template.nii.gz" -in $IN_FILE -out "${IN_FILE_PREFIX}_MNI.nii.gz" -applyxfm -init "${IN_FILE_PATH}/FA_MNI.mat" -dof 6
    /bin/cp "${IN_FILE_PREFIX}.bval" "${IN_FILE_PREFIX}_MNI.bval"
    rotate_bvecs -i "${IN_FILE_PREFIX}.bvec" -t "${IN_FILE_PATH}/FA_MNI.mat" -o "${IN_FILE_PREFIX}_MNI.bvec"

    # Extract shells 0,30450 for tractography
    mrconvert -force "${IN_FILE_PREFIX}_MNI.nii.gz" -fslgrad "${IN_FILE_PREFIX}_MNI.bvec" "${IN_FILE_PREFIX}_MNI.bval" "${IN_FILE_PREFIX}_MNI.mif"
    dwiextract -shells 0,30450 "${IN_FILE_PREFIX}_MNI.mif" "${IN_FILE_PREFIX}_b30000_MNI.mif"
    mrconvert -force "${IN_FILE_PREFIX}_b30000_MNI.mif" -export_grad_fsl "${IN_FILE_PREFIX}_b30000_MNI.bvec" "${IN_FILE_PREFIX}_b30000_MNI.bval" "${IN_FILE_PREFIX}_b30000_MNI.nii.gz"

    # Do tractography with TractSeg
    if $endings_flag; then
        echo "Calculate Tract endings"
        /bin/rm -rf "$IN_FILE_PATH/tracts_endings_MNI"
        TractSeg -i "${IN_FILE_PREFIX}_b30000_MNI.nii.gz" -o "$IN_FILE_PATH/tracts_endings_MNI" --bvals "${IN_FILE_PREFIX}_b30000_MNI.bval" --bvecs "${IN_FILE_PREFIX}_b30000_MNI.bvec" --raw_diffusion_input --output_type=endings_segmentation
    else
        /bin/rm -rf "$IN_FILE_PATH/tracts_MNI"
        TractSeg -i "${IN_FILE_PREFIX}_b30000_MNI.nii.gz" -o "$IN_FILE_PATH/tracts_MNI" --bvals "${IN_FILE_PREFIX}_b30000_MNI.bval" --bvecs "${IN_FILE_PREFIX}_b30000_MNI.bvec" --raw_diffusion_input
    fi
else
    # Extract shells 0,30450 for tractography
    mrconvert -force $IN_FILE -fslgrad "${IN_FILE_PREFIX}.bvec" "${IN_FILE_PREFIX}.bval" "${IN_FILE_PREFIX}.mif"
    dwiextract -force -shells 0,30450 "${IN_FILE_PREFIX}.mif" "${IN_FILE_PREFIX}_b30000.mif"
    mrconvert -force "${IN_FILE_PREFIX}_b30000.mif" -export_grad_fsl "${IN_FILE_PREFIX}_b30000.bvec" "${IN_FILE_PREFIX}_b30000.bval" "${IN_FILE_PREFIX}_b30000.nii.gz"

    # Do tractography with TractSeg
    if $endings_flag; then
        echo "Calculate Tract endings"
        /bin/rm -rf "$IN_FILE_PATH/tracts_endings"
        TractSeg -i "${IN_FILE_PREFIX}_b30000.nii.gz" -o "$IN_FILE_PATH/tracts_endings" --bvals "${IN_FILE_PREFIX}_b30000.bval" --bvecs "${IN_FILE_PREFIX}_b30000.bvec" --raw_diffusion_input --output_type=endings_segmentation
    else
        /bin/rm -rf "$IN_FILE_PATH/tracts"
        TractSeg -i "${IN_FILE_PREFIX}_b30000.nii.gz" -o "$IN_FILE_PATH/tracts" --bvals "${IN_FILE_PREFIX}_b30000.bval" --bvecs "${IN_FILE_PREFIX}_b30000.bvec" --raw_diffusion_input
    fi
fi
