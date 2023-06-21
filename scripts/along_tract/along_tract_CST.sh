#!/bin/bash

set -e

# input
if [ "$#" -lt 5 ]; then
    echo "Not enough input arguments. Expected 4 inputs: DWI_DATA T1_DATA POWDER_AVG_SHELL1 POWDER_AVG_SHELL2 BRAIN_MASK"
    exit 1
else
    IN_FILE="$1"
    T1_FILE="$2"
    PA_FILE1="$3"
    PA_FILE2="$4"
    MASK_FILE="$5"
fi

IN_FILE_PATH=$(dirname $IN_FILE)
IN_FILE_PREFIX=${IN_FILE%%.*}
SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
saving_dir="$IN_FILE_PATH/along_tract"
working_dir="$saving_dir/working_dir"
mkdir -p $saving_dir
mkdir -p $working_dir

# T1 MNI template
MNI_T1_path=$SCRIPTPATH/templates/mni_icbm152_t1_tal_nlin_asym_09c.nii

## ---------------------------------------registered data
dwi_mif="${IN_FILE_PREFIX}.mif"
dwi_nii="${IN_FILE_PREFIX}.nii.gz"
T1_in_dwi_space_nii=$working_dir/T1_in_dwi_space.nii.gz
T1_in_dwi_space_mif=$working_dir/T1_in_dwi_space.mif
dwi_mean_nii="${IN_FILE_PREFIX}_meanb0.nii.gz"

# creating the response function txt files 
WM_txt=$working_dir/wm_response.txt
GM_txt=$working_dir/gm_response.txt
CSF_txt=$working_dir/csf_response.txt
if [ ! -f $WM_txt ]
then
    echo "Starting with dwi2response."
    dwi2response dhollander -mask $MASK_FILE $dwi_mif $WM_txt $GM_txt $CSF_txt
fi

# creating the fiber orientation distributions for each tissue
WM_fod=$working_dir/wm_fod.mif
GM_fod=$working_dir/gm_fod.mif
CSF_fod=$working_dir/csf_fod.mif
if [ ! -f $WM_fod ]
then
    echo "Starting with dwi2fod."
    dwi2fod msmt_csd -mask $MASK_FILE $dwi_mif $WM_txt $WM_fod $GM_txt $GM_fod $CSF_txt $CSF_fod
fi

# ----- Get waypoints(=ROIs) -----
ROIs_folder=$working_dir/ROIs
mkdir -p $ROIs_folder
type_of_transform_matrix=SyN

# Register T1 to DWI
if [ ! -f $T1_in_dwi_space_mif ]
then
    python $SCRIPTPATH/coreg_T1_dwi.py $dwi_mean_nii $T1_FILE $T1_in_dwi_space_nii $type_of_transform_matrix
    mrconvert $T1_in_dwi_space_nii $T1_in_dwi_space_mif
fi

# Inclusion ROIs
CST_1L_orig_file=$SCRIPTPATH/templates/CST_roi1_L.nii.gz
CST_2L_orig_file=$SCRIPTPATH/templates/CST_roi2_L.nii.gz
CST_1R_orig_file=$SCRIPTPATH/templates/CST_roi1_R.nii.gz
CST_2R_orig_file=$SCRIPTPATH/templates/CST_roi2_R.nii.gz
# left 1 
CST_1L_roi=$ROIs_folder/CST_1_left.nii.gz
CST_1L_roi_mif=$ROIs_folder/CST_1_left.mif
if [ ! -f $CST_1L_roi_mif ]
then
    echo "Register first CST ROI to DWI"
    cp $CST_1L_orig_file $CST_1L_roi
    python $SCRIPTPATH/coreg_roi2dwi.py $T1_in_dwi_space_nii $MNI_T1_path $CST_1L_roi $type_of_transform_matrix
    mrconvert $CST_1L_roi $CST_1L_roi_mif
fi
# left 2
CST_2L_roi=$ROIs_folder/CST_2_left.nii.gz
CST_2L_roi_mif=$ROIs_folder/CST_2_left.mif
if [ ! -f $CST_2L_roi_mif ]
then
    echo "Register second CST ROI to DWI"
    cp $CST_2L_orig_file $CST_2L_roi
    python $SCRIPTPATH/coreg_roi2dwi.py $T1_in_dwi_space_nii $MNI_T1_path $CST_2L_roi $type_of_transform_matrix
    mrconvert $CST_2L_roi $CST_2L_roi_mif
fi

# ----- streamlines generation - CST ------
# left
CST_L_folder=$saving_dir/CST_L
mkdir -p $CST_L_folder
CST_L=$CST_L_folder/CST_L.tck
if [ ! -f $CST_L ]
then
    echo "Start generating streamlines."
    tckgen -seed_image $CST_1L_roi_mif -include $CST_1L_roi_mif -include $CST_2L_roi_mif $WM_fod $CST_L -minlength 50 -cutoff 0.1 -angle 30 -select 3000 -seeds 25M -force
fi
# ----- streamlines 1) cleaning; 2) cropping between 1st and last ROIs ------

cleaned_CST_L=$CST_L_folder/CST_L_cleaned.tck

min_sl=1000
clean_rounds=6
distance_threshold=3
length_threshold=3
n_points=100

if [ ! -f $cleaned_CST_L ]
then
    echo "Start cleaning streamlines."
    python3 $SCRIPTPATH/cleaning.py $dwi_nii $CST_L $cleaned_CST_L $min_sl $clean_rounds $distance_threshold $length_threshold $n_points
fi

temporary_txt=$CST_L_folder/clipped.txt
clipped_CST_L=$CST_L_folder/CST_L_clipped.tck
if [ ! -f $clipped_CST_L ]
then
    echo "Start clipping streamlines."
    python3 $SCRIPTPATH/clipping_between_regions.py $dwi_nii $cleaned_CST_L $clipped_CST_L $CST_1L_roi_mif $CST_2L_roi_mif $temporary_txt
fi

CST_L_mapped=$CST_L_folder/CST_L_map.nii
cleaned_CST_L_mapped=$CST_L_folder/CST_L_cleaned_map.nii
clipped_CST_L_mapped=$CST_L_folder/CST_L_clipped_map.nii
if [ ! -f $CST_L_mapped ]
then
    echo "Start tract-density mapping"
    tckmap -template $dwi_nii $CST_L $CST_L_mapped
    tckmap -template $dwi_nii $cleaned_CST_L $cleaned_CST_L_mapped
    tckmap -template $dwi_nii $clipped_CST_L $clipped_CST_L_mapped
fi
# # ----- profiling -----

echo "Start profiling"

# Axon Radius map
axon_path=$IN_FILE_PATH/AxonRadiusMap.nii
path_out0=$CST_L_folder/CST_stats_Axon.txt
if [ ! -f $path_out0 ]
then
    python3 $SCRIPTPATH/profile_fiber.py $axon_path $dwi_nii $cleaned_CST_L $path_out0 $n_points
fi

axon_path=$IN_FILE_PATH/AxonRadiusMap.nii
path_out_reg0=$CST_L_folder/CST_stats_between_regions_Axon.txt
if [ ! -f $path_out_reg0 ]
then
    python3 $SCRIPTPATH/profile_fiber.py $axon_path $dwi_nii $clipped_CST_L $path_out_reg0 $n_points
fi

# powder-averaged shell b=6000
grads6000_path=$IN_FILE_PATH/grads_6000.nii
path_out1=$CST_L_folder/CST_stats_sh6000.txt
path_out_grads1=$CST_L_folder/CST_grads_6000.txt
if [ ! -f $path_out1 ]
then
    python3 $SCRIPTPATH/profile_fiber.py $PA_FILE1 $dwi_nii $cleaned_CST_L $path_out1 $n_points
    python3 $SCRIPTPATH/profile_fiber.py $grads6000_path $dwi_nii $cleaned_CST_L $path_out_grads1 $n_points
fi

path_out_reg1=$CST_L_folder/CST_stats_between_regions_sh6000.txt
path_out_grads_reg1=$CST_L_folder/CST_grads_between_regions_6000.txt
if [ ! -f $path_out_reg1 ]
then
    python3 $SCRIPTPATH/profile_fiber.py $PA_FILE1 $dwi_nii $clipped_CST_L $path_out_reg1 $n_points
    python3 $SCRIPTPATH/profile_fiber.py $grads6000_path $dwi_nii $clipped_CST_L $path_out_grads_reg1 $n_points
fi

# powder-averaged shell b=30450
grads30000_path=$IN_FILE_PATH/grads_30000.nii
path_out2=$CST_L_folder/CST_stats_sh30000.txt
path_out_grads2=$CST_L_folder/CST_grads_30000.txt
if [ ! -f $path_out2 ]
then
    python3 $SCRIPTPATH/profile_fiber.py $PA_FILE2 $dwi_nii $cleaned_CST_L $path_out2 $n_points
    python3 $SCRIPTPATH/profile_fiber.py $grads30000_path $dwi_nii $cleaned_CST_L $path_out_grads2 $n_points
fi

path_out_reg2=$CST_L_folder/CST_stats_between_regions_sh30000.txt
path_out_grads_reg2=$CST_L_folder/CST_grads_between_regions_30000.txt
if [ ! -f $path_out_reg2 ]
then
    python3 $SCRIPTPATH/profile_fiber.py $PA_FILE2 $dwi_nii $clipped_CST_L $path_out_reg2 $n_points
    python3 $SCRIPTPATH/profile_fiber.py $grads30000_path $dwi_nii $clipped_CST_L $path_out_grads_reg2 $n_points
fi
