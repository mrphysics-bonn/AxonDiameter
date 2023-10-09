
# input
if [ "$#" -lt 1 ]; then
    echo "Diffusion file missing"
    exit 1
elif [ "$#" -lt 2 ]; then
    echo "T1 file missing"
    exit 1
else
    IN_FILE="$1"
    T1_FILE="$2"
fi

IN_FILE_PREFIX=${IN_FILE%%.*}
IN_FILE_PATH=$(dirname $IN_FILE)
T1_FILE_PREFIX=${T1_FILE%%.*}
SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

# Create white matter mask from T1 data
if test -f "${T1_FILE_PREFIX}_bet.nii.gz" && test -f "${T1_FILE_PREFIX}_bet_seg.nii.gz"; then
    echo "White matter mask already exists. Skip calculation."
else
    # Brain extraction
    python "${SCRIPTPATH}/t1_processing.py" -m $T1_FILE "${T1_FILE_PREFIX}_bet"

    # Create white matter mask from T1 data
    fast -N "${T1_FILE_PREFIX}_bet.nii.gz"
    fslmaths "${T1_FILE_PREFIX}_bet_pve_2.nii.gz" -thr 0.85 -bin "${T1_FILE_PREFIX}_bet_seg.nii.gz"
fi

# Calculate relative SNR maps
${SCRIPTPATH}/relative_snr.sh "${IN_FILE_PREFIX}_moco_unwarped_bet.nii.gz" "${IN_FILE_PREFIX}_sh_b6000_powderavg.nii.gz" "${IN_FILE_PREFIX}_sh_b30000_powderavg.nii.gz" "${IN_FILE_PREFIX}_noise_map.nii.gz" "${T1_FILE_PREFIX}_bet.nii.gz"

# Quantify axon radius in white matter
${SCRIPTPATH}/wm_axons.sh "${IN_FILE_PREFIX}_moco_unwarped_bet.nii.gz" "${IN_FILE_PATH}/AxonRadiusMap.nii" ${T1_FILE}

# improved tractography with MrTrix & along-fibre quantification (only CST left atm)
${SCRIPTPATH}/along_tract/along_tract_CST.sh -f "${IN_FILE_PREFIX}_moco_unwarped_bet.nii.gz" "${T1_FILE%%.*}_bet.nii.gz" "${IN_FILE_PREFIX}_sh_b6000_powderavg.nii.gz" "${IN_FILE_PREFIX}_sh_b30000_powderavg.nii.gz" "${IN_FILE_PREFIX}_moco_unwarped_meanb0_bet_mask.nii.gz"

# recalculate axon diameters with along fibre quantified powder-averaged shells
PATH_TRACT=${IN_FILE_PATH}/along_tract/CST_L
matlab -nodisplay -r "addpath ${SCRIPTPATH}/../AxonRadiusMapping/;calcAxonAlongTracts('$PATH_TRACT/CST_stats_sh6000.txt', '$PATH_TRACT/CST_stats_sh30000.txt', '$PATH_TRACT/CST_grads_6000.txt', '$PATH_TRACT/CST_grads_30000.txt', 'CST_stats_Axon_afterAFQonShells.txt');exit"
matlab -nodisplay -r "addpath ${SCRIPTPATH}/../AxonRadiusMapping/;calcAxonAlongTracts('$PATH_TRACT/CST_stats_between_regions_sh6000.txt', '$PATH_TRACT/CST_stats_between_regions_sh30000.txt', '$PATH_TRACT/CST_grads_between_regions_6000.txt', '$PATH_TRACT/CST_grads_between_regions_30000.txt', 'CST_stats_between_regions_Axon_afterAFQonShells.txt');exit"
