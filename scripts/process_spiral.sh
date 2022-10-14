#!/bin/bash

# input
if [ "$#" -lt 1 ]; then
    echo "Input file missing"
    exit 1
else
    IN_FILE="$1"
fi

mag_flag=false
while getopts 'm' mag; do
mag_flag=true
done

IN_FILE_PREFIX=${IN_FILE%%.*}
IN_FILE_PATH=$(dirname $IN_FILE)
SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

if $mag_flag; then
    # Convert to magnitude
    python "${SCRIPTPATH}/nifti2mag.py" ${IN_FILE} "${IN_FILE_PREFIX}_mag.nii.gz"
    # Denoising
    dwidenoise -noise "${IN_FILE_PREFIX}_noise_map.nii.gz" "${IN_FILE_PREFIX}_mag.nii.gz" "${IN_FILE_PREFIX}_denoised_mag.nii.gz"
    # create input for rician bias correction
    mrconvert "${IN_FILE_PREFIX}_mag.nii.gz" -fslgrad "${IN_FILE_PREFIX}.bvec" "${IN_FILE_PREFIX}.bval" "${IN_FILE_PREFIX}_ricinput.mif"
else
    # Denoising
    dwidenoise -noise "${IN_FILE_PREFIX}_noise_map.nii.gz" ${IN_FILE} "${IN_FILE_PREFIX}_denoised.nii.gz"
    # Convert to magnitude
    python "${SCRIPTPATH}/nifti2mag.py" "${IN_FILE_PREFIX}_denoised.nii.gz" "${IN_FILE_PREFIX}_denoised_mag.nii.gz"
    # create input for rician bias correction
    mrconvert ${IN_FILE} -fslgrad "${IN_FILE_PREFIX}.bvec" "${IN_FILE_PREFIX}.bval" "${IN_FILE_PREFIX}_ricinput.mif"
fi

# convert nii to mif
mrconvert "${IN_FILE_PREFIX}_denoised_mag.nii.gz" -fslgrad "${IN_FILE_PREFIX}.bvec" "${IN_FILE_PREFIX}.bval" "${IN_FILE_PREFIX}_denoised_mag.mif"

# Low-SNR Rician bias correction (from DESIGNER pipeline)
dwiextract -shell 30450 "${IN_FILE_PREFIX}_ricinput.mif" dwilowb.mif
dwidenoise -noise lowbnoisemap.mif dwilowb.mif dwitmp.mif
mrcalc "${IN_FILE_PREFIX}_denoised_mag.mif" 2 -pow lowbnoisemap.mif 2 -pow -sub -abs -sqrt - | mrcalc - -finite - 0 -if "${IN_FILE_PREFIX}_denoised_mag_rician.mif"
/bin/rm dwilowb.mif lowbnoisemap.mif dwitmp.mif dwirc.mif "${IN_FILE_PREFIX}_ricinput.mif"

# Gibbs-Ringing removal
mrdegibbs "${IN_FILE_PREFIX}_denoised_mag_rician.mif" "${IN_FILE_PREFIX}_denoised_mag_gr.mif"

# this file contains a list of the simultaneously acquired slices in acquisition order
slspec="$SCRIPTPATH/example_slspec.txt"

# Motion correction with eddy
# let mrtrix take care of providing eddy input and output
# readout_time may need to be longer to satisfy eddy's sanity checking (but not so long that it starts expecting large shifts)
# first level model (flm) is set to linear to avoid overfitting as spiral data should already have been corrected for eddy currents
# mporder is recommended to be somewhere between N/4 and N/2, where N is the number of excitations
dwifslpreproc "${IN_FILE_PREFIX}_denoised_mag_gr.mif" "${IN_FILE_PREFIX}_moco.mif" -rpe_none -pe_dir ap -readout_time 0.01 -eddy_slspec $slspec -eddyqc_all $IN_FILE_PATH -eddy_options " --flm=linear --repol --data_is_shelled --mporder=13 --ol_type=both "

# Convert mrtrix output to nii and bvec/bval
mrconvert "${IN_FILE_PREFIX}_moco.mif" -export_grad_fsl "${IN_FILE_PREFIX}_moco.bvec" "${IN_FILE_PREFIX}_moco.bval" "${IN_FILE_PREFIX}_moco.nii.gz"

# Gradient nonlinearity correction
${SCRIPTPATH}/GradientDistortionUnwarp.sh --workingdir="${SCRIPTPATH}/../data/unwarp_wd" --in="${IN_FILE_PREFIX}_moco" --out="${IN_FILE_PREFIX}_moco_unwarped" --coeffs="${SCRIPTPATH}/../connectom_coeff.grad" --owarp="${IN_FILE_PREFIX}_owarp"

# Brain masking
fslsplit "${IN_FILE_PREFIX}_moco_unwarped.nii.gz" "${IN_FILE_PREFIX}_splitted_vol"
bet "${IN_FILE_PREFIX}_splitted_vol0000.nii.gz" "${IN_FILE_PREFIX}_splitted_vol0000.nii.gz" -f 0.3 -m
fslmaths "${IN_FILE_PREFIX}_moco_unwarped.nii.gz" -mul "${IN_FILE_PREFIX}_splitted_vol0000_mask.nii.gz" "${IN_FILE_PREFIX}_moco_unwarped.nii.gz"
/bin/rm "${IN_FILE_PREFIX}_splitted_vol"*

# Spherical harmonic decomposition - disable Rician noise correction as its currently not correct
amp2sh -lmax 6 -shells 0,6000 -normalise -fslgrad "${IN_FILE_PREFIX}_moco.bvec" "${IN_FILE_PREFIX}_moco.bval" "${IN_FILE_PREFIX}_moco_unwarped.nii.gz" "${IN_FILE_PREFIX}_sh_b6000.nii.gz"
amp2sh -lmax 6 -shells 0,30450 -normalise -fslgrad "${IN_FILE_PREFIX}_moco.bvec" "${IN_FILE_PREFIX}_moco.bval" "${IN_FILE_PREFIX}_moco_unwarped.nii.gz" "${IN_FILE_PREFIX}_sh_b30000.nii.gz"

# Divide by sqrt(4pi) to get powder average
fslmaths "${IN_FILE_PREFIX}_sh_b6000.nii.gz" -div 3.5449077018110318 "${IN_FILE_PREFIX}_sh_b6000_powderavg.nii.gz"
fslmaths "${IN_FILE_PREFIX}_sh_b30000.nii.gz" -div 3.5449077018110318 "${IN_FILE_PREFIX}_sh_b30000_powderavg.nii.gz"

# Calculate axon diameters
matlab -nodisplay -r "addpath ${SCRIPTPATH}/../AxonRadiusMapping/;calcAxonMaps('${IN_FILE_PREFIX}_sh_b6000_powderavg.nii.gz', '${IN_FILE_PREFIX}_sh_b30000_powderavg.nii.gz', '${IN_FILE_PREFIX}.bval', '${IN_FILE_PREFIX}.bvec', '${IN_FILE_PATH}/grad_dev.nii.gz');exit"
