# run with MRtrix3 

# input files (please fill in!)
DWI_in=
bvec_in=
bval_in=

# output directory (please fill in!)
outdir=

# path to slspec.txt (please fill in!)
# this file contains a list of the simultaneously acquired slices in acquisition order
slspec= 

# convert nii to mif; note must also provide bvecs and bvals
mrconvert $DWI_in -fslgrad $bvec_in $bval_in $outdir/DWI_in.mif

# let mrtrix take care of providing eddy input and output
# readout_time may need to be longer to satisfy eddy's sanity checking (but not so long that it starts expecting large shifts)
# first level model (flm) is set to linear to avoid overfitting as spiral data should already have been corrected for eddy currents
dwifslpreproc $outdir/DWI_in.mif $outdir/DWI_out.mif -rpe_none -pe_dir ap -readout_time 0.001 -eddy_slspec $slspec -eddyqc_all $outdir -eddy_options " --flm=linear --repol --data_is_shelled "

# convert mrtrix output to nii and bvec/bval
mrconvert $outdir/DWI_out.mif -export_grad_fsl $outdir/DWI_out.bvec $outdir/DWI_out.bval $outdir/DWI_out.nii.gz
