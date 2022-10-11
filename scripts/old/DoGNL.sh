export PATH=$PATH:$scriptpath/gradunwarp/bin/
export PYTHONPATH=$PYTHONPATH:$scriptpath/gradunwarp/lib/python3.6/site-packages/

cwd=$(pwd)

for infile in $rootdir/rawdata/sub*/dti/*_nonDWI-PA.nii.gz; do

	outfolder=$(dirname $infile)
	outwarp=$outfolder/warp
	outref=$outfolder/grad_ref.nii.gz

	WD=$outfolder/WD
	mkdir $WD

	# move (temporarily) into the working directory as gradient_unwarp.py outputs some files directly into pwd
	cd $WD

	# Extract first volume and run gradient distortion correction on this (all others follow suit as scanner coordinate system is unchanged, even with subject motion)
	fslroi $infile $outref 0 1

	# NB: gradient_unwarp.py *must* have the filename extensions written out explicitly or it will crash
	gradient_unwarp.py $outref trilinear.nii.gz siemens -g $(realpath $scriptpath/connectom_coeff.grad) -n --numpoints 128 --interp_order 2 #--nojacobian # default 60, 1

	# Now create an appropriate warpfield output (relative convention) and apply it to all timepoints
	#convertwarp's jacobian output has 8 frames, each combination of one-sided differences, so average them
	convertwarp --abs --ref=$WD/trilinear.nii.gz --warp1=$WD/fullWarp_abs.nii.gz --relout --out=$outwarp --jacobian=${outwarp}_jacobian
	fslmaths ${outwarp}_jacobian -Tmean ${outwarp}_jacobian

	#applywarp --rel --interp=spline -i $InputFile -r $outref -w $outwarp -o $OutputFile

	# Calculate gradent deviations
	GRADDEV=$outfolder/grad_dev.nii 
	GRADDEV_X=$outfolder/grad_dev_x.nii
	GRADDEV_Y=$outfolder/grad_dev_y.nii
	GRADDEV_Z=$outfolder/grad_dev_z.nii
	IMRM=$outfolder/grad_dev_?	

	# Based on:
	# https://github.com/Washington-University/HCPpipelines/blob/master/DiffusionPreprocessing/scripts/eddy_postproc.sh
	calc_grad_perc_dev --fullwarp=$outwarp --out=$GRADDEV
	fslmerge -t $GRADDEV $GRADDEV_X $GRADDEV_Y $GRADDEV_Z
	fslmaths $GRADDEV -div 100 $GRADDEV # Convert from % deviation to absolute
	imrm $IMRM

	cd $cwd
done