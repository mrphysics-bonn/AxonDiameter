# AxonDiameter

Processing for spiral data can be started running scripts/process_spiral.sh #IN_FILE.

## Preprocessing

1. Denoising (MRtrix3: "dwidenoise") preferably on complex data

2. Convert to magnitude data ("nifti2mag.py")

3. EPI:
   - TOPUP (FSL)
   - Eddy (FSL)

4. Spiral
   - Motion correction (Eddy)

5. Gradient nonlinearity correction (is essential due to high b-values/strong gradients)
   - GradientDistortionUnwarp.sh (needs https://github.com/Washington-University/gradunwarp)
   - includes b-vector correction (gradient nonlinearity correction leads to different b-values/b-vectors in different voxels)
   - Bammer, R., et al., (2003), "Analysis and generalized correction of the effect of spatial gradient field distortions in diffusion‐weighted imaging"

6. Spherical harmonic decomposition to get spherical average per shell & per voxel (MRtrix3: amp2sh)
   - correct for Rician noise bias
	
7. Divide the 0th order spherical harmonic by $\sqrt{4\pi}$ to get the powder average 
   - Afzali, et al. "Computing the orientational-average of diffusion-weighted MRI signals: a comparison of different techniques"


## Requirements

- MRtrix3
- FSL
- gradunwarp (included submodule)
- AxonRadiusMapping (included submodule)
- Python (incl. Numpy, Nibabel)