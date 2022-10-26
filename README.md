# AxonDiameter

Processing for spiral data can be started running scripts/process_spiral.sh #IN_FILE.
Processing for EPI data can be started running scripts/process_epi.sh #IN_FILE #IN_FILE_PA.

## Preprocessing & Axon diameter calculation

Steps 1 & 2 can be interchanged, if only magnitude data is available.

1. Denoising (MRtrix3: "dwidenoise") preferably on complex data [1]

2. Convert to magnitude data ("nifti2mag.py")
   
3. Check diffusion gradient table (MRtrix3: "dwigradcheck") [11]
   
4. Gibbs Ringing removal (MRtrix3: "mrdegibbs") [2]

5. EPI:
   - TOPUP (FSL) [3]
   - Eddy (FSL) [4]

6. Spiral
   - Motion correction (FSL Eddy)

7. Gradient nonlinearity correction (is essential due to high b-values/strong gradients) [5]
   - GradientDistortionUnwarp.sh (needs https://github.com/Washington-University/gradunwarp)
   - includes b-vector correction (gradient nonlinearity correction leads to different b-values/b-vectors in different voxels)

8. Brain masking (FSL BET) [6]

9.  Spherical harmonic decomposition to get spherical average per shell & per voxel (MRtrix3: amp2sh) [7]
	
10. Take the 0th order spherical harmonic and divide by $\sqrt{4\pi}$ to get the powder average [8]

11. Calculate axon radius maps [9,10]
    
12. Calculate relative SNR maps
   
13. Tractography using TractSeg [12]
    
14. White matter masking of axon radius maps (FSL FAST [13] & FSL FLIRT [14, 15]) 


## Requirements

- MRtrix3 (needs up to commit 3853c58 from https://github.com/lukeje/mrtrix3, that fixes a bug in the Rician bias correction)
- FSL
- gradunwarp (included submodule)
- AxonRadiusMapping (included submodule)
- TractSeg (included submodule, only for tractography)
- Python (incl. Numpy, Nibabel 3.2.2 (< version 4), Pytorch for TractSeg)
- Matlab

## References

1. Cordero-Grande, L. et. al. Complex diffusion-weighted image estimation via matrix recovery under general noise models. NeuroImage, 2019, 200, 391-404, doi: 10.1016/j.neuroimage.2019.06.039

2. Kellner, E; Dhital, B; Kiselev, V.G & Reisert, M. Gibbs-ringing artifact removal based on local subvoxel-shifts. Magnetic Resonance in Medicine, 2016, 76, 1574–1581

3. Jesper L. R. Andersson and Stamatios N. Sotiropoulos. An integrated approach to correction for off-resonance effects and subject movement in diffusion MR imaging. NeuroImage, 125:1063-1078, 2016. 

4. J.L.R. Andersson, S. Skare, J. Ashburner. How to correct susceptibility distortions in spin-echo echo-planar images: application to diffusion tensor imaging. NeuroImage, 20(2):870-888, 2003. 

5. Bammer, R., et al. Analysis and generalized correction of the effect of spatial gradient field distortions in diffusion‐weighted imaging. Magnetic Resonance in Medicine, 50(3):560-569, 2003

6. S.M. Smith. Fast robust automated brain extraction. Human Brain Mapping, 17(3):143-155, 2002.

7. Afzali, et al. Computing the orientational-average of diffusion-weighted MRI signals: a comparison of different techniques. Scientific Reports, 11:14345, 2021

8. Tournier, J.-D. et. al.. MRtrix3: A fast, flexible and open software framework for medical image processing and visualisation. NeuroImage, 2019, 202, 116137

9. Veraart, J. et. al. Noninvasive quantification of axon radii using diffusion MRI, eLife, 9:e49855, 2020

10. Veraart, J. et. al. The variability of MR axon radii estimates in the human white matter, Human Brain Mapping, 42:2201–2213, 2021
    
11. Jeurissen, B.; Leemans, A.; Sijbers, J. Automated correction of improperly rotated diffusion gradient orientations in diffusion weighted MRI. Medical Image Analysis, 2014, 18(7), 953-962
    
12. Wasserthal, J., Neher, P., Maier-Hein, K. H. TractSeg - Fast and accurate white matter tract segmentation, NeuroImage, 183:239-253, 2018
    
13. Zhang, Y. and Brady, M. and Smith, S. Segmentation of brain MR images through a hidden Markov random field model and the expectation-maximization algorithm. IEEE Trans Med Imag, 20(1):45-57, 2001.
    
14. M. Jenkinson and S.M. Smith. A global optimisation method for robust affine registration of brain images. Medical Image Analysis, 5(2):143-156, 2001. 

15. M. Jenkinson, P.R. Bannister, J.M. Brady, and S.M. Smith. Improved optimisation for the robust and accurate linear registration and motion correction of brain images. NeuroImage, 17(2):825-841, 2002. 
