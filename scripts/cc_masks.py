
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Calculate masks for Corpus Callosum in T1 space
Uses "Hammers_mith-n30r95-MaxProbMap-full-MNI152-SPM12.nii.gz" (CC is label 44)
Hammersmith atlas http://brain-development.org/brain-atlases/adult-brain-atlases/
with MNI T1 template
"""

import os
import argparse
import ants
from scipy.ndimage import binary_erosion
import nibabel as nib
from erode_mask import kernel2d, kernel3d

def main(args):

    scriptpath = os.path.dirname(os.path.abspath(__file__))

    path_to_T1 = args.t1_file
    path_to_mni = os.path.join(scriptpath, "cc_template/MNI152_T1_1mm.nii.gz")
    path_to_mni_cc = os.path.join(scriptpath, "cc_template/cc_mask_mni.nii.gz")

    # Read images
    mni = ants.image_read(path_to_mni)
    cc_mni = ants.image_read(path_to_mni_cc)
    t1 = ants.image_read(path_to_T1)
    t1_mask = ants.image_read(path_to_T1)
    t1_mask[t1_mask>0] = 1

    # Register MNI to T1 and apply transform to CC mask
    mni_to_t1 = ants.registration(fixed=t1, moving=mni, mask=t1_mask, type_of_transform = 'SyNRA')
    cc = ants.apply_transforms(fixed=t1, moving=cc_mni, transformlist=mni_to_t1['fwdtransforms'])
    
    # Threshold and erode
    thresh = 0.8
    cc_numpy = cc.numpy()
    cc_numpy[cc_numpy<thresh] = 0
    cc_numpy[cc_numpy!=0] = 1
    # cc_numpy = binary_erosion(cc_numpy, kernel3d(), iterations=1)
    for j in range(cc_numpy.shape[-1]):
        cc_numpy[...,j] = binary_erosion(cc_numpy[...,j], kernel2d(), iterations=1)

    out_file = os.path.join(os.path.dirname(path_to_T1), "cc_mask_bin.nii.gz")
    t1_nifti = nib.load(path_to_T1)
    nifti_img = nib.Nifti1Image(cc_numpy, affine=t1_nifti.affine, header=t1_nifti.header)
    nib.save(nifti_img, out_file)

    out_file = os.path.join(os.path.dirname(path_to_T1), "cc_mask.nii.gz")
    nifti_img = nib.Nifti1Image(cc.numpy(), affine=t1_nifti.affine, header=t1_nifti.header)
    nib.save(nifti_img, out_file)

    out_file = os.path.join(os.path.dirname(path_to_T1), "mni_in_t1.nii.gz")
    nifti_img = nib.Nifti1Image(mni_to_t1['warpedmovout'].numpy(), affine=t1_nifti.affine, header=t1_nifti.header)
    nib.save(nifti_img, out_file)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Calculate mask for T1.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('t1_file', type=str, help='Input T1 file.')

    args = parser.parse_args()

    main(args)

