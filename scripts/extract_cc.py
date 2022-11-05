#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Extract the Corpus Callosum from tract masks from TractSeg
"""

import argparse
import os
import nibabel as nib
import numpy as np
from skimage.measure import label

def getLargestCC(segmentation):
    """
    gets the largest connected component
    """
    
    labels = label(segmentation)
    assert( labels.max() != 0 ) # assume at least 1 CC
    largestCC = labels == np.argmax(np.bincount(labels.flat)[1:])+1
    return largestCC

def main(args):

    if not os.path.exists(args.out_path):
        os.makedirs(args.out_path)
        print(f"Created folder {args.out_path} for output files")

    cc_segments = 7

    # get header
    file_cc1 = nib.load(args.in_path1+"/CC.nii.gz")
    hdr = file_cc1.header
    affine = file_cc1.affine

    # read data
    data1 = [nib.load(args.in_path1+f"/CC_{k+1}.nii.gz").get_fdata() for k in range(cc_segments)]
    data1.append(file_cc1.get_fdata())
    data_cc_b = nib.load(args.in_path2+"/CC_b.nii.gz").get_fdata()
    data_cc_e = nib.load(args.in_path2+"/CC_e.nii.gz").get_fdata()

    data_cc = [(data1[k] - data_cc_b - data_cc_e) for k in range(len(data1))]
    for k,data in enumerate(data_cc):
        data[data<0] = 0 # endings mask seems to contain gray matter, remove that from the mask
        data = getLargestCC(data).astype(data.dtype) # remove non-connected voxels
        hdr.set_data_dtype(data.dtype)
        img_out = nib.Nifti1Image(data, affine=affine, header=hdr)
        if k == 7:
            nib.save(img_out, args.out_path+"/CC_noEndings.nii.gz")
        else:
            nib.save(img_out, args.out_path+f"/CC_{k+1}_noEndings.nii.gz")

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Convert complex Nifti to magnitude images.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('in_path1', type=str, help='Input path of complete tract segmentation.')
    parser.add_argument('in_path2', type=str, help='Input path of tract endings segmentation.')
    parser.add_argument('out_path', type=str, help='Output path.')

    args = parser.parse_args()

    main(args)
