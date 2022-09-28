#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Convert complex Nifti to magnitude images
"""

import nibabel as nib
import numpy as np
import argparse

def main(args):

    file = nib.load(args.in_file)
    data = file.get_fdata(dtype=np.complex64)
    data = abs(data)
    hdr = file.header
    hdr.set_data_dtype(data.dtype)
    img_new = nib.Nifti1Image(data, affine=file.affine, header=hdr)
    nib.save(img_new, args.out_file)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Convert complex Nifti to magnitude images.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('in_file', type=str, help='Input Nifti file.')
    parser.add_argument('out_file', type=str, help='Output Nifti file.')

    args = parser.parse_args()

    main(args)
