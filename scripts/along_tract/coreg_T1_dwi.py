#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 11:31:31 2022

@author: ferreiram
"""

import argparse
import ants
import nibabel as nib

def main(args):
    path_dwi_image = args.dwi
    path_t1_image = args.t1
    saving_path = args.saving_path
    type_of_transform_matrix = args.type_of_transform_matrix

    print("Start converting the T1 nifti to the dwi space")

    reference = ants.image_read(path_dwi_image)
    subject = ants.image_read(path_t1_image)

    ## NORMALISATION OF REF/TEMPLATE TO SUBJECT SPACE - TO GET TRANSFORMATION MATRIX
    mytx = ants.registration(fixed = reference , moving = subject, type_of_transform = type_of_transform_matrix)
    print(mytx)
    # warped_moving = mytx['warpedmovout']

    ## APLICATION OF TRANSF MATRIX ON THE ROI TO BECOME IN THE SUBJECT SPACE
    mywarpedimage = ants.apply_transforms(fixed = reference, moving = subject, transformlist = mytx['fwdtransforms'])
        
    nib.save(mywarpedimage,saving_path)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Coregister T1 to DWI.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('dwi', type=str, help='Input DWI images.')
    parser.add_argument('t1', type=str, help='Input T1.')
    parser.add_argument('saving_path', type=str, help='Saving path.')
    parser.add_argument('type_of_transform_matrix', type=str, help='Input type_of_transform_matrix.')

    args = parser.parse_args()

    main(args)