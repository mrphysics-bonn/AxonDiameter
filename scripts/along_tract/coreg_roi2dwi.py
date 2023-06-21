#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 11:31:31 2022

@author: ferreiram
"""

import argparse
import ants
import nibabel as nib
from nipype.interfaces.fsl import UnaryMaths

def main(args):

    path_to_T1_in_dwi_space = args.t1_dwi_space
    T1_used_for_segmentation = args.t1_seg
    waypoint_path = args.waypoint_path
    type_of_transform_matrix = args.type_of_transform_matrix

    print("Starting converting the ROI nifti to the dwi space")

    reference = ants.image_read(path_to_T1_in_dwi_space)
    subject = ants.image_read(T1_used_for_segmentation)
    ROI = ants.image_read(waypoint_path)

    ## NORMALISATION OF REF/TEMPLATE TO SUBJECT SPACE - TO GET TRANSFORMATION MATRIX
    mytx = ants.registration(fixed = reference , moving = subject, type_of_transform = type_of_transform_matrix)
    print(mytx)
    # warped_moving = mytx['warpedmovout']

    ## APLICATION OF TRANSF MATRIX ON THE ROI TO BECOME IN THE SUBJECT SPACE
    mywarpedimage = ants.apply_transforms(fixed = reference, moving = ROI, transformlist = mytx['fwdtransforms'])

    nib.save(mywarpedimage,waypoint_path)

    maths_weight = UnaryMaths()
    maths_weight.inputs.in_file = waypoint_path
    maths_weight.inputs.operation = 'bin'
    maths_weight.inputs.out_file = waypoint_path
    print (maths_weight.cmdline)
    maths_weight.run()

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Coregister ROIs to DWI.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('t1_dwi_space', type=str, help='Input T1 in DWI space.')
    parser.add_argument('t1_seg', type=str, help='Input T1 for segmentation.')
    parser.add_argument('waypoint_path', type=str, help='Input waypoint path.')
    parser.add_argument('type_of_transform_matrix', type=str, help='Input type_of_transform_matrix.')

    args = parser.parse_args()

    main(args)