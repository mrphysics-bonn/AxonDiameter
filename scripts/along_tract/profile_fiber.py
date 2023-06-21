#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 17 17:30:02 2022

@author: ferreiram
"""

#sudo pip3 install pyAFQ
import argparse
import numpy as np
import nibabel as nib
from dipy.io.streamline import load_tractogram
from dipy.stats.analysis import afq_profile, gaussian_weights
from dipy.io.stateful_tractogram import Space
from dipy.io.utils import (create_nifti_header, get_reference_info)

def main(args):
    map_path = args.map_path
    dwi_path = args.dwi
    tck_path = args.tck
    path_out = args.path_out
    n_points = args.n_points

    map_nib = nib.load(map_path)
    map_data = map_nib.get_fdata()

    reference_anatomy = nib.load(dwi_path) 

    fileobj_clipped = load_tractogram(tck_path,reference_anatomy, to_space=Space.VOX)
    # dwi_data = reference_anatomy.get_fdata()

    affine, dimensions, voxel_sizes, voxel_order = get_reference_info(reference_anatomy)
    # nifti_header = create_nifti_header(affine, dimensions, voxel_sizes)
                    
    weights = gaussian_weights(fileobj_clipped.streamlines, n_points=n_points)
    profile = afq_profile(map_data, fileobj_clipped.streamlines, np.eye(4), weights = weights, n_points = n_points)

    np.savetxt(path_out, profile)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Do profiling along fibers with specified quantitative map.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('map_path', type=str, help='Input map, which will be used for profiling.')
    parser.add_argument('dwi', type=str, help='Input DWI.')
    parser.add_argument('tck', type=str, help='Input tractogram (.tck).')
    parser.add_argument('path_out', type=str, help='Output path.')
    parser.add_argument('n_points', type=int, help='Number of points for profiling.')

    args = parser.parse_args()

    main(args)