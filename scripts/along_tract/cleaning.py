#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 17 17:30:02 2022

@author: ferreiram
"""

#sudo pip3 install pyAFQ

import argparse
import os.path as op
import numpy as np
import nibabel as nib
from dipy.io.streamline import load_tractogram
                           
import AFQ.segmentation as seg
from dipy.io.utils import create_tractogram_header

from nibabel.streamlines import detect_format
from nibabel.streamlines.tractogram import Tractogram 

def main(args):

    dwi_path = args.dwi
    tck_path = args.tck
    tck_cleaned_path = args.tck_cleaned_path
    min_sl = args.min_sl
    clean_rounds = args.clean_rounds
    distance_threshold = args.distance_threshold
    length_threshold = args.length_threshold
    n_points = args.n_points

    # input example:
    # dwi_path = "/Volumes/cerebnet/ESMI_dwi_tracking/data/ataxSCA3_ESMI_AAC_418-918-735_20170425_BASL/processing/topup_eddy/dwi_eddy_corrected.nii.gz"
    # tck_path = "/Volumes/cerebnet/ESMI_dwi_tracking/2022/Tracts/DRTCT_L/ataxSCA3_ESMI_AAC_418-918-735_20170425_BASL.tck"
    # tck_cleaned_path = "/Volumes/cerebnet/ESMI_dwi_tracking/2022/Tracts/cleaned/DRTCT_L/ataxSCA3_ESMI_AAC_418-918-735_20170425_BASL.tck"

    # cleaning parameters (default)
    # -> n_points=100
    # -> clean_rounds=5           Number of rounds of cleaning based on the Mahalanobis 
    #                             distance from the mean of extracted bundles
    # -> distance_threshold=5     Threshold of cleaning based on the Mahalanobis distance 
    #                             (the units are standard deviations).
    # -> length_threshold=4       Threshold for cleaning based on length (in standard 
    #                             deviations). Lengthof any streamline should not be *more* 
    #                             than this number of stdevs from the mean length.
    # -> min_sl=20                Number of streamlines in a bundle under which we will
    #                             not bother with cleaning outliers. 

    dwi_nib = nib.load(dwi_path)

    tractogram = load_tractogram(tck_path, dwi_nib, trk_header_check = True)
    if not len(tractogram) == 0:
        print(f"Before cleaning: {len(tractogram)} streamlines")
        
        tractogram_cleaned = seg.clean_bundle(tractogram, n_points = n_points, clean_rounds = clean_rounds, 
            distance_threshold = distance_threshold,length_threshold = length_threshold, min_sl = min_sl, stat = 'mean',return_idx = False)
        
        print(f"After cleaning: {len(tractogram_cleaned)} streamlines")
    
        # creating header so that we can save the cleaned fibers properly
        tractogram_type = detect_format(tck_path)
        header = create_tractogram_header(tractogram_type, *tractogram_cleaned.space_attributes)
        new_tractogram = Tractogram(tractogram_cleaned.streamlines, affine_to_rasmm = np.eye(4))
    
    
        new_tractogram.data_per_point = tractogram_cleaned.data_per_point
        new_tractogram.data_per_streamline = tractogram_cleaned.data_per_streamline
    
        fileobj = tractogram_type(new_tractogram, header = header)
        
        nib.streamlines.save(fileobj, tck_cleaned_path)
    else:
        raise ValueError("Empty tractogram. Exiting.")

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Clean tracts.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('dwi', type=str, help='Input DWI images.')
    parser.add_argument('tck', type=str, help='Input tractogram (.tck).')
    parser.add_argument('tck_cleaned_path', type=str, help='Output cleaned tck.')
    parser.add_argument('min_sl', type=int, help='Minimum number of streamlines for cleaning.')
    parser.add_argument('clean_rounds', type=int, help='Number of cleaning rounds.')
    parser.add_argument('distance_threshold', type=int, help='Distance threshold for cleaning.')
    parser.add_argument('length_threshold', type=int, help='Length threshold for cleaning.')
    parser.add_argument('n_points', type=int, help='Number of points for cleaning.')

    args = parser.parse_args()

    main(args)
