#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 17 17:30:02 2022

@author: ferreiram
"""

#sudo pip3 install pyAFQ

import argparse
import os
import numpy as np
import nibabel as nib
from dipy.io.streamline import save_tractogram, load_tractogram
from dipy.io.stateful_tractogram import Space
from dipy.io.stateful_tractogram import StatefulTractogram

def main(args):

    dwi_path = args.dwi
    tck_path = args.tck
    tck_clipped_path = args.tck_clipped_path
    roi1_path = args.roi1
    roi2_path = args.roi2
    temporary_txt = args.temporary_txt

    dwi_nib = nib.load(dwi_path)
    dwi_data = dwi_nib.get_fdata()
    reference_anatomy = nib.load(dwi_path) 

    cc_sft = load_tractogram(tck_path, reference_anatomy, bbox_valid_check=False)

    if not len(cc_sft) == 0:
        # attention: if negative values once in voxel space, these files contain invalid streamlines. use arg bbox_valid_check=False
        cc_sft.to_vox()

        # remove invalid streamlines 
        cc_sft.remove_invalid_streamlines()
        os.system('maskdump ' + roi1_path + ' > ' + temporary_txt)
        txt_roi1 = np.loadtxt(temporary_txt) 
        txt_roi1[:,0] = dwi_data.shape[0] - txt_roi1[:,0]

        os.system('maskdump ' + roi2_path + ' > ' + temporary_txt)
        txt_roi2 = np.loadtxt(temporary_txt) 
        txt_roi2[:,0] = dwi_data.shape[0] - txt_roi2[:,0]
        all_new_streamlines = []
        
        
        if len(cc_sft.streamlines) > 1:
            for idx_streamlines in range(len(cc_sft.streamlines)): # enumerate streamlines
                select_sl = list(cc_sft.streamlines[idx_streamlines]) #one streamline selected
        
                dist_roi1 = np.zeros((len(select_sl),len(txt_roi1)))
                dist_roi2 = np.zeros((len(select_sl),len(txt_roi2)))                      
                
                for streamline_points in range(len(select_sl)):
                    a1 = select_sl[streamline_points]
                    matrix1 = np.ones([len(txt_roi1),3])
                    result1 =  a1 * matrix1
    
                    dist_roi1[streamline_points,:] = np.linalg.norm(result1 - txt_roi1, axis=1)
                    
                    a2 = select_sl[streamline_points]
                    matrix2 = np.ones([len(txt_roi2),3])
                    result2 =  a2 * matrix2
    
                    dist_roi2[streamline_points,:] = np.linalg.norm(result2 - txt_roi2, axis=1)

                minValue1 = np.min(dist_roi1)
                location1 = np.where(dist_roi1==minValue1)        
                minValue2 = np.min(dist_roi2)
                location2 = np.where(dist_roi2==minValue2)  
                
                
                if len(location1[0])>1:
                    location1 = location1[0]
                if len(location2[0])>1:
                    location2 = location2[0]
        
                locat1 = np.int(location1[0])# use int instead of np.int to silence teh warning
                locat2 = np.int(location2[0])
                
                if not locat1 == 0:
                    if locat1<locat2:
                        streamline = np.array(select_sl[locat1 :locat2])
                    elif locat1>locat2:
                        streamline = np.array(select_sl[locat2 :locat1])
                    else:
                        locat1 = locat1 - 1
                        locat2 = locat2 + 1
                        streamline = np.array(select_sl[locat1:locat2])
                
                    all_new_streamlines.append(np.array(streamline))    
                    
            
            new_tract = StatefulTractogram(all_new_streamlines,reference_anatomy,
                                                        Space.VOX)

            save_tractogram(new_tract,tck_clipped_path)
    else:
        raise ValueError("Empty tractogram. Exiting.")
    
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Clip end of tracts.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('dwi', type=str, help='Input DWI images.')
    parser.add_argument('tck', type=str, help='Input tractogram (.tck).')
    parser.add_argument('tck_clipped_path', type=str, help='Output clipped tck.')
    parser.add_argument('roi1', type=str, help='Input ROI 1.')
    parser.add_argument('roi2', type=str, help='Input ROI 2.')
    parser.add_argument('temporary_txt', type=str, help='Path of temporary txt file.')

    args = parser.parse_args()

    main(args)
