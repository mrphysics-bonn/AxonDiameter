#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 17 17:30:02 2022

@author: ferreiram
"""

#sudo pip3 install pyAFQ
import argparse
import numpy as np

import dipy.stats.analysis as dsa
import nibabel as nib
import dipy.tracking.streamline as dts
from dipy.segment.clustering import QuickBundles
from dipy.segment.metric import AveragePointwiseEuclideanMetric
from dipy.segment.featurespeed import ResampleFeature

from dipy.io.streamline import load_tractogram
from dipy.io.image import load_nifti
from dipy.io.streamline import load_trk
from dipy.data.fetcher import get_two_hcp842_bundles
from dipy.data import fetch_bundle_atlas_hcp842

def main(args):
    map_path = args.map_path
    tck_path = args.tck
    path_out = args.path_out
    n_points = args.n_points

    map_nib = nib.load(map_path)
    map_data = map_nib.get_fdata()

    # load tractogram
    tractogram = load_tractogram(tck_path, reference=map_path, bbox_valid_check=False)
    streamlines = tractogram.streamlines

    # calculate centroid line from model bundle
    fetch_bundle_atlas_hcp842()
    _, model_cst_l_file = get_two_hcp842_bundles()
    model_cst_l = load_trk(model_cst_l_file, "same", bbox_valid_check=False).streamlines
    feature = ResampleFeature(nb_points=n_points)
    metric = AveragePointwiseEuclideanMetric(feature)
    qb = QuickBundles(np.inf, metric=metric)
    cluster_tract = qb.cluster(model_cst_l)
    standard_tract = cluster_tract.centroids[0]

    # orient tracts and calculate profile
    oriented_tract = dts.orient_by_streamline(streamlines, standard_tract)
    map_data, affine = load_nifti(map_path)
    weights = dsa.gaussian_weights(oriented_tract, n_points=n_points)
    profile = dsa.afq_profile(map_data, oriented_tract, affine, weights=weights, n_points = n_points)

    np.savetxt(path_out, profile)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Do profiling along fibers with specified quantitative map.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('map_path', type=str, help='Input map, which will be used for profiling.')
    parser.add_argument('tck', type=str, help='Input tractogram (.tck).')
    parser.add_argument('path_out', type=str, help='Output path.')
    parser.add_argument('n_points', type=int, help='Number of points for profiling.')

    args = parser.parse_args()

    main(args)