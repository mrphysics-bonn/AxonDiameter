#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Rotate b-vectors with rotation matrices from McFlirts motion correction
"""

import numpy as np
import argparse
import os

def rot_mat(theta):

    R_x = np.matrix([[ 1,              0,              0],
                     [ 0, np.cos(theta[0]), np.sin(theta[0])],
                     [ 0, -np.sin(theta[0]), np.cos(theta[0])],
                     ])

    R_y = np.matrix([[ np.cos(theta[1]), 0, -np.sin(theta[1])],
                     [ 0,              1,              0],
                     [ np.sin(theta[1]),  0, np.cos(theta[1])],
                     ])

    R_z = np.matrix([[ np.cos(theta[2]),  np.sin(theta[2]),  0],
                     [-np.sin(theta[2]),  np.cos(theta[2]),  0],
                     [ 0,              0,              1],
                     ])

    return R_x*R_y*R_z

def main(args):

    bvecs = np.loadtxt(args.in_vec)

    # read rotation matrix
    basename = os.path.splitext(args.in_mat)[0]
    transparams = np.loadtxt(basename + '.par')[:,:3]
    rotmat = []
    for param in transparams:
        rotmat.append(rot_mat(param))
    rotmat = np.asarray(rotmat)

    if len(bvecs) != len(rotmat):
        raise ValueError("Length of b-vector array does not match the length of the transformation matrix array")

    bvecs_rot = rotmat @ bvecs[:,:,np.newaxis]
    bvecs_rot = bvecs_rot[...,0]
    filename_rot = args.in_vec + '.rot'
    np.savetxt(filename_rot, bvecs_rot, fmt='%1.10f')

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Rotates b-vector acc to McFlirts rotation parameters.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('in_vec', type=str, help='Input bvec file (.bvec).')
    parser.add_argument('in_mat', type=str, help='Input McFlirts transformation parameters (.par).')

    args = parser.parse_args()

    main(args)
