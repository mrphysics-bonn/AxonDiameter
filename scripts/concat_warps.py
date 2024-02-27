#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Concatenate and apply warp fields
"""

import numpy as np
import subprocess
import os
from multiprocessing import Pool
import psutil
import argparse

def run_subprocess(cmd):
    subprocess.run(cmd, check=True)
    return True


def main(args):

    # input parameters
    dfield_dir = os.path.join(args.eddy_dir, "dfields")
    nlgc_file = os.path.join(args.unwarp_dir, "fullWarp_abs.nii.gz")
    flirt_mat = args.flirt_mat
    dwi = os.path.join(args.eddy_dir, "eddy_outlier_free_data.nii.gz")
    bval = args.bval
    out = args.out_file

    # make directories for results
    split_dir = os.path.join(args.eddy_dir, "split_data")
    comb_warp_dir = os.path.join(args.eddy_dir, "comb_warp")
    warp_dir = os.path.join(args.eddy_dir, "warped_data")

    os.makedirs(split_dir, exist_ok=True)
    os.makedirs(comb_warp_dir, exist_ok=True)
    os.makedirs(warp_dir, exist_ok=True)

    # Split DWI
    cmd = ["fslsplit", dwi, split_dir+"/data", "-t"]
    res = run_subprocess(cmd)
    dwi_files = sorted([os.path.join(split_dir, item) for item in os.listdir(split_dir)])

    # Read displacement fields from eddy
    dfield_files = sorted([os.path.join(dfield_dir, item) for item in os.listdir(dfield_dir)])

    # Read b-values
    if bval is None:
        bval = np.zeros(len(dwi_files))
    else:
        bval = np.loadtxt(bval)

    # Concatenate & apply warps
    cmd_list_convert = []
    cmd_list_apply = []
    cmd_list_avg_jac = []
    cmd_list_apply_jac = []
    for i,file in enumerate(dwi_files):
        comb_warp_file = comb_warp_dir + f"/comb_warp.{i}"
        jac_file = comb_warp_dir + f"/jac.{i}"
        warped_file = os.path.join(warp_dir, f"warped_vol.{i}")
        cmd = ["convertwarp",
            "-o", comb_warp_file, 
            "-r", file,
            "-j", jac_file,
            f"--warp1={dfield_files[i]}",
            f"--warp2={nlgc_file}"]
        if flirt_mat is not None and bval[i] > 20000:
            cmd.append(f"--postmat={flirt_mat}")
        cmd_list_convert.append(cmd)
        cmd = ["applywarp", 
            "-i", file, 
            "-r", file, 
            "-o", warped_file, 
            "-w", comb_warp_file,
            "--interp=spline"]
        cmd_list_apply.append(cmd)
        cmd = ["fslmaths",
               jac_file,
               "-Tmean",
               jac_file]
        cmd_list_avg_jac.append(cmd)
        cmd = ["fslmaths",
               warped_file,
               "-mul",
               jac_file,
               warped_file]
        cmd_list_apply_jac.append(cmd)

    cores = psutil.cpu_count(logical = False)
    pool = Pool(processes=cores)
    res = [pool.apply_async(run_subprocess, [cmd]) for cmd in cmd_list_convert]
    res2 = [val.get() for val in res]
    res = [pool.apply_async(run_subprocess, [cmd]) for cmd in cmd_list_apply]
    res2 = [val.get() for val in res]
    res = [pool.apply_async(run_subprocess, [cmd]) for cmd in cmd_list_avg_jac]
    res2 = [val.get() for val in res]
    res = [pool.apply_async(run_subprocess, [cmd]) for cmd in cmd_list_apply_jac]
    res2 = [val.get() for val in res]
    pool.close()

    # Merge DWI
    cmd = ["fslmerge", "-t", out] \
        + [os.path.join(warp_dir, f"warped_vol.{i}") for i in range(len(dwi_files))]
    res = run_subprocess(cmd)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Concatenate warp fields from eddy, gradient nonlinearity correction and optionally a FLIRT transform.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('eddy_dir', type=str, help='Input directory containing results from eddy (eddy_outlier_free_data.nii.gz and displacement fields in folder "dfields").')
    parser.add_argument('unwarp_dir', type=str, help='Input directory containing results from gradient nonlinearity correction (fullWarp_abs.nii.gz).')
    parser.add_argument('out_file', type=str, help='Output file for warped DWIs.')
    parser.add_argument('-f', '--flirt_mat', type=str, default=None, help='Optional transformation from FLIRT applied after eddy and gradient nonlinearity correction warps.')
    parser.add_argument('-b', '--bval', type=str, default=None, help='Optional b-value file as the FLIRT transform is only applied to the b=30000 shell.')

    args = parser.parse_args()

    main(args)
