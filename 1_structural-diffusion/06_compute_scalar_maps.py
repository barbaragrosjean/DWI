#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import logging
import os
import sys
import json
import subprocess
import time
import nibabel as nib
import numpy as np
from datetime import datetime

import itertools

sys.path.append('/home/bgrosjea/mnt/Hummel-Data/TI/mri/51T/barbara/uphummel_imaging_template/1_structural-diffusion')
from tools.registration_ants import *
#from tools.lesionTransplantation_native import * # if lesion not tried 


def buildArgsParser():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter,
        epilog="")
    p._optionals.title = "Generic options"

    p.add_argument('--subj', nargs='+', dest='subj', help="Subject index.")
    p.add_argument('--sess', nargs='+', dest='sess', help="Session folder name.")

    p.add_argument('--data_path', default='/mnt/Hummel-Data/TI/mri/51T', dest='data_path',
        help="Subjects folder path. ['%(default)s']")

    p.add_argument('-f', action='store_true', dest='isForce',
    help='If set, overwrites output file.')


    log_g = p.add_argument_group('Logging options')
    log_g.add_argument(
        '-v', action='store_true', dest='isVerbose',
        help='If set, produces verbose output.')
    
    p.add_argument('--l', action='store_true', dest='lesion', help='True if lesion, default = False')

    return p 
  

def scalar_maps_fnct(data_path:str, subj:str, sess:str, isForce:bool, lesion:bool, isVerbose:bool):
    """Compute scalar maps including FA and MD maps"""
    
    session_folder = os.path.join(data_path, 'derivatives','01_dwi',subj, sess)
    anat_folder = os.path.join(session_folder, "anat")
    
    if lesion:
        lesion_folder = os.path.join(session_folder, "lesion")
        if not os.path.exists(lesion_folder):
            os.makedirs(lesion_folder)
            
    dwi_out_preproc = os.path.join(session_folder, "dwi", "preproc", subj + "_" + sess + "_dwi")
    dwi_out_proc = os.path.join(session_folder, "dwi", "proc", subj + "_" + sess + "_dwi")

    dwi_file = dwi_out_preproc + ".nii.gz"
    #------------------------------------------------------------#
    #### 1 FITS DIFFUSITON TENSOR MODEL AT EACH VOXEL #### 
    #------------------------------------------------------------#

    print('#### Fits diffusion tensor model at each voxel ####')
    bvec_dwi_file = dwi_out_preproc + ".bvec"
    bval_dwi_file = dwi_out_preproc + ".bval"
    FAMaps_file = dwi_out_proc + "_FA.nii.gz"
    bet_file = dwi_out_preproc + "_mean-b0_bet.nii.gz"
    
    if os.path.isfile(FAMaps_file) and not isForce:
        print("dtifit already performed on Subject")
    else:
        dtifit_cmd = "dtifit --data=" + dwi_file + " --out=" + dwi_out_proc + " --mask=" + bet_file + \
                     " --bvecs=" + bvec_dwi_file + " --bvals=" + bval_dwi_file
        print(dtifit_cmd)
        logging.info('dtifit command: "{0}".'.format(dtifit_cmd))
        subprocess.call(dtifit_cmd, shell=True)
        with open(dwi_out_proc + "_FA.json", 'w') as outfile:
            j = {
                "dtifit" : dtifit_cmd,
                'Time' : time.asctime()
                }
            json.dump(j, outfile)
    #------------------------------------------------------------#
    #### 2 REGISTER LESION INTO DWI SPACE #### 
    #------------------------------------------------------------#
    
    if lesion:
        print('#### Register lesion into dwi space ####')
        os.chdir(session_folder)
        lesion_anat_file = os.path.join(anat_folder, subj + "_" + sess + "_acq-mprage_T1w_label-lesion_roi.nii.gz")

        les_dwi_file = os.path.join(lesion_folder, subj + "_" + sess + "_acq-mprage_T1w_label-lesion_roi_dwi.nii.gz")

        warp_folder = os.path.join(output_session_folder, "warps")
        if not os.path.exists(warp_folder):
            os.makedirs(warp_folder)

        input_file = lesion_anat_file
        output_file = les_dwi_file
        json_file = os.path.join(lesion_folder, subj + "_" + sess + "_acq-mprage_T1w_label-lesion_roi_dwi.json")
        warp_name = 'T1w2mean_b0'
        original_file_name = os.path.join(anat_folder, subj + "_" + sess + "_acq-mprage_T1w_bet.nii.gz")
        ref_file = bet_file
        warp_name = 'T1w_bet2mean_b0_bet_ants'

        interp_bool = "NearestNeighbor"
        origin_fnct = 'registerAnts(' + str(input_file) + ', ' + str(output_file) + ', ' + \
                    str(warp_folder) + ', ' + str(warp_name) + ', ' + str(original_file_name) + \
                    ', ' + str(ref_file) + ', ' + str(interp_bool) + ')'
        print(origin_fnct)
        registerAnts(input_file, output_file, warp_folder, warp_name, original_file_name, ref_file, interp_bool)

        with open(json_file, 'w') as outfile:
            j = {
                'Origin function': origin_fnct,
                'Description': 'Register lesion from anatomical to dwi space',
                'DWI_filename': les_dwi_file,
                'Time' : time.asctime()
                }
            json.dump(j, outfile)

    return


if __name__ == "__main__":
    parser = buildArgsParser()
    args = parser.parse_args()

    isForce = args.isForce
    lesion = args.lesion
    
    if args.isVerbose:
        logging.basicConfig(level=logging.DEBUG)
    
    subj_list = [subj for subj in args.subj]

    data_path = args.data_path

    if "all" in subj_list:
        subjects = [s for s in os.listdir(data_path) if os.path.isdir(os.path.join(data_path, s)) if "sub-51T" in s]
    else:
        subjects = ['sub-' + subj for subj in subj_list]
    
    
    sess_list = [sess for sess in args.sess]
    if "all" in sess_list:
       sessions = ["ses-T1", "ses-T2", "ses-T3", "ses-T4"]
    else:
        sessions = ['ses-' + sess for sess in sess_list]
    
    date = datetime.now()
    formatted_datetime = date.strftime("%Y-%m-%d-%H-%M-%S")
    fail_list_filename = f"fail_list_06_compute_scalar_maps{formatted_datetime}.txt"
    for subj, sess in itertools.product(subjects, sessions):        
            try:
                scalar_maps_fnct(data_path, subj, sess, isForce, lesion, args.isVerbose)
            except Exception as e:
                with open(fail_list_filename, "+a") as f:
                    f.write(f"{subj} {sess} \n")
                    f.write(f"{str(e)} \n")
