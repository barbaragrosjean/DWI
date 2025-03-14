#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division

import argparse
import logging
import os
import json
import subprocess
import nibabel as nib
import numpy as np
import os.path
import shutil
import time
from copy import copy
from genericpath import isfile
import itertools
from tools.registration_ants import *
from datetime import datetime


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
    return p 
  

def dwi_processing_func(data_path:str, subj:str, sess:str, isForce:bool):
    ''' Function doing the tractography on the dwi data.
    
        Parameters
        ----------
        data_path :
            Path containing folder with all data
        subj_list :
            Current subject
        sess_list :
            Current session
        isForce :
            Boolean indicating if files have to be overwritten
    '''
    
    session_folder = os.path.join(data_path, "derivatives","01_dwi", subj, sess)
    logging.info('Processing subject/session: {0}.'.format(subj + " / " + sess))

    # Input data are in preproc and output in proc
    preproc_folder = os.path.join(session_folder, "dwi", "preproc")
    proc_folder = os.path.join(session_folder, "dwi", "proc")
    
    #------------------------------------------------------------#
    #### 1 DWI2RESPONSE MSMT 5TT CMD #### 
    #------------------------------------------------------------#
   
    print('#### dwi2response ####')

    # Input
    tt5_file = os.path.join(preproc_folder,subj + "_" + sess + "_acq-mprage_T1wPve5tt_dwi.nii.gz")
    dwi_base_filename = os.path.join(preproc_folder,subj + "_" + sess + "_dwi")
    # Output 
    json_file = dwi_base_filename + "Resp.json"
    respWM_filename = dwi_base_filename + "RespWM.txt"
    respGM_filename = dwi_base_filename + "RespGM.txt"
    respCSF_filename = dwi_base_filename + "RespCSF.txt"

    # If output already existes
    if os.path.isfile(respWM_filename) and not isForce:
        logging.info('dwi2response already done.')
    else:         
        dwi2response_cmd = "dwi2response msmt_5tt " + dwi_base_filename + ".nii.gz " + tt5_file + " " + respWM_filename + " " + respGM_filename + " " + respCSF_filename + " -fslgrad " + dwi_base_filename + ".bvec " + dwi_base_filename + ".bval -nthreads 8"
        logging.info('dwi2response command.')
        subprocess.call(dwi2response_cmd, shell=True)
        with open(json_file, 'w') as outfile:
            j = {
                'Origin function': dwi2response_cmd,
                'Description': 'DWI 2 response',
                'respWM_filename': respWM_filename,
                'Time' : time.asctime()
                }
            json.dump(j, outfile)        

    
    #------------------------------------------------------------#
    #### 2 DWI2FOD CMD #### 
    #------------------------------------------------------------#

    print('#### dwi2fod ####')
    # Input
    t1_mask_filename = os.path.join(preproc_folder,subj + "_" + sess + "_acq-mprage_T1w_dwi.nii.gz")
    # Outputs
    fodWM_filename = os.path.join(proc_folder,subj + "_" + sess + "_fod.nii.gz")
    fodGM_filename = os.path.join(proc_folder,subj + "_" + sess + "_fodGM.nii.gz")
    fodCSF_filename = os.path.join(proc_folder,subj + "_" + sess + "_fodCSF.nii.gz")
    json_file = os.path.join(proc_folder,subj + "_" + sess + "_fod.json")

    if os.path.isfile(fodWM_filename) and not isForce:
        logging.info('dwi2fod already done,')
    else:
        dwi2fod_cmd = "dwi2fod msmt_csd -mask " + t1_mask_filename + " " + dwi_base_filename + ".nii.gz " + respWM_filename + " " + fodWM_filename + " " + respGM_filename + " " + fodGM_filename + " " + respCSF_filename + " " + fodCSF_filename + " " + " -fslgrad " + dwi_base_filename + ".bvec " + dwi_base_filename + ".bval -nthreads 8"
        logging.info('dwi2fod command: "{0}".'.format(dwi2fod_cmd))
        subprocess.call(dwi2fod_cmd, shell=True)
        with open(json_file, 'w') as outfile:
            j = {
                'Origin function': dwi2fod_cmd,
                'Description': 'DWI 2 FOD',
                'fod filename': fodWM_filename,
                'Time' : time.asctime()
                }
            json.dump(j, outfile)

    #------------------------------------------------------------#
    #### 3 STREAMLINE COUNT AND REMOVE UNRELEVANT CONNECTIONS #### 
    #------------------------------------------------------------#
   
    print('#### Streamline count  ####')
    streamlines_count=10000000
    # Input
    wm_pve_filename = os.path.join(preproc_folder,subj + "_" + sess + "_acq-mprage_T1wPveWM_dwi.nii.gz")
    # Outputs
    streamlines_filename = os.path.join(proc_folder,subj + "_" + sess + "_iFOD2.tck")
    json_file = os.path.join(proc_folder,subj + "_" + sess + "_iFOD2.json")

    if os.path.isfile(streamlines_filename) and not isForce:
        logging.info('tckgen iFOD2 already done.')
    else:
        tckgen_cmd = "tckgen " + fodWM_filename + " " + streamlines_filename + \
            " -algorithm iFOD2 -seed_image " + wm_pve_filename + " -select " + \
            str(streamlines_count) + " -force -minlength 1.6 -nthreads 8"
        logging.info('tckgen command: "{0}".'.format(tckgen_cmd))
        subprocess.call(tckgen_cmd, shell=True)
        with open(json_file, 'w') as outfile:
            j = {
                'Origin function': tckgen_cmd,
                'Description': 'Generate tck file',
                'tck filename': streamlines_filename,
                'Time' : time.asctime()
                }
            json.dump(j, outfile)   

    #------------------------------------------------------------#
    #### 4 EXTRACTING WEIGHTS WITH SIFT2 CMD #### 
    #------------------------------------------------------------#

    json_file = os.path.join(proc_folder,subj + "_" + sess + "_sift.json")
    tck_sift_file = os.path.join(proc_folder,subj + "_" + sess + "_sift.txt")
    wm_mask_file = os.path.join(preproc_folder,subj + "_" + sess + "_acq-mprage_T1wPveWM_dwi.nii.gz")

    if os.path.isfile(tck_sift_file) and not isForce:
        logging.info('tcksift2 already done.')
    else: 
        tcksift2_cmd = "tcksift2 -act " + tt5_file + " " + streamlines_filename + " " + fodWM_filename + " " + tck_sift_file + " -proc_mask " + wm_mask_file + " -force " 

        logging.info('tcksift2 command: "{0}".'.format(tcksift2_cmd))

        subprocess.call(tcksift2_cmd, shell=True)

        with open(json_file, 'w') as outfile:
            j = {
                'Origin function': tcksift2_cmd,
                'Description': 'Generate tcksift2 file',
                'tck filename': tck_sift_file,
                'Time' : time.asctime()
                }
            json.dump(j, outfile)    
    return


if __name__ == "__main__":
    
    parser = buildArgsParser()
    args = parser.parse_args()

    isForce = args.isForce
    if args.isVerbose:
        logging.basicConfig(level=logging.DEBUG)

    subj_list = [subj for subj in args.subj]
    sess_list = [sess for sess in args.sess]

    data_path = args.data_path

    if "all" in subj_list:
        subjects = [s for s in os.listdir(data_path) if os.path.isdir(os.path.join(data_path, s)) if "sub-51T" in s]
    else:
        subjects = ['sub-' + subj for subj in subj_list]

    if "all" in sess_list:
        sessions = ["ses-T1", "ses-T2", "ses-T3", "ses-T4"]
    else:
        sessions = ['ses-' + sess for sess in sess_list]

    date = datetime.now()
    formatted_datetime = date.strftime("%Y-%m-%d-%H-%M-%S")
    fail_list_filename = f"fail_list_06_dwi_processing_{formatted_datetime}.txt"
    for subj, sess in itertools.product(subjects, sessions):        
            try:
                dwi_processing_func(data_path, subj, sess, isForce)
            except Exception as e:
                with open(fail_list_filename, "+a") as f:
                    f.write(f"{subj} {sess} \n")
                    f.write(f"{str(e)} \n")    
    
