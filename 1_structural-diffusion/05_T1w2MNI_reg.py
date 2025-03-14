#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division

import argparse
import logging
import os
import json
import time
import nibabel as nib
import numpy as np
import os.path
import itertools
from joblib import Parallel, delayed
from tools.registration_ants import *
from datetime import datetime


def buildArgsParser():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter,
        epilog="")
    p._optionals.title = "Generic options"

    p.add_argument('--subj', nargs='+', dest='subj', help="Subject index.")
    p.add_argument('--sess', nargs='+', dest='sess', help="Session folder name.")

    p.add_argument('--data_path', default='/media/windel/Elements/TiMeS/WP11_MRI/data/times_wp11', dest='data_path',
        help="Subjects folder path. ['%(default)s']")

    p.add_argument('-f', action='store_true', dest='isForce',
    help='If set, overwrites output file.')

    log_g = p.add_argument_group('Logging options')
    log_g.add_argument(
        '-v', action='store_true', dest='isVerbose',
        help='If set, produces verbose output.')
    return p

def T1_reg(data_path:str, subj:str, sess:str, isForce:bool):
    ''' Function computing registration from T1w space to dwi space
    
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

    session_folder = os.path.join(data_path, "derivatives","03_dwi","2_dwi_processing_FW", subj, sess)
    mni_folder = os.path.join(data_path, "derivatives","04_mni", subj, sess)

    if os.path.exists(session_folder):
        logging.info('Processing dataset: "{0}".'.format(session_folder))

        anat_folder = os.path.join(data_path, subj, sess, "anat")
        warp_folder = os.path.join(session_folder, "warps")
        lesion_transpl_folder = os.path.join(data_path, "derivatives", "03_dwi", "0_lesion_transplantations_FW", subj, sess, "anat", "lesion_transplantation")

        
        if not os.path.exists(mni_folder):
            os.makedirs(mni_folder)

        ##################
        ## Registration ##
        ##################
        t1w_filename = os.path.join(anat_folder, subj + "_" + sess + "_acq-mprage_T1w.nii.gz")
        t1w_out_filename = os.path.join(mni_folder, subj + "_" + sess + "_acq-mprage_T1w_mni.nii.gz")
        T1w_transpl_file = os.path.join(lesion_transpl_folder, subj + "_" + sess + "_T1w_with_transplanted_lesion.nii.gz")
        MNI_file = "/usr/local/fsl/data/standard/MNI152_T1_1mm.nii.gz"

        # T1w to mni space 
            # T1w tranplanted used here since lesion is not present in MNI and could lead to errors
        json_file = os.path.join(mni_folder, subj + "_" + sess + "_acq-mprage_T1w_mni.json")
        if (not os.path.isfile(t1w_out_filename)) or isForce:
            warp_name = "T1wtranspl2MNI_ants"
            original_file = T1w_transpl_file
            ref_file = MNI_file
            inv = False

            input_file = t1w_filename
            output_file = t1w_out_filename
            registerAnts(input_file, output_file, warp_folder, warp_name, original_file, ref_file, inv)
            if os.path.isfile(output_file):
                with open(json_file, 'w') as outfile:
                    j = {
                        'Origin function': "registerAnts with registration_ants.py" + input_file + "&" + output_file + " ",
                        'Description': 'Register lesion to MNI space',
                        'mni_filename': output_file,
                        'Time' : time.asctime()
                        }
                    json.dump(j, outfile)
        else:
            with open(fail_list_filename, "+a") as failed_reg:
                mess = 'lesion already in mni space or not existing: "{0}".\n'.format(t1w_out_filename)
                failed_reg.write(mess)

    else:
        raise FileNotFoundError("subj " + subj + ", sess " + sess + " not existing")

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
        subjects = [s for s in os.listdir(data_path) if os.path.isdir(os.path.join(data_path, s)) if "sub-TIMESwp11s" in s]
    else:
        subjects = ['sub-' + subj for subj in subj_list]

    if "all" in sess_list:
        sessions = ["ses-T1", "ses-T2", "ses-T3", "ses-T4"]
    else:
        sessions = ['ses-' + sess for sess in sess_list]

    #Parallel(n_jobs=4)(delayed(anat_reg_dwi)(data_path, subj, sess, isForce)
    #                   for subj, sess in itertools.product(subjects, sessions))
    
    #for subj, sess in itertools.product(subjects, sessions):
           #lesion_reg(data_path, subj, sess, isForce)
        
    date = datetime.now()
    formatted_datetime = date.strftime("%Y-%m-%d-%H-%M-%S")
    fail_list_filename = f"fail_list_T1w2MNI_reg_{formatted_datetime}.txt"
    for subj, sess in itertools.product(subjects, sessions):
        try:
            T1_reg(data_path, subj, sess, isForce)
        except Exception as e:
            with open(fail_list_filename, "+a") as failed_reg:
                failed_reg.write(f"{str(e)} \n")