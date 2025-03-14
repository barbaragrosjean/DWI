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

def lesion_reg(data_path:str, subj:str, sess:str, isForce:bool):
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
    
    if not os.path.exists(session_folder):
        os.makedirs(session_folder)

    # if os.path.exists(session_folder):
    #     logging.info('Processing dataset: "{0}".'.format(session_folder))

        anat_folder = os.path.join(data_path, subj, sess, "anat")
        dwi_folder = os.path.join(session_folder, "dwi")
        warp_folder = os.path.join(session_folder, "warps")
        lesion_transpl_folder = os.path.join(data_path, "derivatives", "03_dwi", "0_lesion_transplantations_FW", subj, sess, "anat", "lesion_transplantation")
        lesion_folder_dwi = os.path.join(data_path, "derivatives","07_lesions", subj, sess, "dwi")
        lesion_folder_mni = os.path.join(data_path, "derivatives","07_lesions", subj, sess, "mni")
        
        if not os.path.exists(lesion_folder_dwi):
            os.makedirs(lesion_folder_dwi)
            
        if not os.path.exists(lesion_folder_mni):
            os.makedirs(lesion_folder_mni)

        ##################
        ## Registration ##
        ##################
        t1_brain_filename = os.path.join(session_folder, "anat", subj + "_" + sess + "_acq-mprage_T1wBrain.nii.gz")
        meanB0bet_filename = os.path.join(dwi_folder, subj + "_" + sess + "_meanB0bet.nii.gz")
        T1w_transpl_file = os.path.join(lesion_transpl_folder, subj + "_" + sess + "_T1w_with_transplanted_lesion.nii.gz")
        MNI_file = "/usr/local/fsl/data/standard/MNI152_T1_1mm.nii.gz"

        lesion_list = ("acute", "combined", "old")
        for les in lesion_list:
            lesion_in = os.path.join(anat_folder,subj + "_" + sess + "_T1w_label-" + les + "lesion_roi.nii.gz")
            lesion_out_dwi = os.path.join(lesion_folder_dwi, subj + "_" + sess + "_T1w_label-" + les + "lesion_roi_dwi.nii.gz")
            lesion_out_mni = os.path.join(lesion_folder_mni, subj + "_" + sess + "_T1w_label-" + les + "lesion_roi_mni.nii.gz")

            # # lesion to dwi space (meanB0)
            #     # T1 brain file (BET) used here because DWI also has lesion
            # json_file = os.path.join(lesion_folder_dwi, subj + "_" + sess + "_T1w_label-" + les + "lesion_roi_dwi.json")
            # if os.path.isdir(os.path.dirname(lesion_in)) and (not os.path.isfile(lesion_out_dwi)) or isForce:
            #     warp_name = "T1w2meanB0_ants"
            #     original_file = t1_brain_filename
            #     ref_file = meanB0bet_filename
            #     inv = False
            #     input_file = lesion_in
            #     output_file = lesion_out_dwi
            #     interp_meth = "MultiLabel"
            #     registerAnts(input_file, output_file, warp_folder, warp_name, original_file, ref_file, inv, interp_meth)
            #     if os.path.isfile(output_file):
            #         with open(json_file, 'w') as outfile:
            #             j = {
            #                 'Origin function': "registerAnts with registration_ants.py" + input_file + "&" + output_file + " ",
            #                 'Description': 'Register lesion to DWI space',
            #                 'Anat_filename': lesion_in,
            #                 'Time' : time.asctime()
            #                 }
            #             json.dump(j, outfile)                 
            # else:
            #     with open(fail_list_filename, "+a") as failed_reg:
            #         mess = 'lesion already in dwi space or not existing: "{0}".\n'.format(lesion_out_dwi)
            #         failed_reg.write(mess)
            
            # lesion to mni space 
                # T1w tranplanted used here since lesion is not present in MNI and could lead to errors
            json_file = os.path.join(lesion_folder_mni, subj + "_" + sess + "_T1w_label-" + les + "lesion_roi_mni.json")
            if (os.path.isdir(os.path.dirname(lesion_in)) and (not os.path.isfile(lesion_out_mni))) or isForce:
                warp_name = "T1wtranspl2MNI_ants"
                original_file = T1w_transpl_file
                ref_file = MNI_file
                inv = False
                input_file = lesion_in
                output_file = lesion_out_mni
                interp_meth = "MultiLabel"
                registerAnts(input_file, output_file, warp_folder, warp_name, original_file, ref_file, inv, interp_meth)
                if os.path.isfile(output_file):
                    with open(json_file, 'w') as outfile:
                        j = {
                            'Origin function': "registerAnts with registration_ants.py" + input_file + "&" + output_file + " ",
                            'Description': 'Register lesion to MNI space',
                            'Anat_filename': lesion_in,
                            'Time' : time.asctime()
                            }
                        json.dump(j, outfile)
            else:
                with open(fail_list_filename, "+a") as failed_reg:
                    mess = 'lesion already in mni space or not existing: "{0}".\n'.format(lesion_out_mni)
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
    fail_list_filename = f"fail_list_04_lesion_registration_{formatted_datetime}.txt"
    for subj, sess in itertools.product(subjects, sessions):
        try:
            lesion_reg(data_path, subj, sess, isForce)
        except Exception as e:
            with open(fail_list_filename, "+a") as failed_reg:
                failed_reg.write(f"{str(e)} \n")


        # if (subj == "sub-TIMESwp11s017" and sess == "ses-T2") or \
        #     (subj == "sub-TIMESwp11s027") or \
        #     (subj == "sub-TIMESwp11s050" and sess == "ses-T1") or \
        #     (subj == "sub-TIMESwp11s050" and sess == "ses-T3") or \
        #     (subj == "sub-TIMESwp11s057"):