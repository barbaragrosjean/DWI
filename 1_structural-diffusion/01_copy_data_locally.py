#!/usr/bin/env python
# -*- coding: utf-8 -*-
# copies all files needed for processing pipeline from server to local computer

from __future__ import division

import argparse
from copy import copy
from genericpath import isfile
import logging
import os
import shutil
import itertools
from datetime import datetime

# Build Parser
def buildArgsParser():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter,
        epilog="")
    p._optionals.title = "Generic options"

    p.add_argument('--subj', nargs='+', dest='subj', help="Subject index.")
    p.add_argument('--sess', nargs='+', dest='sess', help="Session folder name.")
    
    p.add_argument('--data_path', default='/mnt/Hummel-Data/TiMeS/WP11_MRI/data/times_wp11', dest='data_path',
                   help="Data folder path. ['%(default)s']")
    p.add_argument('--out_path', default='/media/windel/Elements/TiMeS/WP11_MRI/data/times_wp11', dest='out_path',
                   help="Output folder path. ['%(default)s']")
    # p.add_argument('--out_path', default='/media/elena/EXT/data/subjects/patients/TiMeS11_EB', dest='out_path',
    #                help="Output folder path. ['%(default)s']")
    p.add_argument('-f', action='store_true', dest='isForce',
    help='If set, overwrites output file.')

    log_g = p.add_argument_group('Logging options')
    log_g.add_argument(
        '-v', action='store_true', dest='isVerbose',
        help='If set, produces verbose output.')
    return p

def transfer_local(data_path:str, out_path:str, subj:str, sess:str, isForce:bool):
    ''' Copy files from original preprocessing path to modeling path to 
    create a folder with all rquired data for modeling
    
        Parameters
        ----------
        data_path :
            Path containing folder with all data
        out_path :
            Path containing folder where to place data
        subj_list :
            Current subject
        sess_list :
            Current session
        isForce :
            Boolean indicating if files have to be overwritten
    '''
    
    # Define server folders & file prefixes
    anat_raw_folder = os.path.join(data_path, subj, sess, "anat")
    anat_raw_prefix = os.path.join(anat_raw_folder, subj + '_' + sess)

    dwi_raw_folder = os.path.join(data_path, subj, sess, "dwi")
    dwi_raw_prefix = os.path.join(dwi_raw_folder, subj + '_' + sess)

    # Define local folders & file prefixes
    anat_out = os.path.join(out_path, subj, sess, "anat")
    dwi_out = os.path.join(out_path, 'derivatives', '03_dwi', '2_dwi_processing_FW', subj, sess, "dwi")
    dwi_raw_out = os.path.join(out_path, subj, sess, "dwi")
    
    # copy T1w anat raw
    files_to_copy = (anat_raw_prefix + '_acq-mprage_T1w.nii.gz', \
        anat_raw_prefix + '_acq-mprage_T1w.json', \
        anat_raw_prefix + '_T1w_label-acutelesion_roi.nii.gz', \
        anat_raw_prefix + '_T1w_label-oldlesion_roi.nii.gz', \
        anat_raw_prefix + '_T1w_label-combinedlesion_roi.nii.gz')
    
    # For each of the files you want to copy
    for orig_file in files_to_copy:
        # Does the file you want to copy exist?
        if (isfile(orig_file)):
            # If yes, we proceed

            # Does the folder where you want to copy the file exist?
            if not os.path.exists(anat_out):
                # If not, then create it
                os.makedirs(anat_out)
            
            # Once you are ready with the existing file and folder

            # Check if the file is not copied yet or if you want to force the copy anyway
            if (not isfile(os.path.join(anat_out, os.path.basename(orig_file)))) | (isForce):
                
                # If one of the two is true (not existing or force), display what I am copying if verbose mode is on 
                logging.info('Copying: "{0}/{1}/{2}".'.format(subj, sess, os.path.basename(orig_file)))

                # And then copy
                shutil.copy(orig_file, anat_out)

    
   # copy raw dwi
    files_to_copy = (dwi_raw_prefix + '_dir-AP_dwi.nii.gz', 
        dwi_raw_prefix + '_dir-AP_dwi.bval', 
        dwi_raw_prefix + '_dir-AP_dwi.bvec',
        dwi_raw_prefix + '_dir-AP_dwi.json',
        dwi_raw_prefix + '_dir-PA_dwi.bval', 
        dwi_raw_prefix + '_dir-PA_dwi.bvec',
        dwi_raw_prefix + '_dir-PA_dwi.json',          
        dwi_raw_prefix + '_dir-PA_dwi.nii.gz')

    for orig_file in files_to_copy:
        if (isfile(orig_file)):
            if not os.path.exists(dwi_raw_out):
                os.makedirs(dwi_raw_out)

            if (not isfile(os.path.join(dwi_raw_out, os.path.basename(orig_file)))) | (isForce):
                logging.info('Copying: "{0} / {1} / {2}".'.format(subj, sess, os.path.basename(orig_file)))
                shutil.copy(orig_file, dwi_raw_out)


if __name__ == "__main__":
    
    parser = buildArgsParser()
    args = parser.parse_args()

    isForce = args.isForce
    if args.isVerbose:
        logging.basicConfig(level=logging.DEBUG)

    data_path = args.data_path
    out_path = args.out_path

    subj_list = [subj for subj in args.subj]
    sess_list = [sess for sess in args.sess]
    
    if "all" in subj_list:
        subjects = [s for s in os.listdir(data_path) if os.path.isdir(os.path.join(data_path, s)) if "sub-TIMESwp11s" in s]
    else:
        subjects = ['sub-' + subj for subj in subj_list]

    if "all" in sess_list:
        sessions = ["ses-T1", "ses-T2", "ses-T3", "ses-T4"]
    else:
        sessions = ['ses-' + sess for sess in sess_list]
    
    for subj, sess in itertools.product(subjects, sessions):
            transfer_local(data_path, out_path, subj, sess, isForce)

    # date = datetime.now()
    # formatted_datetime = date.strftime("%Y-%m-%d-%H-%M-%S")
    # fail_list_filename = f"fail_list_01_copy_data_locally_{formatted_datetime}.txt"
    # for subj, sess in itertools.product(subjects, sessions):
    #     try:
    #         transfer_local(data_path, out_path, subj, sess, isForce)
    #     except Exception as e:
    #         with open(fail_list_filename, "+a") as f:
    #             f.write(f"{subj} {sess} \n")
    #             f.write(f"{str(e)} \n")

