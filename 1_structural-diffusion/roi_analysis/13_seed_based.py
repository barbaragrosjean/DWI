#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This file provide a way to do a seed based tractography from the full-brain tractograme previously obtained. 
# Must be used after the dwi pipeline once the full brain tractography is extracted.

from __future__ import division

import argparse
import os
import sys
import time
import itertools
import json
import subprocess
import logging
from datetime import datetime
import nibabel as nib
import numpy as np

sys.path.insert(1,'/home/bgrosjea/mnt/Hummel-Data/TI/mri/51T/barbara/uphummel_imaging_template/1_structural-diffusion')
from tools.tck2conn4stream_measures import * 
from tools.formate_data import formate_seed_based

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


def seed_based(data_path:str, subj:str, sess:str, isVerbose:bool, isForce:bool): 
    ''' 
    Track extraction using tckedit 
    '''
    
    tract_folder =  os.path.join(data_path, "derivatives", "01_tracts", subj, sess)
    session_folder = os.path.join(data_path, "derivatives", "01_dwi", subj, sess) 

    if not os.path.exists(tract_folder):
        os.makedirs(tract_folder)
    
    if os.path.exists(session_folder): 
        #------------------------------------------------------------#
        # SETUP:  TO CHANGE IF WANT TO RUN THE ANALYSIS FOR OTHER ROIS
        #------------------------------------------------------------#
        seed_folder = os.path.join(tract_folder, 'striat')
        rois = ['v_d_Ca_L', 'v_d_Ca_R', 'vm_dl_PU_L', 'vm_dl_PU_R']
        #------------------------------------------------------------#

        #------------------------------------------------------------#
        #### SEED-BASED : TRACTS EXTRACTION FOR ALL THE SEED #### 
        #------------------------------------------------------------#

        # All brain tractograme file
        dwi_out = os.path.join(session_folder, 'dwi', "proc", subj + "_" + sess)
        tck_file = dwi_out + "_iFOD2.tck"
        sift_file = dwi_out + "_sift.txt"
        
        for roi in rois: 
            print('Processing roi :', roi)

            # Select the right file 
            roi_file = os.path.join(seed_folder,  subj + "_" + sess + "_roi_" + str(roi) + "_dwi_ants.nii.gz")
                
            if not os.path.isfile(roi_file) : 
                logging.info('must perform the registration step (11_)')
                return 

            # Apply to the tractogram with tckedit command
            tck_out_path = os.path.join(seed_folder, 'tracts_tckedit')

            if not os.path.exists(tck_out_path):
                os.makedirs(tck_out_path)
            
            tck_out_file = os.path.join(tck_out_path, subj + "_" + sess + "_" + str(roi) + ".tck")
            json_out = os.path.join(tck_out_path, subj + "_" + sess + "_" + str(roi) + ".json")
            include_options =  "-include " + roi_file #"-mask " 

            if os.path.isfile(tck_out_file) and not isForce:
                print(f'%s file already existing' %tck_out_file)
            else:  
                sift_outpath = os.path.join(tck_out_path, subj + "_" + sess + "_" + str(roi) + "_sift2.txt")
                tckedit_cmd = "tckedit " + tck_file + " " + tck_out_file + " " + include_options + " -tck_weights_in " + sift_file + " -tck_weights_out " + sift_outpath + ' -force'
                
                logging.info('tckedit command: "{0}".'.format(tckedit_cmd))
                subprocess.call(tckedit_cmd, shell=True)

                with open(json_out, 'w') as outfile:
                    j = {
                        'Origin function': tckedit_cmd,
                        'Description': 'Applied rois selection to the tractogram with tckedit command',
                        'Anat_filename': tck_out_file,
                        'Time' : time.asctime()
                        }
                    json.dump(j, outfile)

            # tck2connectome command         
            # Ouputs 
            connectome = os.path.join(seed_folder, subj + "_" + sess +"_"+ roi +"_metric.csv")
            json_out = os.path.join(seed_folder,  subj + "_" + sess + "_"+ roi +"_metric.json")
            
            if os.path.isfile(connectome) and not isForce: 
                logging.info('connectom already done: "{0}".'.format(connectome))
            else:
                tck2connectome_cmd = 'tck2connectome ' + tck_file + ' ' + roi_file + ' ' + connectome + ' -tck_weights_in ' +  sift_file + ' -force' 

                logging.info('tck2connectome command: "{0}".'.format(tck2connectome_cmd))
                subprocess.call(tck2connectome_cmd, shell=True)
                
                with open(json_out, 'w') as outfile:
                    j = {
                        'Origin function': tck2connectome_cmd,
                        'Description': 'extract csv connectome',
                        'Anat_filename': connectome,
                        'Time' : time.asctime()
                        }
                    json.dump(j, outfile)

            

    else:
        raise FileNotFoundError("subj " + subj + ", sess " + sess + " not existing")
        return
      

if __name__ == "__main__":
    print('Starting roi-to-roi : extracting seed base tracts...')

    parser = buildArgsParser()
    args = parser.parse_args()

    isForce = args.isForce

    if args.isVerbose:
        logging.basicConfig(level=logging.DEBUG)
    
    subj_list = [subj for subj in args.subj]

    data_path = args.data_path

    if "all" in subj_list:
        subjects = [s for s in os.listdir(data_path) if os.path.isdir(os.path.join(data_path,'derivatives', '01_dwi', s)) if "sub-51T" in s]
    else:
        subjects = ['sub-' + subj for subj in subj_list]
    
    
    sess_list = [sess for sess in args.sess]
    if "all" in sess_list:
       sessions = ["ses-T1", "ses-T2", "ses-T3", "ses-T4"]
    else:
        sessions = ['ses-' + sess for sess in sess_list]
    
    date = datetime.now()
    formatted_datetime = date.strftime("%Y-%m-%d-%H-%M-%S")
    fail_list_filename = f"fail_list_13_seed_based{formatted_datetime}.txt"
    for subj, sess in itertools.product(subjects, sessions):        
            try:
                seed_based(data_path, subj, sess, args.isVerbose, isForce)
                
            except Exception as e:
                with open(fail_list_filename, "+a") as f:
                    f.write(f"{subj} {sess} \n")
                    f.write(f"{str(e)} \n")
    