#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This file provide a way to extract specific tract from the all-brain tractography. 
# Must be used after the dwi pipeline once the full tractograme is extracted,  as a second step of the roi-to-roi analysis.

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
import pandas as pd

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

def create_parc(subj:str, sess:str,data_path:str, isForce:bool):  
    ''' 
        Create individual mask AND global mask for connectivity analysis.
    ''' 

    #------------------------------------------------------------#
    # SETUP:  TO CHANGE IF WANT TO RUN THE ANALYSIS FOR OTHER ROIS
    #------------------------------------------------------------#
    study = 'fMRI_study'
    rois = ['Loc_NA_Postcentral_L', 'Loc_NA_Cerebellum', 'Thal_IL_R', 'Precentral_L', 'Supp_Motor_Area_R']
    index = [2, 4, 5, 6, 7] # from cluster files
    
    roiClusters_file = os.path.join(data_path, "derivatives", "01_tracts", subj, sess, 'roi2roi', study, subj + "_" + sess +'_roi_Clusters_dwi_ants.nii.gz')
    striat = ['v_d_Ca_L', 'v_d_Ca_R', 'vm_dl_PU_L', 'vm_dl_PU_R']
    #------------------------------------------------------------#
    
    tract_folder =  os.path.join(data_path, "derivatives", "01_tracts", subj, sess)
    roi_folder = os.path.join(tract_folder, 'roi2roi')

    #------------------------------------------------------------#
    #### 1 CREATE INDIVIDUAL MASK FOR EACH ROI - FOR TRACT EXTRACTION #### 
    #------------------------------------------------------------#

    for idx_roi, roi in zip(index, rois) : 
        # Apply fslmath with threashold [v1, v2] and save the result into masks
        if not os.path.exists(os.path.join(roi_folder, study, 'masks')) : 
            os.makedirs(os.path.join(roi_folder, study, 'masks'))

        mask_file = os.path.join(roi_folder, study, 'masks', subj + '_' + sess + "_roi_" + str(roi) +'_mask.nii.gz')
        json_out = os.path.join(roi_folder, study, 'masks', subj + '_' + sess + "_roi_" + str(roi) +'_mask.json')

        # set the tresholds 
        v1 = idx_roi - 0.1
        v2 = idx_roi + 0.1

        if os.path.isfile(mask_file) and not isForce: 
            logging.info('Individual mask already done: "{0}".'.format(mask_file))
            print('Individual mask already done: "{0}".'.format(mask_file))
        else:
            fslmath_cmd_thr = "fslmaths " + roiClusters_file + " -thr " + str(v1) + " -uthr " + str(v2) + " -bin " + mask_file 
            logging.info('fslmaths threshold command: "{0}".'.format(fslmath_cmd_thr))
            subprocess.call(fslmath_cmd_thr, shell=True)
                
            with open(json_out, 'w') as outfile:
                j = {
                    'Origin function': fslmath_cmd_thr,
                    'Description': 'Apply fslmath with threashold [' + str(v1) + ',' + str(v2) +']',
                    'Anat_filename': mask_file,
                    'Time' : time.asctime()
                    }
                json.dump(j, outfile)
          

    # add the striatum to the label list
    tot_rois = striat + rois

    #------------------------------------------------------------#
    #### 2 CREATE GLOBAL MASK FOR WITH ALL THE ROI - FOR CONNECTOME #### 
    #------------------------------------------------------------#

    parcellation_file = os.path.join(roi_folder, study, 'masks', subj + "_" + sess + "_global_mask.nii.gz")
    json_file = os.path.join(roi_folder, study, 'masks', subj + "_" + sess + "_global_mask.json")

    if os.path.isfile(parcellation_file) and not isForce:
        logging.info('global mask already done: "{0}".'.format(parcellation_file))
        print('global mask already done: "{0}".'.format(parcellation_file))
    else : 
        parc = nib.load(roiClusters_file)
        parc_data = parc.get_fdata()
        reduced = np.zeros(parc_data.shape)

        for i, roi in enumerate(tot_rois): 
            # Select the right individual mask
            if  roi in striat : 
                mask_file = os.path.join(tract_folder, 'striat',  subj + "_" + sess + "_roi_" + roi + "_dwi_ants.nii.gz")
            else: 
                mask_file = os.path.join(roi_folder, study, 'masks', subj + '_' + sess + "_roi_" + str(roi) +'_mask.nii.gz')
            
            mask = nib.load(mask_file)
            reduced[mask.get_fdata() == 1] = i + 1


        nib.Nifti1Image(reduced, parc.affine, parc.header).to_filename(parcellation_file)

        with open(json_file, 'w') as outfile:
            j = {
                'Origin function': './processing_scripts/diffusion/ROIs_functions/2_create_parc.py',
                'Description': 'Parcellation with ROIs and Striat',
                'DWI_filename': parcellation_file,
                'Time': time.asctime()
                }
            json.dump(j, outfile)
                    
        # Save the labels into a csv file
        label_file = os.path.join(roi_folder, study, 'masks', subj + "_" + sess + "_global_mask.csv")
        lab = {'roi' : tot_rois}
        df = pd.DataFrame(lab, index=(np.arange(1,len(tot_rois)+1)))
        df.to_csv(label_file)


if __name__ == "__main__":
    print('Starting roi-to-roi : extracting roi masks...')
    
    parser = buildArgsParser()
    args = parser.parse_args()

    isForce = args.isForce
        
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
    fail_list_filename = f"fail_list_12_create_parc{formatted_datetime}.txt"
    for subj, sess in itertools.product(subjects, sessions):        
            try:
                create_parc(subj, sess,data_path, isForce)
            except Exception as e:
                with open(fail_list_filename, "+a") as f:
                    f.write(f"{subj} {sess} \n")
                    f.write(f"{str(e)} \n")
