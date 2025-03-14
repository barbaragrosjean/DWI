
#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This script is using T1 (the transplanted lesion or not) and applys freesurfer's recon all function = all parts of the cortical reconstruction process
# Reconstructing a two dimensional cortical surface from a three dimensional volume (tranplanted T1 in this case)
# General freesurfer (fs): any 3D volumes are stored in "mri" folder
# First step: skull strip to create brainmask.mgz (specific fs extension)
# Second step: estimation of interface between white and grey matter for the two hemispheres > saved as lh.orig and rh.orig
# Refined and saved under lh.white and rh.white
# Edge of grey matter detected and saved under lh.pial and rh.pial (all if the above can be visualised with freeview)
# The pial files can be expanded (to better determine locations along the banks and ridges of the gyri) and are saved under lh.inflated and rh.inflated
# Inflated surfaces can be inflated again into a sphere, normalised to a template image (fsaverage, 40 subjects)
# Once individual surface map is normalised to the template, two atlases can be used for parcellation: Desikan-Killiany and Destrieux (more detailed)
# The -all flag instructs fs to run all processing steps 
# -openmp 12 specifies the number of OpenMP (multi-platform shared-memory multiprocessing programming) threads to be used -> here 12 to speed up the processing time
# -brainstem-structures tells fs to include segmentation of braistem structures in the processing pipeline

#EDIT 31.07.2024 major fix: 
#many output files in the surf folder (subjxxx_ses-Tx/surf) were damaged because the a folder was missing in the freesurfer directory
#before running this script, the fsaverage folder from the root freesurfer folder needs to be copied into the output directory: 
#i.e. copy the whole folder "fsaverage" from /usr/local/freesurfer/subjects/fsaverage to the output folder of this script /media/windel/Elements/TiMeS/WP11_MRI/data/times_wp11/derivatives/03_dwi/1_freesurfer_FW/
#go to this website: https://figshare.com/articles/dataset/HCP-MMP1_0_projected_on_fsaverage/3498446, download "lh.hcpmmp1.annot" & "rh.hcpmmp1.annot" and "hcpmmp1_ordered.txt" from the git
#and add them in the label folder i.e. "/media/windel/Elements/TiMeS/WP11_MRI/data/times_wp11/derivatives/03_dwi/1_freesurfer_FW/fsaverage/label/"


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

    p.add_argument('--l', action='store_false', dest='lesion', help='True if lesion, default = False')

    return p

def freesurfer_func(data_path:str, subj:str, sess:str, isForce:bool, lesion:bool=False):    
    logging.info('Running freesurfer: {0}.'.format(subj + " / " + sess))

    freesurfer_folder = os.path.join(data_path, 'derivatives', '01_freesurfer') 
    os.environ["SUBJECTS_DIR"] = freesurfer_folder

    # Create the folder mri/orig
    if not os.path.exists(os.path.join(freesurfer_folder, subj + "-" + sess, "mri", "orig")):
        os.makedirs(os.path.join(freesurfer_folder, subj + "-" + sess, "mri", "orig"))

    t1_for_fs_pre = os.path.join(freesurfer_folder, subj + "-" + sess, "mri", "orig", subj + "_" + sess + "_acq-mprage_T1w.nii.gz")
    t1 = os.path.join(data_path, subj, sess, "anat", subj + "_" + sess + "_T1w.nii.gz")
    
    if lesion==True: # Fabienne paths not tested
        session_folder = os.path.join(data_path, "derivatives", "03_dwi")
        lesion_transplantation_folder = os.path.join(session_folder, "0_lesion_transplantations_FW", subj, sess, "anat", "lesion_transplantation")
        lesion_transplantation_prefix = os.path.join(lesion_transplantation_folder, subj + '_' + sess)
        T1_transplanted_lesion = (lesion_transplantation_prefix + '_T1w_with_transplanted_lesion.nii.gz')
        if (os.path.isfile(T1_transplanted_lesion)):
            t1 = T1_transplanted_lesion
            
    # Copy the T1 file for fs to free surfer folder in t1 for fs
    if (not os.path.isfile(t1_for_fs_pre)) and not (isForce):
        logging.info('Copying {0} into fs folder.'.format(os.path.basename(t1)))
        shutil.copyfile(t1, t1_for_fs_pre)
   
    #------------------------------------------------------------#
    #### 1 MRI CONVERT CMD #### 
    #------------------------------------------------------------#

    print('#### mri convert ####')
    t1_for_fs = os.path.join(freesurfer_folder,subj + "-" + sess, 'mri', "orig", "001.mgz")
    json_file = os.path.join(freesurfer_folder,subj + "-" + sess, 'mri', "orig", "001.json")
    
    if os.path.isfile(t1_for_fs) and not isForce:
        logging.info('T1 for freesurfer already converted')
    else: 
        mr_convert_cmd = "mrconvert " + t1_for_fs_pre + " " + t1_for_fs 
        logging.info('mr convert command: "{0}".'.format(mr_convert_cmd))
        subprocess.call(mr_convert_cmd, shell=True)
        with open(json_file, 'w') as outfile:
            j = {
                'Origin function': mr_convert_cmd,
                'Description': 'Conversion of T1 for freesurfer',
                'Anat_filename': t1_for_fs,
                'Time' : time.asctime()
                }
            json.dump(j, outfile)

    #------------------------------------------------------------#
    #### 2 RECON ALL CMD #### 
    #------------------------------------------------------------#

    print('#### recon all ####')
    fs_output = os.path.join(freesurfer_folder, "brain.mgz")
    if os.path.isfile(fs_output) and not isForce:
         logging.info('Freesurfer already run')
    else: 
         reconall_cmd = "recon-all -all -subjid " + subj + "-" + sess + " -openmp 12 -brainstem-structures" ## PB HERE
         logging.info('recon all command: "{0}".'.format(reconall_cmd))
         subprocess.call(reconall_cmd, shell=True)

    #------------------------------------------------------------#
    #### 3 SEGMENT BS CMD #### 
    #------------------------------------------------------------#
    print('#### segment BS ####')
    fs_bs_output = os.path.join(freesurfer_folder, "brainstemSsLabels.v13.mgz")
    if os.path.isfile(fs_bs_output) and not isForce:
         logging.info('Freesurfer already run')
    else: 
         #segment_bs_cmd = "segmentBS.sh " + subj + "-" + sess + " " + freesurfer_folder #FOR FREESURFER v 7
         segment_bs_cmd = "recon-all -s " + subj + "-" + sess + " -brainstem-structures" #FOR FREESURFER v 6
         logging.info('segment BS command: "{0}".'.format(segment_bs_cmd))
         subprocess.call(segment_bs_cmd, shell=True)

if __name__ == "__main__":  
    print('Starting free surfer module...')
    parser = buildArgsParser()
    args = parser.parse_args()

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
    fail_list_filename = f"fail_list_03_freesurfer_{formatted_datetime}.txt"
    for subj, sess in itertools.product(subjects, sessions):
        
        if (not (subj in ("sub-TIMESwp11s036"))) & \
            (not (subj in ("sub-TIMESwp11s063"))):
            try:
                freesurfer_func(data_path, subj, sess, args.isForce, args.lesion)
            except Exception as e:
                with open(fail_list_filename, "+a") as f:
                    f.write(f"{subj} {sess} \n")
                    f.write(f"{str(e)} \n")

