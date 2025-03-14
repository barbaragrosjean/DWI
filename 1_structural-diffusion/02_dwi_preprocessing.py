#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division

import argparse
import logging
import os
import json
import shutil
import subprocess
import time
import os.path
import nibabel as nib
import numpy as np
import itertools
from datetime import datetime
from tools.registration_ants import *
import resource

def buildArgsParser():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter,
        epilog="")
    p._optionals.title = "Generic options"

    p.add_argument('--subj', nargs='+', dest='subj', help="Subject index.")
    p.add_argument('--sess', nargs='+', dest='sess', help="Session folder name.")

    p.add_argument('--data_path', default='/mnt/Hummel-Data/TI/mri/51T', dest='data_path',
        help="Subjects folder path. ['%(default)s']")

    current_working_directory = os.getcwd()

    p.add_argument('--acqparams_file', default= current_working_directory + '/eddy/acqparams.txt', 
        help="Aquisition paramters file. ['%(default)s']") 

    p.add_argument('--index_file', default= current_working_directory + '/eddy/eddy_index.txt', 
        help="Aquisition paramters file. ['%(default)s']")

    p.add_argument('-f', action='store_true', dest='isForce', 
    help='If set, overwrites output file.')

    log_g = p.add_argument_group('Logging options')
    log_g.add_argument(
        '-v', action='store_true', dest='isVerbose',
        help='If set, produces verbose output.')
    return p

def pre_proc(data_path:str, subj:str, sess:str, isForce:bool):
    ''' Function doing the preprocessing of the image data, according the folowing steps: 
    degibbs, Extracting b0s, BET, TOPUP, Eddy, and DEBIAS.
    
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

    session_folder = os.path.join(data_path, "derivatives")
    
    if os.path.exists(session_folder):

        logging.info('Processing: "{0}".'.format(session_folder))

        # Input and output data folder
        dwi_raw_folder = os.path.join(data_path, subj, sess, "dwi")
        dwi_target_folder = os.path.join(session_folder, "01_dwi", subj, sess, "dwi", "preproc")
        dwi_out = os.path.join(dwi_target_folder,subj + "_" + sess + "_dwi")

        if os.path.isfile(dwi_out + ".nii.gz") and not isForce:
            logging.info('Preprocessing already done.')
            return    

        # Path to the dwi file
        dwi_ap_filename = os.path.join(dwi_raw_folder, subj + '_' + sess + "_dwi_AP_1.nii.gz")
        dwi_pa_filename = os.path.join(dwi_raw_folder, subj + '_' + sess + "_dwi_PA_1.nii.gz")

        #------------------------------------------------------------#
        #### 1 DEGGIBS #### 
        #------------------------------------------------------------#
        print('#### Degibbs ####')
        dwi_ap_degibbs_filename = os.path.join(dwi_target_folder,subj + "_" + sess + "_dir-AP_degibbsDwi.nii.gz")
        json_file = os.path.join(dwi_target_folder,subj + "_" + sess + "_dir-AP_degibbsDwi.json")
        if os.path.isfile(dwi_ap_degibbs_filename) and not isForce:
            logging.info('mrdegibbs AP already done.')
        else:
            mrdegibbs_cmd = "mrdegibbs " + dwi_ap_filename + " " + dwi_ap_degibbs_filename
            logging.info('mrdegibbs command: "{0}".'.format(mrdegibbs_cmd))
            subprocess.call(mrdegibbs_cmd, shell=True)
            with open(json_file, 'w') as outfile:
                j = {
                    'Origin function': mrdegibbs_cmd,
                    'Description': 'degibbs ap',
                    'dwi_ap_filename': dwi_ap_degibbs_filename,
                    'Time' : time.asctime()
                    }
                json.dump(j, outfile)

        dwi_pa_degibbs_filename = os.path.join(dwi_target_folder,subj + "_" + sess + "_dir-PA_degibbsDwi.nii.gz")
        json_file = os.path.join(dwi_target_folder,subj + "_" + sess + "_dir-PA_degibbsDwi.json")
            
        if os.path.isfile(dwi_pa_degibbs_filename) and not isForce:
            logging.info('mrdegibbs PA already done.')
        else:
            mrdegibbs_cmd = "mrdegibbs " + dwi_pa_filename + " " + dwi_pa_degibbs_filename
            logging.info('mrdegibbs command: "{0}".'.format(mrdegibbs_cmd))
            subprocess.call(mrdegibbs_cmd, shell=True)
            with open(json_file, 'w') as outfile:
                j = {
                    'Origin function': mrdegibbs_cmd,
                    'Description': 'degibbs pa',
                    'dwi_pa_filename': dwi_pa_degibbs_filename,
                    'Time' : time.asctime()
                    }
                json.dump(j, outfile)        

        #------------------------------------------------------------#
        #### 2 EXTRACTING B0s #### 
        #------------------------------------------------------------#
    
        print('#### extracting b0s ####')
        dwi_ap = nib.load(dwi_ap_degibbs_filename)
        dwi_pa = nib.load(dwi_pa_degibbs_filename)

        b0s_filename = os.path.join(dwi_target_folder, subj + "_" + sess + "_dir-APPA_b0s.nii.gz")
        b0s = np.stack([dwi_ap.get_fdata()[:,:,:,0], dwi_pa.get_fdata()[:,:,:,0]],axis=3)
        json_file = os.path.join(dwi_target_folder,subj + "_" + sess + "_dir-APPA_b0s.json")
        if os.path.isfile(b0s_filename) and not isForce:
            logging.info('B0s extracted already done.')
        else:
            nib.Nifti1Image(b0s, dwi_pa.affine, dwi_pa.header).to_filename(b0s_filename)
            orig_func_label = "nib.Nifti1Image(b0s, dwi_pa.affine, dwi_pa.header).to_filename(" + b0s_filename + ")"
            with open(json_file, 'w') as outfile:
                j = {
                    'Origin function': orig_func_label,
                    'Description': 'extracting b0',
                    'b0_filename': b0s_filename,
                    'Time' : time.asctime()
                    }
                json.dump(j, outfile)
            

        b0s_mean_filename = os.path.join(dwi_target_folder, subj + "_" + sess + "_dir-APPA_meanB0.nii.gz")
        b0_mean = np.mean(b0s,axis=3)
        json_file = os.path.join(dwi_target_folder, subj + "_" + sess + "_dir-APPA_mean0.json")
        if os.path.isfile(b0s_mean_filename) and not isForce:
            logging.info('B0 mean extracted already done.')
        else:
            nib.Nifti1Image(b0_mean,dwi_pa.affine, dwi_pa.header).to_filename(b0s_mean_filename)
            orig_func_label = "nib.Nifti1Image(b0_mean,dwi_pa.affine, dwi_pa.header).to_filename(" + b0s_mean_filename + ")"
            with open(json_file, 'w') as outfile:
                j = {
                    'Origin function': orig_func_label,
                    'Description': 'mean b0',
                    'meanb0_filename': b0s_mean_filename,
                    'Time' : time.asctime()
                    }
                json.dump(j, outfile)

        #------------------------------------------------------------#
        #### 3 BET #### 
        #------------------------------------------------------------#
        print('#### Bet ####')
        b0s_mean_brain_filename = os.path.join(dwi_target_folder,subj + "_" + sess + "_dir-APPA_meanB0brain.nii.gz")
        json_file = os.path.join(dwi_target_folder,subj + "_" + sess + "_dir-APPA_meanB0brain.json")
        if os.path.isfile(b0s_mean_brain_filename) and not isForce:
            logging.info('Bet already done.')
        else:
            bet_cmd = "bet " + b0s_mean_filename + " " + b0s_mean_brain_filename + " -f 0.4 -g 0"
            logging.info('Bet command: "{0}".'.format(bet_cmd))
            subprocess.call(bet_cmd, shell=True)
            with open(json_file, 'w') as outfile:
                j = {
                    'Origin function': bet_cmd,
                    'Description': 'bet on mean b0',
                    'bet meanb0_filename': b0s_mean_brain_filename,
                    'Time' : time.asctime()
                    }
                json.dump(j, outfile)

        #------------------------------------------------------------#
        #### 4 TOPUP #### 
        #------------------------------------------------------------# 
        # /!\ time consuming

        print('#### Topup ####')
        topup_out =  os.path.join(dwi_target_folder,subj + "_" + sess + "_topup")
        json_file = os.path.join(dwi_target_folder,subj + "_" + sess + "_topup_fieldcoef.json")
        if os.path.isfile(topup_out + "_fieldcoef.nii.gz") and not isForce:
                logging.info('Topup already done.')
        else:
            topup_cmd = "topup --imain=" + b0s_filename + \
            " --datain=" + args.acqparams_file + \
            " --out=" + topup_out +\
            " --config=b02b0.cnf --subsamp=1" 
            logging.info('Topup command: "{0}".'.format(topup_cmd))
            subprocess.call(topup_cmd, shell=True)
            with open(json_file, 'w') as outfile:
                j = {
                    'Origin function': topup_cmd,
                    'Description': 'topup on b0',
                    'topup field coefficient file': (topup_out + "_fieldcoef.nii.gz"),
                    'Time' : time.asctime()
                    }
                json.dump(j, outfile)        

        #------------------------------------------------------------#
        #### 5 EDDY #### 
        #------------------------------------------------------------#
        
        # /!\ time consuming
        print('#### Eddy ####')
        eddy_out=os.path.join(dwi_target_folder,subj + "_" + sess + "_eddy")
        json_file = os.path.join(dwi_target_folder,subj + "_" + sess + "_eddy.json")
        if os.path.isfile(eddy_out + ".nii.gz") and not isForce:
            logging.info('Eddy already done.')
        else:   
            eddy_cmd ="eddy_openmp --imain=" + dwi_ap_degibbs_filename + " --mask=" + b0s_mean_brain_filename + " --index=" + args.index_file + " --mb=2 --acqp=" + args.acqparams_file +  " --bvals=" + os.path.join(dwi_raw_folder,subj + "_" + sess + "_dwi_AP_1.bval") + " --topup=" + topup_out + " --bvecs=" + os.path.join(dwi_raw_folder,subj + "_" + sess + "_dwi_AP_1.bvec") + " --out=" + eddy_out  + " --data_is_shelled" 
            logging.info('Eddy command: "{0}".'.format(eddy_cmd))
            subprocess.call(eddy_cmd, shell=True)
            # Copy(source, destination) data and right -> corrected bvec and bval : everything we need will be in derivatives
            shutil.copy(eddy_out + ".eddy_rotated_bvecs", dwi_out + ".bvec") 
            shutil.copy(os.path.join(dwi_raw_folder, subj + "_" + sess + "_dwi_AP_1.bval"), dwi_out + ".bval")
            with open(json_file, 'w') as outfile:
                j = {
                    'Origin function': eddy_cmd,
                    'Description': 'eddy on b0',
                    'eddy file': (eddy_out + ".nii.gz"),
                    'Time' : time.asctime()
                    }
                json.dump(j, outfile)        

        #------------------------------------------------------------#
        #### 6 DEBIAS #### 
        #------------------------------------------------------------#

        print('#### Debias ####')
        biasField_out=os.path.join(dwi_target_folder,subj + "_" + sess + "_meanB0brain")
        json_file = os.path.join(dwi_target_folder,subj + "_" + sess + "_meanB0brain_bias.json")
        if os.path.isfile(biasField_out + "_bias.nii.gz") and not isForce:
            logging.info('Debias already done.')
        else:
            fast_cmd="fast -t 2 -n 3 -H 0.1 -I 4 -l 20.0 -b -o " + biasField_out + " " + " " + b0s_mean_brain_filename
            logging.info('Fast debias command: "{0}".'.format(fast_cmd))
            subprocess.call(fast_cmd, shell=True)

            fslmaths_cmd = "fslmaths " + eddy_out + " -div " + biasField_out + "_bias.nii.gz " + dwi_out
            logging.info('Apply debias command: "{0}".'.format(fslmaths_cmd))
            subprocess.call(fslmaths_cmd, shell=True)
            with open(json_file, 'w') as outfile:
                j = {
                    'Origin function': fslmaths_cmd, 
                    'Description': 'debias field on b0',
                    'eddy file': (biasField_out + "_bias.nii.gz"),
                    'Time' : time.asctime()
                    }
                json.dump(j, outfile)       

        #------------------------------------------------------------#
        #### CLEANNING FOLDER #### 
        #------------------------------------------------------------#
         
        print('#### MOVE TO TRASH ####')
        files_to_move = ["_eddy.eddy_command_txt" , "_eddy.eddy_movement_rms", "_eddy.eddy_outlier_map",  "_eddy.eddy_outlier_n_sqr_stdev_map", 
                        "_eddy.eddy_outlier_n_stdev_map", "_eddy.eddy_outlier_report",  "_eddy.eddy_parameters", "_eddy.eddy_parameters", 
                        "_eddy.eddy_parameters", "_eddy.eddy_post_eddy_shell_alignment_parameters", "_eddy.eddy_post_eddy_shell_PE_translation_parameters", 
                        "_eddy.eddy_restricted_movement_rms", "_eddy.eddy_rotated_bvecs", "_eddy.eddy_values_of_all_input_parameters", "_topup_movpar.txt"]
        
        folder_trash = os.path.join(dwi_target_folder, 'trash')

        # If doesn't existe creat a directory trash
        if not os.path.isdir(folder_trash):
            cmd_trash_dir = 'mkdir ' + folder_trash
            subprocess.call(cmd_trash_dir, shell=True)

        # Put the file on it if exist
        for file in files_to_move:
            file_to_move = os.path.join(dwi_target_folder, subj + "_" + sess + file)
            if os.path.exists(file_to_move):
                cmd_move = "mv " + file_to_move + " " + folder_trash
                subprocess.call(cmd_move, shell=True)
        
        if os.path.exists(b0s_filename[:-6]+"topup_log"):
            cmd_move = "mv " + b0s_filename[:-6]+"topup_log" + " " + folder_trash
            subprocess.call(cmd_move, shell=True)

        
        #------------------------------------------------------------#
        #### 7 MEAN B0 #### 
        #------------------------------------------------------------#

        print('#### MEAN B0 ####')
        meanB0_filename = os.path.join(dwi_target_folder, subj + "_" + sess + "_dwi_mean-b0.nii.gz")
        meanB0bet_filename = os.path.join(dwi_target_folder, subj + "_" + sess + "_dwi_mean-b0_bet.nii.gz")
        json_file = os.path.join(dwi_target_folder, subj + "_" + sess + "_dwi_mean-b0.json")
        if os.path.isfile(meanB0_filename) and not isForce:
            logging.info('mean b0 already extracted: "{0}".'.format(meanB0_filename))
        else:
            dwi = nib.load(dwi_out + ".nii.gz")
            bval = np.loadtxt(dwi_out + ".bval")
            b0s = dwi.get_fdata()[:,:,:,bval==0]
            meanB0 = np.mean(b0s, axis=3)
            nib.Nifti1Image(meanB0,dwi.affine,dwi.header).to_filename(meanB0_filename)

            orig_func_label = "nib.Nifti1Image(meanB0,dwi.affine,dwi.header).to_filename(" + meanB0_filename + ") with dwi --> dwi = nib.load(" + dwi_out + ".nii.gz)"
            with open(json_file, 'w') as outfile:
                j = {
                    'Origin function': orig_func_label,
                    'Description': 'create mean b0',
                    'b0_filename': meanB0_filename,
                    'Time' : time.asctime()
                    }
                json.dump(j, outfile)

        #------------------------------------------------------------#
        #### 8 MEAN B0 BET #### 
        #------------------------------------------------------------#

        print('#### MEAN B0 BET ####')
        # Code taken from script https://gitlab.epfl.ch/ebeanato/mcgrase_times/-/blob/main/functions/reg_mcGRASE_proc.py line 194
        json_file = os.path.join(dwi_target_folder, subj + "_" + sess + "_dwi_mean-b0_bet.json")
        if os.path.isfile(meanB0bet_filename):
            logging.info('mean b0 bet already done.: "{0}".'.format(meanB0bet_filename))
        else:
            bet_cmd = "bet " + meanB0_filename + " " + meanB0bet_filename + " -f 0.4 -g 0 -m"
            logging.info('Bet command: "{0}".'.format(bet_cmd))
            print(bet_cmd)
            subprocess.call(bet_cmd, shell=True)
            with open(json_file, 'w') as outfile:
                j = {
                    'Origin function': bet_cmd,
                    'Description': 'bet on mean b0',
                    'b0bet_filename': meanB0bet_filename,
                    'Time' : time.asctime()
                    }
                json.dump(j, outfile)

            
if __name__ == "__main__":
    print('Starting the preprocessing module...')

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
    fail_list_filename = f"fail_list_02_dwi_preprocessing_{formatted_datetime}.txt"
    for subj, sess in itertools.product(subjects, sessions):        
            try:
                pre_proc(data_path, subj, sess, isForce)
            except Exception as e:
                with open(fail_list_filename, "+a") as f:
                    f.write(f"{subj} {sess} \n")
                    f.write(f"{str(e)} \n")
