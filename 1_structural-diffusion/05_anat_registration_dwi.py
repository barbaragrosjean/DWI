#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division

import argparse
import logging
import os
import json
import time
import subprocess
import nibabel as nib
import numpy as np
import os.path
import itertools
from datetime import datetime
from tools.registration_ants import *

def buildArgsParser():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter,
        epilog="")
    p._optionals.title = "Generic options"

    p.add_argument('--subj', nargs='+', dest='subj', help="Subject index.")
    p.add_argument('--sess', nargs='+', dest='sess', help="Session folder name.")

    p.add_argument('--data_path', default='/mnt/Hummel-Data/TI/mri/51T', dest='data_path',
        help="Subjects folder path. ['%(default)s']")

    p.add_argument('-f', action='store_false', dest='isForce',
    help='If set, overwrites output file.')

    log_g = p.add_argument_group('Logging options')
    log_g.add_argument('-v', action='store_false', dest='isVerbose', help='If set, produces verbose output.')
    return p

def anat_reg_dwi(data_path:str, subj:str, sess:str, isForce:bool):
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

    session_folder = os.path.join(data_path, "derivatives", "01_dwi", subj, sess)
    freesurfer_folder = os.path.join(data_path, "derivatives", "01_freesurfer", subj + "-" + sess)
    dwi_folder = os.path.join(data_path, subj, sess, "dwi") #raw data
    
    os.environ["SUBJECTS_DIR"] = freesurfer_folder

    if os.path.exists(session_folder):
        logging.info('Processing dataset: "{0}".'.format(session_folder))
        logging.info('Freesurfer dataset: "{0}".'.format(freesurfer_folder))

        anat_folder = os.path.join(session_folder, "anat")
        preproc_folder = os.path.join(session_folder, "dwi", "preproc")
        dwi_base_filename = os.path.join(preproc_folder, subj + "_" + sess + "_dwi")

        if not os.path.exists(os.path.join(anat_folder)):
            os.makedirs(os.path.join(anat_folder))
            logging.info('anat folder created: "{0}".'.format(anat_folder))

        #------------------------------------------------------------#
        #### 1 DWI #### 
        #------------------------------------------------------------#

        print('------ dwi -------')
        print('#### mean b0 ####')
        # Mean b0: in preproc check if exists otherwise do it
        meanB0_filename = os.path.join(preproc_folder, subj + "_" + sess + "_meanB0.nii.gz") 
        json_file = os.path.join(preproc_folder, subj + "_" + sess + "_meanB0.json")
        if os.path.isfile(meanB0_filename) and not isForce:
            logging.info('mean b0 already extracted: "{0}".'.format(meanB0_filename))
        else: 
            dwi = nib.load(dwi_base_filename + ".nii.gz") 
            bval = np.loadtxt(dwi_base_filename + ".bval")
            b0s = dwi.get_fdata()[:,:,:,bval==0]
            meanB0 = np.mean(b0s, axis=3)
            nib.Nifti1Image(meanB0,dwi.affine,dwi.header).to_filename(meanB0_filename)

            orig_func_label = "nib.Nifti1Image(meanB0,dwi.affine,dwi.header).to_filename(" + meanB0_filename + ") with dwi --> dwi = nib.load(" + dwi_base_filename + ".nii.gz)"
            with open(json_file, 'w') as outfile:
                j = { 
                    'Origin function': orig_func_label,
                    'Description': 'create mean b0',
                    'b0_filename': meanB0_filename,
                    'Time' : time.asctime()
                    }
                json.dump(j, outfile)

        # Bet mean b0 in preproc check if exists otherwise do it
        print('#### bet meanB0 ####')
        # Code taken from script https://gitlab.epfl.ch/ebeanato/mcgrase_times/-/blob/main/functions/reg_mcGRASE_proc.py line 194
        meanB0bet_filename = os.path.join(preproc_folder, subj + "_" + sess + "_dwi_mean-b0_bet.nii.gz")
        json_file = os.path.join(preproc_folder, subj + "_" + sess + "_dwi_mean-b0_bet.json")
        if os.path.isfile(meanB0bet_filename) and not isForce:
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


        #------------------------------------------------------------#
        #### 2 T1 #### 
        #------------------------------------------------------------#
        print('------ T1 -------')
        t1_raw = os.path.join(data_path, subj, sess, "anat", subj + "_" + sess + "_T1w.nii.gz")
        t1_base_filename = os.path.join(anat_folder, subj + "_" + sess + "_acq-mprage_T1w") # in derivatives/anat
        
        print('#### T1 bet ####')
        t1_brain_filename = t1_base_filename + "brain.nii.gz"
        json_file = t1_base_filename + "brain.json"
    
        if os.path.isfile(t1_brain_filename) and not isForce:
            logging.info('T1w brain already extracted: "{0}".'.format(t1_brain_filename))
        else:
            bet_cmd = "bet " + t1_raw +" " + t1_brain_filename + " -B -f 0.2 -g -0.2 -o -m -s -v"
            logging.info('t1 Bet command: "{0}".'.format(bet_cmd))
            subprocess.call(bet_cmd, shell=True)
            with open(json_file, 'w') as outfile:
                j = {
                    'Origin function': bet_cmd,
                    'Description': 'bet on T1',
                    'Anat_filename': t1_brain_filename,
                    'Time' : time.asctime()
                    }
                json.dump(j, outfile)            

        # Fast on anat/T1w - segmentation into tissue types
        print('#### fast on anat/T1w ####') 
        if os.path.isfile(t1_brain_filename[:-7] +"PveWM.nii.gz")  and not isForce:
            logging.info('Fast already performed: "{0}".'.format(t1_brain_filename[:-7] +"PveWM.nii.gz"))
        else:
            fast_cmd = "fast -n 3 -t 1 -g -v -o " + t1_brain_filename[:-7] + " " + t1_brain_filename
            logging.info('Fast command: "{0}".'.format(fast_cmd))
            subprocess.call(fast_cmd, shell=True)
            
            print('Cleaning ...')
            # Place file in a folder trash if not used
            files_to_move = ["_seg_0.nii.gz", "_seg_1.nii.gz", "_seg_2.nii.gz", "_seg.nii.gz", "_pveseg.nii.gz", "_mixeltype.nii.gz", "_mask.nii.gz", "_overlay.nii.gz", "_skull.nii.gz"]
            folder_trash = os.path.join(anat_folder, 'trash')
            if not os.path.isdir(folder_trash):
                cmd_trash_dir = 'mkdir ' + folder_trash
                subprocess.call(cmd_trash_dir, shell=True)
            for file in files_to_move:
                file_to_move = os.path.join(anat_folder, t1_brain_filename[:-7] + file) 
                if os.path.exists(file_to_move):
                    cmd_move = "mv " + t1_brain_filename[:-7] + file + " " + folder_trash
                    subprocess.call(cmd_move, shell=True)
            
            # Rename the ones that will be used
            os.rename(t1_brain_filename[:-7] + "_pve_0.nii.gz", t1_brain_filename[:-7] +"PveCSF.nii.gz")
            os.rename(t1_brain_filename[:-7] + "_pve_1.nii.gz", t1_brain_filename[:-7] +"PveGM.nii.gz")
            os.rename(t1_brain_filename[:-7] + "_pve_2.nii.gz", t1_brain_filename[:-7] +"PveWM.nii.gz")

            for label_pve in ["CSF", "GM", "WM"]:
                json_file = t1_base_filename + "Pve" + label_pve + ".json"
                with open(json_file, 'w') as outfile:
                    j = {
                        'Origin function': fast_cmd,
                        'Description': 'fast on T1, segmentation ' + label_pve + ' file',
                        'Anat_filename': t1_brain_filename,
                        'Time' : time.asctime()
                        }
                    json.dump(j, outfile)             

        #------------------------------------------------------------#
        #### 3 REGISTRATION #### 
        #------------------------------------------------------------#
        
        print('------ Registration -------')
        # registration brain from t1 to b0
        print('#### T1 raw-> b0 ####')
        output_file = os.path.join(preproc_folder,subj + "_" + sess + "_acq-mprage_T1w_dwi.nii.gz")
        json_file = os.path.join(preproc_folder,subj + "_" + sess + "_acq-mprage_T1w_dwi.json")
        if os.path.isfile(output_file) and not isForce: 
            logging.info('ANTS already performed: "{0}".'.format(output_file))
        else:
            input_file = t1_raw
            warp_folder = os.path.join(session_folder, "warps")
            warp_name = "T1w2meanB0_ants"
            original_file = t1_brain_filename
            ref_file = meanB0bet_filename
            registerAnts(input_file, output_file, warp_folder, warp_name, original_file, ref_file)
            orig_func_label = "registerAnts(" + input_file + "," +  output_file + "," + warp_folder + "," + warp_name + "," + original_file + "," + ref_file + ")"
            with open(json_file, 'w') as outfile:
                j = {
                    'Origin function': orig_func_label,
                    'Description': 'register T1 to b0',
                    'Anat_filename': output_file,
                    'Time' : time.asctime()
                    }
                json.dump(j, outfile)


        # register tissue maps from t1 to b0
        print('#### T1 CSF -> b0 ####')
        if os.path.isfile(os.path.join(preproc_folder, subj + "_" + sess + "_acq-mprage_T1wPveCSF_dwi.nii.gz")) and not isForce:
            logging.info('WARP already aplied: "{0}".'.format(os.path.join(preproc_folder,subj + "_" + sess + "_acq-mprage_T1wPveCSF.nii.gz")))
        else:
            warp_folder = os.path.join(session_folder, "warps")
            warp_name = "T1w2meanB0_ants"
            original_file = t1_brain_filename
            ref_file = meanB0bet_filename
            for label_pve in ["CSF", "GM", "WM"]:
                input_file = t1_brain_filename[:-7] + "Pve" + label_pve + ".nii.gz"
                output_file = os.path.join(preproc_folder, subj + "_" + sess + "_acq-mprage_T1wPve" + label_pve + "_dwi.nii.gz")
                registerAnts(input_file, output_file, warp_folder, warp_name, original_file, ref_file) 
                orig_func_label = "registerAnts(" + input_file + "," +  output_file + "," + warp_folder + "," + warp_name + "," + original_file + "," + ref_file + ")"

                json_file = os.path.join(preproc_folder, subj + "_" + sess + "_acq-mprage_T1wPve" + label_pve + ".json")
                with open(json_file, 'w') as outfile:
                    j = {
                        'Origin function': orig_func_label,
                        'Description': 'registering tissue types from T1 to b0 ' + label_pve + ' file',
                        'Anat_filename': output_file,
                        'Time' : time.asctime()
                        }
                    json.dump(j, outfile)  
           

        # 5 tissue types file
        print('#### Create 5 tissue type file ####')
        tt5_file = os.path.join(preproc_folder,subj + "_" + sess + "_acq-mprage_T1wPve5tt_dwi.nii.gz")
        json_file = os.path.join(preproc_folder,subj + "_" + sess + "_acq-mprage_T1wPve5tt_dwi.json")
        if os.path.isfile(tt5_file) and not isForce:
            logging.info('5TT file already generated: "{0}".'.format(tt5_file))
        else:
            csf = nib.load(os.path.join(preproc_folder,subj + "_" + sess + "_acq-mprage_T1wPveCSF_dwi.nii.gz"))
            gm = nib.load(os.path.join(preproc_folder,subj + "_" + sess + "_acq-mprage_T1wPveGM_dwi.nii.gz"))
            wm = nib.load(os.path.join(preproc_folder,subj + "_" + sess + "_acq-mprage_T1wPveWM_dwi.nii.gz"))
            empty = np.zeros(csf.shape)
            tt5 = np.stack([gm.get_fdata(),empty,wm.get_fdata(),csf.get_fdata(),empty],axis=3)
            nib.Nifti1Image(tt5, gm.affine, gm.header).to_filename(tt5_file)

            orig_func_label = "nib.Nifti1Image(tt5, gm.affine, gm.header).to_filename(" + tt5_file + ")"
            with open(json_file, 'w') as outfile:
                j = {
                    'Origin function': orig_func_label,
                    'Description': 'create 5 tissue types file',
                    'tt5_filename': tt5_file,
                    'Time' : time.asctime()
                    }
                json.dump(j, outfile)
    

    ########################################################################################
    ######################## Commented whilst Freesurfer is running ########################
    ######################## To be run once Freesurfer is completed ########################
    ##################### Warps need to be updated to T1w2meanB0_ants ######################
    ########################################################################################

        #------------------------------------------------------------#
        #### 5 FROM FS TO T1 TO DWI #### 
        #------------------------------------------------------------#

        # aparc+aseg+bss (brain stem segmentation) to T1 space
        #--> to pass from freesurfer space to DWI space, it's best to do it in this order: freesurfer->T1->DWI
        print('#### fs -> T1 -> dwi ####')
        aparcaseg = os.path.join(freesurfer_folder, "mri", "aparc.a2009s+aseg")
        aparcasegbsssub = os.path.join(os.path.join(anat_folder, subj + "_" + sess + "_acq-mprage_T1w") + "AparcA2009sAsegBSS.nii.gz")
        aparcasegsub = os.path.join(os.path.join(anat_folder, subj + "_" + sess + "_acq-mprage_T1w") + "AparcA2009sAseg.nii.gz")
        bsssub = os.path.join(os.path.join(anat_folder, subj + "_" + sess + "_acq-mprage_T1w") + "BrainstemSsLabels.nii.gz") #pb here 
        json_file = os.path.join(os.path.join(anat_folder, subj + "_" + sess + "_acq-mprage_T1w") + "BrainstemSsLabels.json")
        if os.path.isfile(aparcasegbsssub) and not isForce:
            logging.info('aparc+aseg already in t1 space: "{0}".'.format(aparcasegbsssub))
        else:
            vol2vol_cmd = "mri_vol2vol --targ " + t1_brain_filename + " --mov " + aparcaseg + ".mgz --o " + aparcasegsub + " --regheader --interp nearest"
            logging.info('VOL2VOL command: "{0}".'.format(vol2vol_cmd))
            subprocess.call(vol2vol_cmd, shell=True)

            vol2vol_cmd = "mri_vol2vol --mov " + os.path.join(freesurfer_folder, "mri", "brainstemSsLabels.v10.FSvoxelSpace.mgz") + " --targ " + os.path.join(freesurfer_folder, "mri", "rawavg.mgz") + " --regheader --o "+ bsssub + " --no-save-reg --interp nearest"
            subprocess.call(vol2vol_cmd, shell=True) 

            aparcasegsub_img = nib.load(aparcasegsub)
            bsssub_data = nib.load(bsssub).get_fdata()
            aparcasegbsss_data = aparcasegsub_img.get_fdata()
            aparcasegbsss_data[aparcasegbsss_data==16] = 170 # brainstem
            aparcasegbsss_data[bsssub_data==171] = 171 # DCG
            aparcasegbsss_data[bsssub_data==172] = 172 # Vermis
            aparcasegbsss_data[bsssub_data==173] = 173 # Midbrain
            aparcasegbsss_data[bsssub_data==174] = 174 # Pons
            aparcasegbsss_data[bsssub_data==175] = 175 # Medulla
            aparcasegbsss_data[bsssub_data==177] = 177 # Vermis-White-Matter
            aparcasegbsss_data[bsssub_data==178] = 178 # SPC
            aparcasegbsss_data[bsssub_data==179] = 179 # Floculus
            nib.Nifti1Image(aparcasegbsss_data, aparcasegsub_img.affine, aparcasegsub_img.header).to_filename(aparcasegbsssub)
            with open(json_file, 'w') as outfile:
                j = {
                    'Origin function': vol2vol_cmd,
                    'Description': 'brain stem segmentation',
                    'Anat_filename': bsssub,
                    'Time' : time.asctime()
                    }
                json.dump(j, outfile)

        # aparc+aseg to dwi space
        dwi_aparc_filename = os.path.join(preproc_folder, subj + "_" + sess + "_acq-mprage_T1wAparcA2009sAseg_dwi.nii.gz")
        json_file = os.path.join(preproc_folder, subj + "_" + sess + "_acq-mprage_T1wAparcA2009sAseg_dwi.json")
        if os.path.isfile(dwi_aparc_filename) and not isForce:
            logging.info('aparc+aseg already in dwi space: "{0}".'.format(dwi_aparc_filename))
        else:
            warp_folder = os.path.join(session_folder, "warps")
            warp_name = "T1w2meanB0_ants"
            original_file = t1_brain_filename
            ref_file = meanB0bet_filename
            inv = False
            input_file = aparcasegsub
            output_file = dwi_aparc_filename
            interp_meth = "MultiLabel"
            registerAnts(input_file, output_file, warp_folder, warp_name, original_file, ref_file, inv, interp_meth)
            orig_func_label = "registerAnts(" + input_file + "," +  output_file + "," + \
                warp_folder + "," + warp_name + "," + original_file + "," + ref_file + "," + interp_meth + ")"            
            with open(json_file, 'w') as outfile:
                j = {
                    'Origin function': orig_func_label,
                    'Description': 'register aparc+aseg to b0',
                    'Anat_filename': dwi_aparc_filename,
                    'Time' : time.asctime()
                    }
                json.dump(j, outfile)            


        # aparc+aseg+bss to dwi space
        dwi_aparcbss_filename = os.path.join(preproc_folder,subj + "_" + sess + "_acq-mprage_T1wAparcA2009sAsegBSS_dwi.nii.gz")
        json_file = os.path.join(preproc_folder,subj + "_" + sess + "_acq-mprage_T1wAparcA2009sAsegBSS_dwi.json")        
        if os.path.isfile(dwi_aparcbss_filename) and not isForce:
            logging.info('aparc+aseg+bss already in dwi space: "{0}".'.format(dwi_aparcbss_filename))
        else:
            warp_folder = os.path.join(session_folder, "warps")
            warp_name = "T1w2meanB0_ants"
            original_file = t1_brain_filename
            ref_file = meanB0bet_filename
            inv = False
            input_file = aparcasegbsssub
            output_file = dwi_aparcbss_filename
            interp_meth = "MultiLabel"
            registerAnts(input_file, output_file, warp_folder, warp_name, original_file, ref_file, inv, interp_meth)
            orig_func_label = "registerAnts(" + input_file + "," +  output_file + "," + \
                warp_folder + "," + warp_name + "," + original_file + "," + ref_file + "," + interp_meth + ")"            
            with open(json_file, 'w') as outfile:
                j = {
                    'Origin function': orig_func_label,
                    'Description': 'register aparc+aseg+bss to b0',
                    'Anat_filename': dwi_aparcbss_filename,
                    'Time' : time.asctime()
                    }
                json.dump(j, outfile) 

        # wmparc to t1 space
        wmparc = os.path.join(anat_folder, subj + "_" + sess + "_acq-mprage_T1wWmparc.nii.gz")
        wmparc_filename_bss = os.path.join(anat_folder, subj + "_" + sess + "_acq-mprage_T1wWmparcBSS.nii.gz")
        json_file = os.path.join(anat_folder, subj + "_" + sess + "_acq-mprage_T1wWmparcBSS.json")
        if os.path.isfile(wmparc) and not isForce:
            logging.info('wmparc in t1 space: "{0}".'.format(aparcasegbsssub))
        else:
            vol2vol_cmd = "mri_vol2vol --targ " + t1_brain_filename + " --mov " + os.path.join(freesurfer_folder, "mri", "wmparc") + ".mgz --o " + wmparc + " --regheader --interp nearest"
            logging.info('VOL2VOL command: "{0}".'.format(vol2vol_cmd))
            subprocess.call(vol2vol_cmd, shell=True)
            bsssub_data = nib.load(bsssub).get_fdata()

            wmparc_img = nib.load(wmparc)
            wmparc_img_data = wmparc_img.get_fdata()
            wmparc_img_data[wmparc_img_data==16] = 170 # brainstem
            wmparc_img_data[bsssub_data==171] = 171 # DCG
            wmparc_img_data[bsssub_data==172] = 172 # Vermis
            wmparc_img_data[bsssub_data==173] = 173 # Midbrain
            wmparc_img_data[bsssub_data==174] = 174 # Pons
            wmparc_img_data[bsssub_data==175] = 175 # Medulla
            wmparc_img_data[bsssub_data==177] = 177 # Vermis-White-Matter
            wmparc_img_data[bsssub_data==178] = 178 # SPC
            wmparc_img_data[bsssub_data==179] = 179 # Floculus
            nib.Nifti1Image(wmparc_img_data, wmparc_img.affine, wmparc_img.header).to_filename(wmparc_filename_bss)
            with open(json_file, 'w') as outfile:
                j = {
                    'Origin function': vol2vol_cmd,
                    'Description': 'wm parcellation in t1 space',
                    'Anat_filename': wmparc_filename_bss,
                    'Time' : time.asctime()
                    }
                json.dump(j, outfile)

        # wmparc & wmparc+bss to dwi space
        dwi_wmparc_filename = os.path.join(preproc_folder,subj + "_" + sess + "_acq-mprage_T1wWmparc_dwi.nii.gz")
        dwi_wmparc_filename_bss = os.path.join(preproc_folder,subj + "_" + sess + "_acq-mprage_T1wWmparcBSS_dwi.nii.gz")
        json_file = os.path.join(preproc_folder,subj + "_" + sess + "_acq-mprage_T1wWmparcBSS_dwi.json")
        if os.path.isfile(dwi_wmparc_filename) and os.path.isfile(dwi_wmparc_filename_bss) and not isForce:
            logging.info('wmparc already in dwi space: "{0}".'.format(dwi_wmparc_filename))
        else:
            warp_folder = os.path.join(session_folder, "warps")
            warp_name = "T1w2meanB0_ants"
            original_file = t1_brain_filename
            ref_file = meanB0bet_filename
            inv = False

            input_file = wmparc
            output_file = dwi_wmparc_filename
            interp_meth = "MultiLabel"
            registerAnts(input_file, output_file, warp_folder, warp_name, original_file, ref_file, inv, interp_meth)
            orig_func_label = "registerAnts(" + input_file + "," +  output_file + "," + \
                warp_folder + "," + warp_name + "," + original_file + "," + ref_file + "," + interp_meth + ")"            
            with open(json_file, 'w') as outfile:
                j = {
                    'Origin function': orig_func_label,
                    'Description': 'register wmparc to b0',
                    'Anat_filename': dwi_wmparc_filename,
                    'Time' : time.asctime()
                    }
                json.dump(j, outfile)  

            input_file = wmparc_filename_bss
            output_file = dwi_wmparc_filename_bss
            interp_meth = "MultiLabel"
            registerAnts(input_file, output_file, warp_folder, warp_name, original_file, ref_file, inv, interp_meth)
            orig_func_label = "registerAnts(" + input_file + "," +  output_file + "," + \
                warp_folder + "," + warp_name + "," + original_file + "," + ref_file + "," + interp_meth + ")"            
            with open(json_file, 'w') as outfile:
                j = {
                    'Origin function': orig_func_label,
                    'Description': 'register wmparc+bss to b0',
                    'Anat_filename': dwi_aparcbss_filename,
                    'Time' : time.asctime()
                    }
                json.dump(j, outfile)         
    
    else:
        print(["subj " + subj + ", sess " + sess + " not existing"])

    return
        

if __name__ == "__main__":
    parser = buildArgsParser()
    args = parser.parse_args()
    isForce = False
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
    fail_list_filename = f"fail_list_04_anat_registration_dwi_{formatted_datetime}.txt"
    for subj, sess in itertools.product(subjects, sessions):
        if (not (subj in ("sub-TIMESwp11s036"))) & \
            (not (subj in ("sub-TIMESwp11s063"))):        
            try:
                anat_reg_dwi(data_path, subj, sess, isForce)
            except Exception as e:
                with open(fail_list_filename, "+a") as f:
                    f.write(f"{subj} {sess} \n")
                    f.write(f"{str(e)} \n")

