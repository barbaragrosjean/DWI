#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This file provide a way to register specific voxel from MNI space to dwi space. 
# Must be used after the dwi pipeline once the full tractogram is extracted, as a first step for roi-to-roi analysis.
# Here the voxel used are from the fMRI study + Striatum

from __future__ import division

import argparse
import os
import sys
import time
import itertools
import logging
import json
from datetime import datetime
import nibabel as nib

sys.path.insert(1,'/home/bgrosjea/mnt/Hummel-Data/TI/mri/51T/barbara/uphummel_imaging_template/1_structural-diffusion')
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

    p.add_argument('-f', action='store_true', dest='isForce', 
    help='If set, overwrites output file.')

    log_g = p.add_argument_group('Logging options')
    log_g.add_argument(
        '-v', action='store_true', dest='isVerbose',
        help='If set, produces verbose output.')
    return p


def reg_MNI2B0(data_path:str, subj:str, sess:str, isForce:bool):
    #------------------------------------------------------------#
    # SETUP:  TO CHANGE IF WANT TO RUN THE ANALYSIS FOR OTHER ROIS
    #------------------------------------------------------------#
    study = 'fMRI_study'

    # Path to file which contains the rois
    mni_folder = os.path.join(data_path, "derivatives","01_mni")
    clusterMNI = os.path.join(mni_folder, 'iTBS_vs_HF_control_FDR_001_n_clusters.nii')
    #------------------------------------------------------------#


    session_folder = os.path.join(data_path, "derivatives","01_dwi", subj, sess)
    tract_folder =  os.path.join(data_path, "derivatives", "01_tracts", subj, sess)

    if not os.path.exists(tract_folder) :
        os.makedirs(tract_folder) 
    
    # Usefull folders
    warp_folder = os.path.join(session_folder, "warps") # contains the transformations

    # Anat files
    t1_raw = os.path.join(data_path, subj, sess,'anat',   subj + "_" + sess + "_T1w.nii.gz")
    t1_brain_filename = os.path.join(session_folder, "anat", subj + "_" + sess + "_acq-mprage_T1wbrain.nii.gz")

    if os.path.exists(session_folder): 
        logging.info('Processing dataset: "{0}".'.format(session_folder))

        #------------------------------------------------------------#
        #### 1.1 REGISTRATION CLUSTERS MNI TO T1 SPACE #### 
        #------------------------------------------------------------#
        
        print('#### Register voxel from MNI template to T1 SPACE ####')
        # Inputs
        MNI_file = os.path.join("/usr/local/fsl/data/standard/MNI152_T1_1mm.nii.gz")

        if not os.path.isfile(MNI_file) : 
            logging.info('MNI file must be in: "{0}".'.format(MNI_file))
            print('MNI file must be in: "{0}".'.format(MNI_file))
            return 

        # Ouputs 
        out_folder = os.path.join(tract_folder)
        if not os.path.exists(os.path.join(out_folder, 'roi2roi', study)) :
            os.makedirs(os.path.join(out_folder, 'roi2roi', study)) 

        MNI2tw1 = os.path.join(out_folder, 'roi2roi', study, subj + "_" + sess + "_roi_Clusters_Tw1_ants.nii.gz")
        MNI2tw1_json = os.path.join(out_folder,'roi2roi', study, subj + "_" + sess + "_roi_Clusters_Tw1_ants.json")

        if os.path.isfile(MNI2tw1) and not isForce: 
            logging.info('ANTS already performed: "{0}".'.format(MNI2tw1))
        else:
            input_file = clusterMNI
            warp_name = "MNI2Tw1_ants"
            original_file = MNI_file 
            ref_file = t1_raw # Use not brain extracted
            interp_meth = "MultiLabel"
            registerAnts(input_file, MNI2tw1, warp_folder, warp_name, original_file, ref_file, False, interp_meth)
            orig_func_label = "registerAnts(" + input_file + "," +  MNI2tw1 + "," + warp_folder + "," + warp_name + "," + original_file + "," + ref_file + ")"
            
            with open(MNI2tw1_json, 'w') as outfile:
                j = {
                    'Origin function': orig_func_label,
                    'Description': 'register MNI to tw1',
                    'Anat_filename': MNI2tw1,
                    'Time' : time.asctime()
                    }
                json.dump(j, outfile)

        # CTRL 
        MNItemplate2tw1 = os.path.join(out_folder, 'roi2roi', study, subj + "_" + sess + "_templateMNI_Tw1_ants.nii.gz")
        MNItemplate2tw1_json = os.path.join(out_folder,'roi2roi', study, subj + "_" + sess + "_templateMNI_Tw1_ants.json")

        if os.path.isfile(MNItemplate2tw1) and not isForce: 
            logging.info('ANTS already performed: "{0}".'.format(MNI2tw1))
        else:
            input_file = MNI_file
            warp_name = "MNI2Tw1brain_ants"
            original_file = MNI_file 
            ref_file = t1_raw # Use not brain extracted
            interp_meth = "MultiLabel"
            registerAnts(input_file, MNItemplate2tw1, warp_folder, warp_name, original_file, ref_file, False, interp_meth)
            orig_func_label = "registerAnts(" + input_file + "," +  MNI2tw1 + "," + warp_folder + "," + warp_name + "," + original_file + "," + ref_file + ")"
            
            with open(MNItemplate2tw1_json, 'w') as outfile:
                j = {
                    'Origin function': orig_func_label,
                    'Description': 'register MNI to tw1',
                    'Anat_filename': MNI2tw1,
                    'Time' : time.asctime()
                    }
                json.dump(j, outfile)

        #------------------------------------------------------------#
        #### 1.2 REGISTRATION CLUSTERS MNIT1 TO B0 SPACE  #### 
        #------------------------------------------------------------#
        print('#### Register voxel from MNI template to B0 ####')
        # Voxel found from fMRI study, namely: Loc_NA_Postcentral_L, Rolandic_Oper_R, Loc_NA_Cerebellum, Thal_IL_R, Precentral_L, Supp_Motor_Area_R, Temporal_Sup_L, Insula_L
        # area found in iTBS_vs_HF_control_FDR_001_n_clusters.txt, iTBS_vs_HF_control_FDR_001_n_clusters.nii
        meanB0bet_filename = os.path.join(session_folder, "dwi", "preproc",subj + "_" + sess + "_dwi_mean-b0_bet.nii.gz")

        if not os.path.isfile(meanB0bet_filename) :
            logging.info('meanB0bet file must be process previously')
            return   

        # Ouputs  
        MNI_Tw12B0 = os.path.join(out_folder, 'roi2roi', study, subj + "_" + sess + "_roi_Clusters_dwi_ants.nii.gz")
        MNI_Tw12B0_json = os.path.join(out_folder, 'roi2roi', study, subj + "_" + sess + "_roi_Clusters_dwi_ants.json")

        if not os.path.exists(out_folder) : 
            os.makedirs(out_folder)

        if os.path.isfile(MNI_Tw12B0) and not isForce: 
            logging.info('ANTS already performed: "{0}".'.format(MNI_Tw12B0))
        else:
            input_file = MNI2tw1
            warp_name = "T1w2meanB0_ants"
            original_file = t1_brain_filename
            ref_file = meanB0bet_filename
            inv = False
            interp_meth = "MultiLabel"
            registerAnts(input_file, MNI_Tw12B0, warp_folder, warp_name, original_file, ref_file, inv, interp_meth)
            orig_func_label = "registerAnts(" + input_file + "," +  MNI_Tw12B0 + "," + warp_folder + "," + warp_name + "," + original_file + "," + ref_file + ")"
            
            with open(MNI_Tw12B0_json, 'w') as outfile:
                j = {
                    'Origin function': orig_func_label,
                    'Description': 'register MNItw1 to b0',
                    'Anat_filename': MNI_Tw12B0,
                    'Time' : time.asctime()
                    }
                json.dump(j, outfile)

        # CTRL
        # Ouputs  
        MNItemplatetw12B0 = os.path.join(out_folder, 'roi2roi', study, subj + "_" + sess + "_templateMNI_dwi_ants.nii.gz")
        MNItemplatetw12B0_json = os.path.join(out_folder, 'roi2roi', study, subj + "_" + sess + "_templateMNI_dwi_ants.json")

        if not os.path.exists(out_folder) : 
            os.makedirs(out_folder)

        if os.path.isfile(MNItemplatetw12B0) and not isForce: 
            logging.info('ANTS already performed: "{0}".'.format(MNItemplatetw12B0))
        else:
            input_file = MNItemplate2tw1
            warp_name = "T1w2meanB0_ants"
            original_file = t1_brain_filename
            ref_file = meanB0bet_filename
            interp_meth = "MultiLabel"
            registerAnts(input_file, MNItemplatetw12B0, warp_folder, warp_name, original_file, ref_file, False, interp_meth)
            orig_func_label = "registerAnts(" + input_file + "," +  MNI_Tw12B0 + "," + warp_folder + "," + warp_name + "," + original_file + "," + ref_file + ")"
            
            with open(MNItemplatetw12B0_json, 'w') as outfile:
                j = {
                    'Origin function': orig_func_label,
                    'Description': 'register MNItw1 to b0',
                    'Anat_filename': MNI_Tw12B0,
                    'Time' : time.asctime()
                    }
                json.dump(j, outfile)
        #------------------------------------------------------------#
        #### 2.1 REGISTRATION STRIATUM TO Tw1 and B0 SPACE  #### 
        #------------------------------------------------------------#
        print('#### Register Striatum from ABI atlas to tw1 and to dwi ####')
        # Inputs parcellation of the striatum 
        striatum_files = ['roi_v_d_Ca_L_roi.nii', 'roi_v_d_Ca_R_roi.nii', 'roi_vm_dl_PU_L_roi.nii', 'roi_vm_dl_PU_R_roi.nii']
        
        for file in striatum_files:
            #Input
            MNI_striat = os.path.join(mni_folder, file)

            # Ouputs 
            if not os.path.exists(os.path.join(out_folder, 'striat')) :
                os.makedirs(os.path.join(out_folder, 'striat')) 

            MNIstriat2Tw1 = os.path.join(out_folder, 'striat', subj + "_" + sess + "_" + file[:-8] + "_Tw1_ants.nii.gz")
            MNIstriat2Tw1_json = os.path.join(out_folder,'striat', subj + "_" + sess + "_" + file[:-8] + "_Tw1_ants.json")

            if os.path.isfile(MNIstriat2Tw1) and not isForce: 
                logging.info('ANTS already performed: "{0}".'.format(MNIstriat2Tw1))
            else:
                input_file = MNI_striat
                warp_name = "MNI2Tw1brain_ants"
                original_file = MNI_file 
                ref_file = t1_brain_filename
                interp_meth = "MultiLabel"
                registerAnts(input_file, MNIstriat2Tw1, warp_folder, warp_name, original_file, ref_file, False, interp_meth)
                orig_func_label = "registerAnts(" + input_file + "," +  MNIstriat2Tw1 + "," + warp_folder + "," + warp_name + "," + original_file + "," + ref_file + ")"
                
                with open(MNIstriat2Tw1_json, 'w') as outfile:
                    j = {
                        'Origin function': orig_func_label,
                        'Description': 'register MNI to tw1',
                        'Anat_filename': MNIstriat2Tw1,
                        'Time' : time.asctime()
                        }
                    json.dump(j, outfile)

            #------------------------------------------------------------#
            #### 2.2 REGISTRATION STRIAT Tw1 TO B0 SPACE #### 
            #------------------------------------------------------------#

            # Ouputs 
            MNIstriat2dwi = os.path.join(out_folder,'striat', subj + "_" + sess + "_" + file[:-8] + "_dwi_ants.nii.gz")
            MNIstriat2dwi_json = os.path.join(out_folder,'striat', subj + "_" + sess + "_" + file[:-8] + "_dwi_ants.json")

            if os.path.isfile(MNIstriat2dwi) and not isForce: 
                logging.info('ANTS already performed: "{0}".'.format(MNIstriat2dwi))
            else:
                input_file = MNIstriat2Tw1
                warp_name = "T1w2meanB0_ants"
                original_file = t1_brain_filename 
                ref_file = meanB0bet_filename
                interp_meth = "MultiLabel"
                registerAnts(input_file, MNIstriat2dwi, warp_folder, warp_name, original_file, ref_file, False, interp_meth)
                orig_func_label = "registerAnts(" + input_file + "," +  MNIstriat2dwi + "," + warp_folder + "," + warp_name + "," + original_file + "," + ref_file + ")"
                
                with open(MNIstriat2dwi_json, 'w') as outfile:
                    j = {
                        'Origin function': orig_func_label,
                        'Description': 'register MNItw1 to B0',
                        'Anat_filename': MNIstriat2dwi,
                        'Time' : time.asctime()
                        }
                    json.dump(j, outfile)


    else:
        raise FileNotFoundError("subj " + subj + ", sess " + sess + " not existing")

    return
        

if __name__ == "__main__":
    print('Starting roi-to-roi : registration...')

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
    fail_list_filename = f"fail_list_11_register_rois_MNI2B0{formatted_datetime}.txt"
    for subj, sess in itertools.product(subjects, sessions):        
            try:
                reg_MNI2B0(data_path, subj, sess, isForce)
            except Exception as e:
                with open(fail_list_filename, "+a") as f:
                    f.write(f"{subj} {sess} \n")
                    f.write(f"{str(e)} \n")
