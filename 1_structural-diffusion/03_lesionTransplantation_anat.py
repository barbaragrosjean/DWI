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
import itertools
import os.path
import nibabel as nib
import numpy as np
from utils import *
from tools.registration_ants import *
from datetime import datetime

def buildArgsParser():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter,
        epilog="")
    p._optionals.title = "Generic options"

    p.add_argument('--subj', nargs='+', dest='subj', help="Subject index.")
    p.add_argument('--sess', nargs='+', dest='sess', help="Session folder name.")

    p.add_argument('--data_path', default='/media/windel/Elements/TiMeS/WP11_MRI/data/times_wp11',
                   dest='data_path', help="Data folder path. ['%(default)s']")
    p.add_argument('--output_path', default='/media/windel/Elements/TiMeS/WP11_MRI/data/times_wp11/derivatives/03_dwi/0_lesion_transplantations_FW',
                   dest='output_path', help="Output folder path. ['%(default)s']")

    p.add_argument('-f', action='store_true', dest='isForce',
    help='If set, overwrites output file.')

    log_g = p.add_argument_group('Logging options')
    log_g.add_argument(
        '-v', action='store_true', dest='isVerbose',
        help='If set, produces verbose output.')
    return p


def lesionTransplantation_anat(data_path:str, output_path:str, subj:str, sess:str, isForce:bool):

    data_session_folder = os.path.join(data_path, subj, sess)
    data_anat_folder = os.path.join(data_session_folder, "anat")

    output_session_folder = os.path.join(output_path, subj, sess)
    output_anat_folder = os.path.join(output_session_folder, "anat")
    output_lesion_folder = os.path.join(output_session_folder, "lesion")
    output_lesion_transpl_folder = os.path.join(output_anat_folder, "lesion_transplantation")

    if not os.path.exists(output_anat_folder):
        os.makedirs(output_anat_folder)
    if not os.path.exists(output_lesion_folder):
        os.makedirs(output_lesion_folder)
    if not os.path.exists(output_lesion_transpl_folder):
        os.makedirs(output_lesion_transpl_folder)

    logging.info('Brain transplantation in T1w: "{0}".'.format(output_lesion_folder))
    print(" BRAIN TRANSPLANTATION IN T1w ")
    
    anat_file = os.path.join(data_anat_folder, subj + "_" + sess + '_acq-mprage_T1w.nii.gz')
    if not os.path.isfile(anat_file):
        raise FileNotFoundError("T1 not existing: " + anat_file)

    ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ###
    ### %%%%%%%%%%%%%%%%%%% Move lesion to the output folder %%%%%%%%%%%%%%%%%%% ###
    ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ###
    data_acutelesion_file = os.path.join(data_anat_folder, subj + "_" + sess + "_T1w_label-acutelesion_roi.nii.gz")
    data_combinedlesion_file = os.path.join(data_anat_folder, subj + "_" + sess + "_T1w_label-combinedlesion_roi.nii.gz")
    if os.path.isfile(data_combinedlesion_file):
        data_lesion_file = data_combinedlesion_file
        logging.info('Using combined lesion file')
        print("Using combined lesion file")
    else:
        data_lesion_file = data_acutelesion_file
        logging.info('Using acute lesion file')
        print("Using acute lesion file")

    output_lesion_file = os.path.join(output_lesion_folder, subj + "_" + sess + "_T1w_label-lesion_roi.nii.gz")

    if os.path.isfile(data_lesion_file):
        shutil.copyfile(data_lesion_file, output_lesion_file)
    else:
        raise FileNotFoundError("You need to manually draw a lesion mask and save it in: " + data_lesion_file)

    ### %%%%%%%%%%%%%%%%%%%%%%%%%%%% Flip anat image %%%%%%%%%%%%%%%%%%%%%%%%%%%% ###
    anat_file_flipped = os.path.join(output_lesion_transpl_folder, subj + "_" + sess + "_acq-mprage_T1w_flipped.nii.gz")
    json_file = os.path.join(output_lesion_transpl_folder, subj + "_" + sess + "_acq-mprage_T1w_flipped.json")
    if os.path.isfile(anat_file_flipped):
        print("Non affected hemi already flipped")
    else:
        fslswapdim_cmd = "fslswapdim " + anat_file + " -x y z " + anat_file_flipped
        print(fslswapdim_cmd)
        logging.info('fslswapdim command: "{0}".'.format(fslswapdim_cmd))
        subprocess.call(fslswapdim_cmd, shell=True)
        with open(json_file, 'w') as outfile:
            j = {
                'Origin function': fslswapdim_cmd,
                'Description': 'Flipped non affected hemisphere',
                'Anat_filename': anat_file_flipped,
                'Time' : time.asctime()
                }
            json.dump(j, outfile)

    ### %%%%%%%%%%%%%%%% Select voxels of the flipped anatomical %%%%%%%%%%%%%%% ###
    ### %%%%%%%%%%%%%%%%%%%%% image within the lesion mask %%%%%%%%%%%%%%%%%%%%% ###
    lesion_extracted_T1w_file = os.path.join(output_lesion_folder, subj + "_" + sess + "_T1w_label-lesion_extracted_T1w.nii.gz")
    json_file = os.path.join(output_lesion_folder, subj + "_" + sess + "_T1w_label-lesion_extracted_T1w.json")
    if os.path.isfile(lesion_extracted_T1w_file):
        print("Voxels within the mask already extracted")
    else:
        fslmaths_cmd = "fslmaths " + anat_file_flipped + " -mul " + output_lesion_file + " " + lesion_extracted_T1w_file
        print(fslmaths_cmd)
        logging.info('fslmaths command: "{0}".'.format(fslmaths_cmd))
        subprocess.call(fslmaths_cmd, shell=True)
        with open(json_file, 'w') as outfile:
            j = {
                'Origin function': fslmaths_cmd,
                'Description': 'Select voxels within mask',
                'Anat_filename': lesion_extracted_T1w_file,
                'Time' : time.asctime()
                }
            json.dump(j, outfile)

    ### %%%%%%%%%% Substitute original voxels with the extracted ones %%%%%%%%%% ###
    ### %%%%% in order to obtain the reference image for the coregistration %%%% ###
    inversed_lesion_mask_file = os.path.join(output_lesion_folder, subj + "_" + sess + "_T1w_label-lesion_roi_inversed.nii.gz")
    json_file = os.path.join(output_lesion_folder, subj + "_" + sess + "_T1w_label-lesion_roi_inversed.json")
    if os.path.isfile(inversed_lesion_mask_file):
        print(inversed_lesion_mask_file)
        print("Lesion mask already inverted")
    else:
        fslmaths_cmd = "fslmaths " + output_lesion_file + " -sub 1 -abs " + inversed_lesion_mask_file
        print(fslmaths_cmd)
        logging.info('fslmaths command: "{0}".'.format(fslmaths_cmd))
        subprocess.call(fslmaths_cmd, shell=True)
        with open(json_file, 'w') as outfile:
            j = {
                'Origin function': fslmaths_cmd,
                'Description': 'Lesion mask inversion',
                'Anat_filename': inversed_lesion_mask_file,
                'Time' : time.asctime()
                }
            json.dump(j, outfile)

    tmp_anat_without_lesion_mask_file = os.path.join(output_lesion_transpl_folder, subj + "_" + sess + '_T1w_without_lesion_mask.nii.gz')
    json_file = os.path.join(output_lesion_transpl_folder, subj + "_" + sess + '_T1w_without_lesion_mask.json')
    if os.path.isfile(tmp_anat_without_lesion_mask_file):
        print("Voxels within the mask already extracted")
    else:
        fslmaths_cmd = "fslmaths " + anat_file + " -mul " + inversed_lesion_mask_file + " " + tmp_anat_without_lesion_mask_file
        print(fslmaths_cmd)
        logging.info('fslmaths command: "{0}".'.format(fslmaths_cmd))
        subprocess.call(fslmaths_cmd, shell=True)
        with open(json_file, 'w') as outfile:
            j = {
                'Origin function': fslmaths_cmd,
                'Description': 'Lesion removed from the anatomical image',
                'Anat_filename': tmp_anat_without_lesion_mask_file,
                'Time' : time.asctime()
                }
            json.dump(j, outfile)

    anat_no_coregistered_lesion_file = os.path.join(output_lesion_transpl_folder, subj + "_" + sess + '_T1w_with_no_coregistered_lesion.nii.gz')
    json_file = os.path.join(output_lesion_transpl_folder, subj + "_" + sess + '_T1w_with_no_coregistered_lesion.json')
    if os.path.isfile(anat_no_coregistered_lesion_file):
        print("Voxels within the mask already extracted")
    else:
        fslmaths_cmd = "fslmaths " + tmp_anat_without_lesion_mask_file + " -add " + lesion_extracted_T1w_file + " " + anat_no_coregistered_lesion_file
        print(fslmaths_cmd)
        logging.info('fslmaths command: "{0}".'.format(fslmaths_cmd))
        subprocess.call(fslmaths_cmd, shell=True)
        with open(json_file, 'w') as outfile:
            j = {
                'Origin function': fslmaths_cmd,
                'Description': 'Anatomical image with transplanted region',
                'Anat_filename': anat_no_coregistered_lesion_file,
                'Time' : time.asctime()
                }
            json.dump(j, outfile)

    ### %%%%%%%%%%%%%%%%%%%%%%%%%%%% Flip lesion mask %%%%%%%%%%%%%%%%%%%%%%%%%%% ###
    inversed_lesion_mask_flipped = os.path.join(output_lesion_folder, subj + "_" + sess + "_T1w_label-lesion_roi_inversed_flipped.nii.gz")
    json_file = os.path.join(output_lesion_folder, subj + "_" + sess + "_T1w_label-lesion_roi_inversed_flipped.json")
    if os.path.isfile(inversed_lesion_mask_flipped):
        print("Non affected hemi already flipped")
    else:
        fslswapdim_cmd = "fslswapdim " + inversed_lesion_mask_file + " -x y z " + inversed_lesion_mask_flipped
        print(fslswapdim_cmd)
        logging.info('fslswapdim command: "{0}".'.format(fslswapdim_cmd))
        subprocess.call(fslswapdim_cmd, shell=True)
        with open(json_file, 'w') as outfile:
            j = {
                'Origin function': fslswapdim_cmd,
                'Description': 'Flipped inversed lesion mask',
                'Anat_filename': inversed_lesion_mask_flipped,
                'Time' : time.asctime()
                }
            json.dump(j, outfile)

    ### %%%%%%%%%%%% Flip extracted healthy tissue within lesion mask %%%%%%%%%%% ###
    lesion_extracted_T1w_flipped = os.path.join(output_lesion_folder, subj + "_" + sess + "_T1w_label-lesion_extracted_T1w_flipped.nii.gz")
    json_file = os.path.join(output_lesion_folder, subj + "_" + sess + "_T1w_label-lesion_extracted_T1w_flipped.json")
    if os.path.isfile(lesion_extracted_T1w_flipped):
        print("Non affected hemi already flipped")
    else:
        fslswapdim_cmd = "fslswapdim " + lesion_extracted_T1w_file + " -x y z " + lesion_extracted_T1w_flipped
        print(fslswapdim_cmd)
        logging.info('fslswapdim command: "{0}".'.format(fslswapdim_cmd))
        subprocess.call(fslswapdim_cmd, shell=True)
        with open(json_file, 'w') as outfile:
            j = {
                'Origin function': fslswapdim_cmd,
                'Description': 'Flipped extracted healthy tissue within lesion mask',
                'Anat_filename': lesion_extracted_T1w_flipped,
                'Time' : time.asctime()
                }
            json.dump(j, outfile)

    ### %%%%%%%%%%%%%%%% Substitute lesion voxels in the flipped %%%%%%%%%%%%%%% ###
    ### %%%%%%%%%%%%%%% anatomical image with the extracted ones %%%%%%%%%%%%%%% ###
    ### %%%%%%% in order to obtain the input image for the coregistration %%%%%% ###
    tmp_anat_flipped_without_lesion_mask_file = os.path.join(output_lesion_transpl_folder, subj + "_" + sess + '_T1w_flipped_without_lesion_mask.nii.gz')
    json_file = os.path.join(output_lesion_transpl_folder, subj + "_" + sess + '_T1w_flipped_without_lesion_mask.json')
    if os.path.isfile(tmp_anat_flipped_without_lesion_mask_file):
        print("Voxels within the mask already extracted")
    else:
        fslmaths_cmd = "fslmaths " + anat_file_flipped + " -mul " + inversed_lesion_mask_flipped + " " + tmp_anat_flipped_without_lesion_mask_file
        print(fslmaths_cmd)
        logging.info('fslmaths command: "{0}".'.format(fslmaths_cmd))
        subprocess.call(fslmaths_cmd, shell=True)
        with open(json_file, 'w') as outfile:
            j = {
                'Origin function': fslmaths_cmd,
                'Description': 'Lesion removed from the flipped anatomical image',
                'Anat_filename': tmp_anat_flipped_without_lesion_mask_file,
                'Time' : time.asctime()
                }
            json.dump(j, outfile)

    anat_flipped_no_coregistered_lesion_file = os.path.join(output_lesion_transpl_folder, subj + "_" + sess + '_T1w_flipped_with_no_coregistered_lesion.nii.gz')
    json_file = os.path.join(output_lesion_transpl_folder, subj + "_" + sess + '_T1w_flipped_with_no_coregistered_lesion.json')
    if os.path.isfile(anat_flipped_no_coregistered_lesion_file):
        print("Voxels within the mask already extracted")
    else:
        fslmaths_cmd = "fslmaths " + tmp_anat_flipped_without_lesion_mask_file + " -add " + lesion_extracted_T1w_flipped + " " + anat_flipped_no_coregistered_lesion_file
        print(fslmaths_cmd)
        logging.info('fslmaths command: "{0}".'.format(fslmaths_cmd))
        subprocess.call(fslmaths_cmd, shell=True)
        with open(json_file, 'w') as outfile:
            j = {
                'Origin function': fslmaths_cmd,
                'Description': 'Flipped anatomical image with transplanted region',
                'Anat_filename': anat_flipped_no_coregistered_lesion_file,
                'Time' : time.asctime()
                }
            json.dump(j, outfile)

    ### %%%%%%%%%%%%%%%%%%%%% Coregistration of hemispheres %%%%%%%%%%%%%%%%%%%%% ###
    nonAffected2affected_brain_hemi_file = os.path.join(output_lesion_transpl_folder, subj + "_" + sess + "_T1w_flipped2T1w_without_les.nii.gz")
    json_file = os.path.join(output_lesion_transpl_folder, subj + "_" + sess + "_T1w_flipped2T1w_without_les.json")
    if os.path.isfile(nonAffected2affected_brain_hemi_file) and not isForce:
        logging.info('ANTS already performed: "{0}".'.format(nonAffected2affected_brain_hemi_file))
    else:
        input_file = anat_flipped_no_coregistered_lesion_file
        output_file = nonAffected2affected_brain_hemi_file
        warp_folder = os.path.join(output_session_folder, 'warps')
        warp_name = 'T1w_flipped2T1w_without_les'
        original_file = os.path.join(output_lesion_transpl_folder, subj + "_" + sess + '_T1w_flipped_with_no_coregistered_lesion.nii.gz')
        ref_file = anat_no_coregistered_lesion_file
        inv = False
        registerAnts(input_file, output_file, warp_folder, warp_name, original_file, ref_file, inv)
        with open(json_file, 'w') as outfile:
            j = {
                'Origin function': "registerAnts with registration_ants.py" + input_file + "&" + output_file + " ",
                'Description': 'Non affected hemisphere coregistered on the affected one',
                'Anat_filename': nonAffected2affected_brain_hemi_file,
                'Time' : time.asctime()
                }
            json.dump(j, outfile)

    ### %%%%%%%%%%%%%%%%%%%%%%% Select voxels within mask %%%%%%%%%%%%%%%%%%%%%% ###
    lesion_extracted_T1w_coregistered_file = os.path.join(output_lesion_folder, subj + "_" + sess + "_T1w_label-lesion_extracted_T1w_coregistered.nii.gz")
    json_file = os.path.join(output_lesion_folder, subj + "_" + sess + "_T1w_label-lesion_extracted_T1w_coregistered.json")
    if os.path.isfile(lesion_extracted_T1w_coregistered_file):
        print("Voxels within the mask already extracted")
    else:
        fslmaths_cmd = "fslmaths " + nonAffected2affected_brain_hemi_file + " -mul " + output_lesion_file + " " + lesion_extracted_T1w_coregistered_file
        print(fslmaths_cmd)
        logging.info('fslmaths command: "{0}".'.format(fslmaths_cmd))
        subprocess.call(fslmaths_cmd, shell=True)
        with open(json_file, 'w') as outfile:
            j = {
                'Origin function': fslmaths_cmd,
                'Description': 'Select voxels within mask',
                'Anat_filename': lesion_extracted_T1w_file,
                'Time' : time.asctime()
                }
            json.dump(j, outfile)

    ### %%%%%%%%%% Substitute original voxels with the extracted ones %%%%%%%%%% ###
    anat_transplanted_lesion_file = os.path.join(output_lesion_transpl_folder, subj + "_" + sess + '_T1w_with_transplanted_lesion.nii.gz')
    json_file = os.path.join(output_lesion_transpl_folder, subj + "_" + sess + '_T1w_with_transplanted_lesion.json')
    if os.path.isfile(anat_transplanted_lesion_file):
        print("Voxels within the mask already extracted")
    else:
        fslmaths_cmd = "fslmaths " + tmp_anat_without_lesion_mask_file + " -add " + lesion_extracted_T1w_coregistered_file + " " + anat_transplanted_lesion_file
        print(fslmaths_cmd)
        logging.info('fslmaths command: "{0}".'.format(fslmaths_cmd))
        subprocess.call(fslmaths_cmd, shell=True)
        with open(json_file, 'w') as outfile:
            j = {
                'Origin function': fslmaths_cmd,
                'Description': 'Anatomical image with transplanted region',
                'Anat_filename': anat_transplanted_lesion_file,
                'Time' : time.asctime()
                }
            json.dump(j, outfile)

    logging.info('Brain transplantation in T1w done')


if __name__ == "__main__":

    parser = buildArgsParser()
    args = parser.parse_args()

    isForce = args.isForce
    if args.isVerbose:
        logging.basicConfig(level=logging.DEBUG)

    subj_list = [subj for subj in args.subj]
    sess_list = [sess for sess in args.sess]

    data_path = args.data_path
    output_path = args.output_path

    if "all" in subj_list:
        subjects = [s for s in os.listdir(data_path) if os.path.isdir(os.path.join(data_path, s)) if "sub-TIMESwp11s" in s]
    else:
        subjects = ['sub-' + subj for subj in subj_list]

    if "all" in sess_list:
        sessions = ["ses-T1", "ses-T2", "ses-T3", "ses-T4"]
    else:
        sessions = ['ses-' + sess for sess in sess_list]

    #Parallel(n_jobs=4)(delayed(lesionTransplantation_anat)(subj, sess, data_path, output_path, isForce)
                          #for subj, sess in itertools.product(subjects, sessions))
    
    #for subj, sess in itertools.product(subjects, sessions):
        #lesionTransplantation_anat(data_path, output_path, subj, sess, isForce)

    # for subj, sess in itertools.product(subjects, sessions):
    #     try:
    #         lesionTransplantation_anat(data_path, output_path, subj, sess, isForce)
    #     except:
    #         with open("fail_list_02_lesionTransplantation.txt", "+a") as f:
    #             f.write(f"{subj} {sess} \n")

    date = datetime.now()
    formatted_datetime = date.strftime("%Y-%m-%d-%H-%M-%S")
    fail_list_filename = f"fail_list_02_lesionTransplantation_{formatted_datetime}.txt"
    for subj, sess in itertools.product(subjects, sessions):
            try:
                lesionTransplantation_anat(data_path, output_path, subj, sess, isForce)
            except Exception as e:
                with open(fail_list_filename, "+a") as f:
                    f.write(f"{subj} {sess} \n")
                    f.write(f"{str(e)} \n")


        # if (subj == "sub-TIMESwp11s027" and sess == "ses-T2") or \
        #     (subj == "sub-TIMESwp11s036") or \
        #     (subj == "sub-TIMESwp11s039" and sess == "ses-T2") or \
        #     (subj == "sub-TIMESwp11s040") or \
        #     (subj == "sub-TIMESwp11s041") or \
        #     (subj == "sub-TIMESwp11s060") or \
        #     (subj == "sub-TIMESwp11s068") or \
        #     (subj == "sub-TIMESwp11s070" and sess == "ses-T1") or \
        #     (subj == "sub-TIMESwp11s070" and sess == "ses-T3") or \
        #     (subj == "sub-TIMESwp11s071" and sess == "ses-T2") or \
        #     (subj == "sub-TIMESwp11s075") or \
        #     (subj == "sub-TIMESwp11s076" and sess == "ses-T1") or \
        #     (subj == "sub-TIMESwp11s076" and sess == "ses-T2") or \
        #     (subj == "sub-TIMESwp11s076" and sess == "ses-T4") or \
        #     (subj == "sub-TIMESwp11s077" and sess == "ses-T2") or \
        #     (subj == "sub-TIMESwp11s077" and sess == "ses-T3") or \
        #     (subj == "sub-TIMESwp11s077" and sess == "ses-T4") or \
        #     (subj == "sub-TIMESwp11s080" and sess == "ses-T2") or \
        #     (subj == "sub-TIMESwp11s083" and sess == "ses-T1") or \
        #     (subj == "sub-TIMESwp11s083" and sess == "ses-T2"):
                    
            #     if (subj == "sub-TIMESwp11s017" and sess == "ses-T2") or \
            # (subj == "sub-TIMESwp11s027") or \
            # (subj == "sub-TIMESwp11s050" and sess == "ses-T1") or \
            # (subj == "sub-TIMESwp11s050" and sess == "ses-T3") or \
            # (subj == "sub-TIMESwp11s057"):