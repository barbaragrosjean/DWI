#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division

import argparse
import logging
import os
import json
import subprocess
import shutil
import itertools
from joblib import Parallel, delayed
import glob

import nibabel as nib
import numpy as np

def buildArgsParser():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter,
        epilog="")
    p._optionals.title = "Generic options"

    p.add_argument('--subjects', nargs='+', help="Subject folder name.")
    p.add_argument('--sessions', nargs='+', help="Session folder name.")

    p.add_argument('--data_path', default='data_dir', dest='data_path',
        help="Subjects folder path. ['%(default)s']")
    p.add_argument('-f', action='store_true', dest='isForce',
    help='If set, overwrites output file.')

    log_g = p.add_argument_group('Logging options')
    log_g.add_argument(
        '-v', action='store_true', dest='isVerbose',
        help='If set, produces verbose output.')
    return p

def reorganize_data_fnct(data_path:str, subject:str, session:str, isVerbose:bool, isForce:bool):
    """Reorganization script

    Run tractography

    Parameters
    ----------
    data_path :
        Path of the general BIDS structure
    subject :
        Current subject
    session :
        Current session
    output_path :
        Path where to save output files
    isVerbose :
        Boolean indicating if verbose mode should be on
    isForce :
        Boolean indicating if files have to be overwritten

    Examples: python processing_scripts/2_reorganize_data.py --subject P206 --session pre1
    --data_path /mnt/Hummel-data/data_Elena/subjects/patients/ReOrch4Vision -v"""

    source_folder = os.path.join(data_path, "sourcedata", "BIDS", subject, session)
    target_folder = os.path.join(data_path, subject, session)
    if not os.path.exists(target_folder):
        os.makedirs(target_folder)

    derivatives_folder = os.path.join(data_path, "derivatives", subject, session)
    if not os.path.exists(derivatives_folder):
        os.makedirs(derivatives_folder)

    if isVerbose:
        logging.basicConfig(filename=os.path.join(target_folder,'2_reorganize_data.log'), format='%(asctime)s %(message)s', level=logging.DEBUG)

    logging.info('Reorganization started')
    logging.info('Creating new folder structure: "{0}".'.format(target_folder))
    target_folder_anat = os.path.join(target_folder, "anat")
    if not os.path.exists(target_folder_anat):
        os.makedirs(target_folder_anat)
    target_folder_dwi = os.path.join(target_folder, "dwi")
    if not os.path.exists(target_folder_dwi):
        os.makedirs(target_folder_dwi)


    ### ANAT
    anat_folder = os.path.join(source_folder, "anat")
    logging.info('Reorganizing folder: "{0}".'.format(anat_folder))

    if os.path.exists(anat_folder):
        anat_json_files = np.sort([f for f in os.listdir(anat_folder)
                                if os.path.isfile(os.path.join(anat_folder, f))
                                and f.endswith(".json")])
        anat_nii_files = np.sort([f for f in os.listdir(anat_folder)
                                if os.path.isfile(os.path.join(anat_folder, f))
                                and f.endswith(".nii.gz")])

        logging.info('There is {0} .json files in "{1}".'
                    .format(len(anat_json_files), anat_folder))
        logging.info('There is {0} .nii.gz files in "{1}".'
                    .format(len(anat_nii_files), anat_folder))

        if len(anat_json_files) != len(anat_nii_files):
            raise ValueError('Different number of .json and .nii.gz in:"{0}".'
                            .format(anat_folder))

        anat_json = []
        for f in anat_json_files:
            with open(os.path.join(anat_folder, f)) as json_file:
                anat_json.append(json.load(json_file))
        anat_nii = []
        for f in anat_nii_files:
            anat_nii.append(nib.load(os.path.join(anat_folder, f)))

        ## mcgrase
        sessions = dict()
        for j, n, j_f, n_f in zip(anat_json, anat_nii, anat_json_files,
                                anat_nii_files):
            if j['ProtocolName'] == 'mcGRASE_1p6iso_84_AF3x2':
                if 'ShimSetting' in j.keys():
                    session_id = hash((tuple(j['ImageOrientationPatientDICOM']),tuple(j['ShimSetting'])))
                else:
                    session_id = hash(tuple(j['ImageOrientationPatientDICOM']))
                if not session_id in sessions.keys():
                    sessions[session_id] = []
                sessions[session_id].append([j['EchoTime'],
                                            j['AcquisitionTime'],
                                            j,
                                            n,
                                            j_f,
                                            n_f])

        logging.info('There is {0} mcGRASE acquisitions.'
                    .format(len(sessions.keys())))
        for no, k in enumerate(sessions.keys()):
            logging.info('Processing mcGRASE acquisitions: {0}.'
                        .format(sessions[k][0][2]['ImageOrientationPatientDICOM']))
            logging.info('Founded {0} 3d volumes.'
                        .format(len(sessions[k])))
            if len(sessions[k]) == 32:
                json_mcgrass = sessions[k][0][2]
                echo_times = []
                acquisition_times = []
                nii_mcgrass = np.zeros(list(sessions[k][0][3].shape) + [32])
                nii_affine = sessions[k][0][3].affine
                nii_header = sessions[k][0][3].header
                for i in range(32):
                    echo_times.append(sessions[k][i][0])
                    acquisition_times.append(sessions[k][i][1])

                #ordering echo time and images
                idx = np.argsort(echo_times)
                echo_times = [echo_times[e] for e in idx]
                acquisition_times = [acquisition_times[e] for e in idx]
                for i in range(32):
                    nii_mcgrass[:, :, :, i] = sessions[k][idx[i]][3].get_fdata()

                if no==0:
                    mcGRASE_nii_filename = subject + "_" + session + "_mcGRASE.nii.gz"
                    mcGRASE_nii_filename = os.path.join(target_folder_anat, mcGRASE_nii_filename)
                    mcGRASE_json_filename = subject + "_" + session + "_mcGRASE.json"
                    mcGRASE_json_filename = os.path.join(target_folder_anat, mcGRASE_json_filename)
                else:
                    raise ValueError('More than 1 mcGRASE dataset')

                if not os.path.exists(mcGRASE_nii_filename):
                    nib.Nifti1Image(nii_mcgrass, nii_affine, nii_header) \
                        .to_filename(mcGRASE_nii_filename)

                    json_mcgrass['EchoTime'] = echo_times
                    json_mcgrass['AcquisitionTime'] = acquisition_times
                    with open(mcGRASE_json_filename, 'w') as outfile:
                        json.dump(json_mcgrass, outfile)
            else:
                logging.info('Skipping, invalid number of volume.')


        ## t1/t2
        logging.info('Copying anat images.')
        files = ["*_T1w", "*_T2w"]
        for f in files:
            nii_filename = os.path.join(anat_folder, subject + "_" + session + f + ".nii.gz")
            json_filename = os.path.join(anat_folder, subject + "_" + session + f + ".json")

            nii_files = glob.glob(nii_filename)
            json_files = glob.glob(json_filename)

            for idx, nii_file in enumerate(nii_files):
                target_nii_filename = os.path.join(target_folder_anat, os.path.basename(nii_file))
                target_json_filename = os.path.join(target_folder_anat, os.path.basename(json_files[idx]))

                if os.path.exists(nii_file):
                    if not os.path.exists(target_nii_filename):
                        shutil.copyfile(json_files[idx], target_json_filename)
                        subprocess.call("mrconvert -strides 1,2,3 -quiet -force " + nii_file + ".nii.gz " + target_nii_filename + ".nii.gz", shell=True)

    ## DWI
    dwi_folder = os.path.join(source_folder, "dwi")
    logging.info('Reorganizing folder: "{0}".'.format(dwi_folder))

    if os.path.exists(dwi_folder):
        AP_filename = os.path.join(dwi_folder, subject + "_" + session + "_dir-AP_dwi")
        PA_filename = os.path.join(dwi_folder, subject + "_" + session + "_dir-PA_dwi")

        target_AP_filename = os.path.join(target_folder_dwi, subject + "_" + session + "_dir-AP_dwi")
        target_PA_filename = os.path.join(target_folder_dwi, subject + "_" + session + "_dir-PA_dwi")

        if not os.path.exists(target_PA_filename + ".nii.gz"):
            shutil.copyfile(os.path.join(PA_filename + ".json"), target_PA_filename  + ".json")
            shutil.copyfile(os.path.join(PA_filename + ".bval"), target_PA_filename  + ".bval")
            shutil.copyfile(os.path.join(PA_filename + ".bvec"), target_PA_filename  + ".bvec")
            subprocess.call("mrconvert -strides 1,2,3 -quiet -force " + PA_filename + ".nii.gz " + target_PA_filename + ".nii.gz", shell=True)
        if not os.path.exists(target_AP_filename + ".nii.gz"):
            shutil.copyfile(os.path.join(AP_filename + ".json"), target_AP_filename  + ".json")
            shutil.copyfile(os.path.join(AP_filename + ".bval"), target_AP_filename  + ".bval")
            shutil.copyfile(os.path.join(AP_filename + ".bvec"), target_AP_filename  + ".bvec")
            subprocess.call("mrconvert -strides 1,2,3 -quiet -force " + AP_filename + ".nii.gz " + target_AP_filename + ".nii.gz", shell=True)

    logging.info('Reorganization done')
    return


if __name__ == "__main__":
    parser = buildArgsParser()
    args = parser.parse_args()

    data_path = args.data_path
    subjects = ['sub-' + subj for subj in args.subjects]
    sessions = ['ses-' + sess for sess in args.sessions]
    isVerbose = args.isVerbose
    isForce = args.isForce

    Parallel(n_jobs=4)(delayed(reorganize_data_fnct)(data_path, subj, sess, isVerbose, isForce) for subj, sess in itertools.product(subjects, sessions))

