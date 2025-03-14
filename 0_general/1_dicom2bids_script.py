#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import logging
import os
import subprocess
import shutil

def buildArgsParser():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter,
        epilog="")
    p._optionals.title = "Generic options"

    p.add_argument('--subjects', nargs='+', dest='subjects', help="Subject number ID (i.e. P206)")
    p.add_argument('--sessions', nargs='+', dest='sessions', help="Session name (i.e. pre1)")

    p.add_argument('--data_path', default='dimom_dir', dest='data_path',
        help="Subjects folder path. ['%(default)s']")
    p.add_argument('--output_path', default='out_dir', dest='output_path',
                   help="Subjects folder path. ['%(default)s']")
    p.add_argument('-f', action='store_true', dest='isForce',
    help='If set, overwrites output file.')

    log_g = p.add_argument_group('Logging options')
    log_g.add_argument(
        '-v', action='store_true', dest='isVerbose',
        help='If set, produces verbose output.')
    return p

def dcm2bids_fnct(data_path, output_data_folder, subj, sess):
    dicom_folder = os.path.join(data_path, subj + '_' + sess + "/")

    config_file = "path_to_config"
    
    if os.path.isdir(dicom_folder):
        dcm2bids_cmd = "dcm2bids -d " + dicom_folder + " -p " + subj + " -s " + sess + " -o " + output_data_folder + " -c " + config_file
        logging.info('dcm2bids command: "{0}".'.format(dcm2bids_cmd))
        subprocess.call(dcm2bids_cmd, shell=True)
    return

def main():

    parser = buildArgsParser()
    args = parser.parse_args()

    data_path = args.data_path
    subjects = [subj for subj in args.subjects]
    sessions = [sess for sess in args.sessions]
    if args.isVerbose:
        logging.basicConfig(level=logging.DEBUG)

    output_path = args.output_path

    for subj in subjects:
        for sess in sessions:
            dcm2bids_fnct(data_path, output_path, subj, sess)

    
if __name__ == "__main__":
    """Dicom to BIDS conversion
    Examples: python processing_scripts/1_dicom2bids_script.py --subjects P206 --sessions pre1
    --data_path /mnt/hummelarch/AVANCER/Data/raw/stroke/ 
    --output_path /mnt/hummelarch/AVANCER/Data/analysis/stroke/group_analysis/mri/sourcedata/BIDS"""
    main()
