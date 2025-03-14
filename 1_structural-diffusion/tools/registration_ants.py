#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import subprocess
import os
import sys

sys.path
sys.path.append('/opt/ants-2.4.3')


def registerAnts(input_file:str, output_file:str, warp_folder:str, warp_name:str, original_file:str, ref_file:str, inv:bool=False, interp:str="Linear", dim_add:str=""):
    '''Registration from one space to another
        
        Parameters
        ----------
        input_file :
            File to register to a target space for antsApplyTransforms_cmd
        output_file :
            File registered to the target space
        warp_folder :
            Folder path where to save the warp
        warp_name :
            Name of the warp
        original_file :
            File in the same space of the input file for antsRegistrationSyN_cmd
        ref_file :
            File in the target space
        inv :
            Specify if you want to use the inversed warp
        interp :
            Type of interpolation method to use - Please refer to "antsApplyTransforms --help" 
            (i.e. for parcellation use "NearestNeighbor")
        dim_add :
            Specify if the dimension of the input file is different than the 
            dimension usd to create the warp (i.e. " -e 3" fro 4 dim images)
    '''
    if not os.path.exists(warp_folder):
        os.makedirs(warp_folder)

    warp_file = os.path.join(warp_folder, warp_name)

    # Set the right path for the sh file --> to add to the readme
    antsRegistrationSyN_cmd = "bash ./tools/antsRegistrationSyN.sh -d 3 -r 2 -f " + \
        ref_file + " -m " + original_file + " -o " + warp_file
    
    if inv:
        complete_warp_file = warp_file + "1InverseWarp.nii.gz"
        affine_mat = "[" + warp_file + "0GenericAffine.mat, 1 ]"
    else:
        complete_warp_file = warp_file + "1Warp.nii.gz"
        affine_mat = warp_file + "0GenericAffine.mat"

    antsApplyTransforms_cmd = "antsApplyTransforms -d 3" + dim_add + " -t " + complete_warp_file + " -t " + \
        affine_mat + " -r " + ref_file + " -i " + input_file + " -o " + output_file + " -n " + interp
    
    mri_binarize_cmd = "mri_binarize --i " + output_file + " --o " + output_file + " --min 0.00001"

    if os.path.exists(warp_file + "1Warp.nii.gz"):
        print("antsRegistrationSyN already run")
    else:
        print(antsRegistrationSyN_cmd)
        logging.info('antsRegistrationSyN command: "{0}".'.format(antsRegistrationSyN_cmd))
        subprocess.call(antsRegistrationSyN_cmd, shell=True)

    
    if os.path.isfile(output_file):
        print("File already registered")
    else:
        print(antsApplyTransforms_cmd)
        logging.info('c command: "{0}".'.format(antsApplyTransforms_cmd))
        subprocess.call(antsApplyTransforms_cmd, shell=True)


    return
