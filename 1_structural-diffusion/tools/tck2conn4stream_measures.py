#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division

import argparse
from copy import copy
from genericpath import isfile
import logging
import os
import nibabel as nib
import numpy as np
import pandas as pd
from pandas import array
from turtle import isvisible
import itertools
from datetime import datetime
import json
import subprocess
import time
import csv

# to uncommant if want visualization and graph
#import matplotlib.pyplot as plt
#import matplotlib as mpl


def extract_weights_sum(weights_file):
    # opening the file in read mode
    weights_to_read = open(weights_file, "r")
                    
    # reading the file
    weights = weights_to_read.read()
                    
    # replacing end of line('/n') with ' ' and
    # splitting the text it further when '.' is seen.
    weights_list_tmp = weights.split(" ")
                    
    # close the file
    weights_to_read.close()
    weights_list = [float(i) for i in weights_list_tmp[:-1]]

    # sum 
    weights_stream_sum = np.sum(weights_list)

    return weights_stream_sum

def get_sum_of_weights(subj:str, sess:str, data_path:str, tracts:list, current_folder:str, isVerbose:bool, isForce:bool):
    rois = ['Loc_NA_Postcentral_L', 'Loc_NA_Cerebellum', 'Thal_IL_R', 'Precentral_L', 'Supp_Motor_Area_R']
    striat = ['v_d_Ca_L', 'v_d_Ca_R', 'vm_dl_PU_L', 'vm_dl_PU_R']
    tot= striat + rois

    # Define streamlines file
    tract_folder_path = os.path.join(current_folder, "tracts_tckedit")
    
    sum_of_weight = np.zeros([9,9])
        
    for tract in tracts :
        weights_file = os.path.join(tract_folder_path, subj + "_" + sess + '_' + tract[0] + '-'+ tract[1] + '_sift2.txt')
        if isfile(weights_file):
            i = tot.index(tract[0])
            j = tot.index(tract[1])

            sum_of_weight[i,j] += extract_weights_sum(weights_file)

        else : 
            print('problem of file location')

    output_path = current_folder
    df = pd.DataFrame(sum_of_weight, columns=tot)
    df.to_csv(os.path.join(output_path, subj + "_" + sess + '_' + 'metrics_sumofweights.csv'), index = False)

def get_mean_metric_over_streamlines(streamlines, weights_list):
    metric_values = []
    count_stream = len(streamlines)
    count_weights = len(weights_list)

    if count_stream != count_weights:
        raise ValueError("Streamlines and weights are not of the same length")
        
    weights_stream_sum = np.sum(weights_list)
    weights_stream_avg = np.sum(weights_list) / count_stream

    return weights_stream_sum, weights_stream_avg

def extract_stream_metrics(subjects:list, session:list, data_path:str, tracts:list, current_folder:str, isVerbose:bool, isForce:bool):
    ''' Copy files from original processing path to modeling path to 
    create a folder with all rquired data for modeling
    
        Parameters
        ----------
        subjects :
            Current subject
        session :
            Current session
        data_path :
            Path containing folder with all data
        isVerbose :
            Boolean indicating if descrption of running commands should be displayed
        isForce :
            Boolean indicating if files have to be overwritten
    '''
    
    subj_vec = []
    sess_vec = []
    tracts_vec = []
    weights_stream_sum_vec = []
    weights_stream_avg_vec = []
    FA_tcksample_vec = []
    for subj, sess in itertools.product(subjects, session):
        if isVerbose:
            logging.info('Extracting measures over tracts in subject {0}'.format(subj))

        logging.info('Subj: "{0}" - Sess: "{1}"'.format(subj, sess))

        proc_folder = os.path.join(data_path, 'derivatives', '01_dwi', subj, sess, 'dwi', 'proc')
        tract_folder_path = os.path.join(current_folder, "tracts_tckedit")
        FAMaps_file = os.path.join(proc_folder, subj + "_" + sess + "_dwi_FA.nii.gz")
        
        output_path = os.path.join(current_folder)

        if not os.path.exists(output_path) : 
            os.makedirs(output_path)
    

        if isfile(FAMaps_file):
            for tract in tracts:

                # Define streamlines file
                stramlines_scanner_file = os.path.join(tract_folder_path, \
                    subj + "_" + sess + '_' + tract + '.tck')

                if isfile(stramlines_scanner_file):
                    # Define weights file
                    weights_file = os.path.join(tract_folder_path, \
                        subj + "_" + sess + '_' + tract + '_sift2.txt')

                    # opening the file in read mode
                    weights_to_read = open(weights_file, "r")
                    
                    # reading the file
                    weights = weights_to_read.read()
                    
                    # replacing end of line('/n') with ' ' and
                    # splitting the text it further when '.' is seen.
                    weights_list_tmp = weights.split(" ")
                    
                    # close the file
                    weights_to_read.close()
                    weights_list = [float(i) for i in weights_list_tmp[:-1]]

                    # Define output FA file
                    FA_val = os.path.join(output_path, 'FA_csv', \
                        subj + "_" + sess + '_' + tract + '_FA.csv')
                    
                    if not os.path.exists(os.path.join(output_path, 'FA_csv')) :
                        os.makedirs(os.path.join(output_path, 'FA_csv'))
                    
                    weights_stream_sum = 0
                    weights_stream_avg = 0
                
                    tcksample_cmd = "tcksample " + stramlines_scanner_file + " " + FAMaps_file + \
                        " " + FA_val + " -stat_tck mean -force"
                    logging.info('tcksample command: "{0}".'.format(tcksample_cmd))
                    print(tcksample_cmd)
                    subprocess.call(tcksample_cmd, shell=True)

                    with open(FA_val, newline='') as csvfile:
                        try :
                            FA_df = list(csv.reader(csvfile)) # PB HERE error fail when too havy
                        except Exception as e:
                            print("The error is: ",e, 'on tract ', tract)
                            FA_df = 'NaN'
                            pass  

                    if FA_df :
                        if FA_df == 'NaN':
                            weights_stream_sum = 'out of bound'
                            weights_stream_avg = 'out of bound'
                            FA_tcksample = 'out of bound'
                        else : 
                            # reading the file                    
                            FA_df = FA_df[0][0].split()
                            FA_tcksample = np.mean([float(i) for i in FA_df])
                            streamlines_tck = nib.streamlines.load(stramlines_scanner_file)

                            #if len([float(i) for i in FA_df[1]]) != len(streamlines_tck.streamlines):
                            if len([float(i) for i in FA_df]) != len(streamlines_tck.streamlines):
                                raise ValueError("Streamlines and extracted FA not of the same size")
                                
                            logging.info('# of streamlines: "{0}".'.format(len(streamlines_tck.streamlines)))
                            
                            # Derive FA values along the streamlines by using tck file
                            [weights_stream_sum, weights_stream_avg] = get_mean_metric_over_streamlines(streamlines_tck.streamlines, \
                                weights_list)
                    else:
                        # printstreamlines_tck + ( "does not exist")
                        weights_stream_sum = 0
                        weights_stream_avg = 0
                        FA_tcksample = 0

                else:
                   # printstreamlines_tck + ( "does not exist")
                    weights_stream_sum = 0
                    weights_stream_avg = 0
                    FA_tcksample = 0
                
                subj_vec.append(subj)
                sess_vec.append(sess)
                tracts_vec.append(tract)
                weights_stream_sum_vec.append(weights_stream_sum)
                weights_stream_avg_vec.append(weights_stream_avg)
                FA_tcksample_vec.append(FA_tcksample)


        else:
            print('Subject {0}, session {1}'.format(subj, sess))

    data = {'subj': subj_vec,
    'sess': sess_vec,
    'tract': tracts_vec,
    'weights_stream_sum': weights_stream_sum_vec,
    'weights_stream_avg': weights_stream_avg_vec,
    'FA_tcksample_means': FA_tcksample_vec}

    df = pd.DataFrame(data, columns=['subj', 'sess', 'tract', 'weights_stream_sum', 'weights_stream_avg', 'FA_tcksample_means'])
    df.to_csv(os.path.join(output_path, subj + "_" + sess + '_' + 'metrics.csv'), index = False)

    return  

### ---------------test tck2connectom-----------------
def tckfortract(subj:str, sess:str, data_path:str, tracts:list, current_folder:str, isVerbose:bool, isForce:bool):
    for tract in tracts:
        
        # Input 
        tract_folder = os.path.join(current_folder, "tracts_tckedit")

        parcellation_file = os.path.join(current_folder, 'masks', subj + "_" + sess + '_global_mask.nii.gz')
        sift_file = os.path.join(tract_folder, subj + "_" + sess + '_' + tract[0] + '-'+ tract[1] + '_sift2.txt')
        tck_file = os.path.join(tract_folder, subj + "_" + sess + '_' + tract[0] + '-'+ tract[1] + '.tck')
        
        # Ouputs 
        connectome = os.path.join(current_folder,'test_tck2connectom', subj + "_" + sess + '_' + tract[0] + '-'+ tract[1] + "_tck2connectom_matrix.csv")
        
        if not os.path.exists(os.path.join(current_folder,'test_tck2connectom')):
            os.mkdir(os.path.join(current_folder,'test_tck2connectom'))
        
        if os.path.isfile(connectome) and not isForce: 
            logging.info('connectom already done: "{0}".'.format(connectome))
        else:
            tck2connectome_cmd = 'tck2connectome ' + tck_file + ' ' + parcellation_file+ ' ' + connectome + ' -tck_weights_in ' +  sift_file + ' -force' 

            logging.info('tck2connectome command: "{0}".'.format(tck2connectome_cmd))
            subprocess.call(tck2connectome_cmd, shell=True)
        
### ---------------visualization-----------------
# Don't work because don't have Matplotlib on the server - Run them locally! 

def connectivity_matrix(df, lim_sup=6000, title = 'connectivity matrix') : 
    
    connectivity = np.zeros([9,9])

    for i, row in enumerate(df.to_numpy()):
        if isinstance(row[0], str): row = row[0].split(' ')
            
        for j, el in enumerate(row): 
            connectivity[i,j] = el

    #load the labels
    df_lab = pd.read_csv('sub-51T01_ses-baseline_global_mask.csv')

    labels = [lab[1] for lab in df_lab.to_numpy()]


    # To saturate the heatmap
    def zlim(x,lim): 
        x[x < min(lim)] = min(lim)
        x[x > max(lim)] = max(lim)
        return(x)

    lim = (0,lim_sup)

    fig, ax = plt.subplots(figsize=(7, 7))
    im = ax.imshow(zlim(connectivity, lim))

    ax.set_xticks(np.arange(len(labels)), labels=labels)
    ax.set_yticks(np.arange(len(labels)), labels=labels)
    
    fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(0, lim_sup)),
             ax=ax, orientation='vertical', location= 'right', shrink=0.6)
    
    for i in range(len(labels)):
        for j in range(len(labels)):
            text = ax.text(j, i, "{:.1f}".format(connectivity[i, j]), ha="center", va="center", color="w")


    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")
    plt.setp(ax.get_yticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")


    ax.set_title(title)
    fig.tight_layout()
    plt.show()

def print_connect_matrix(subj:str, sess:str, data_path:str, tracts:list, isVerbose:bool, isForce:bool):
    file= os.path.join(data_path, '01_tracts', subj, sess,'roi2roi','fMRI_study', subj + "_" + sess + '_connect_matrix.csv')
    df = pd.read_csv(file)

    connectivity_matrix(df, lim_sup = 5000, title = f'Connectivity matrix: {subj}, {sess}')
    



