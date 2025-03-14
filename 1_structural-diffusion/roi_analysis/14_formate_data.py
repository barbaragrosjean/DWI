# The file is used to formate the data used for the statistical analysis analysis. it manipulate csv file. You will find :
# - formate_seed_based to formate the csv from the seed based that contains one unique value into a matrix.
# - formate_behav to formate the behaviral gain into one vector column.
# - formate_roi2roi_Pu and formate_roi2roi_Ca to have the connectivity metric of the tract that link the putamen right and the caudate respectively to the rest of the network.
# - formate_roi2roi_Pu_network to have the connectivity metric for the network selected included the putamen.

# Depending call the one needed in the main.

import numpy as np
import tools
import os
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

    
def formate_seed_based(subjects:list, sess:str, data_path:str, isForce:bool): 

    print('Formating the seed based data for analysis ...')
    output = os.path.join(data_path, 'derivatives', '01_analysis', 'seed_metric_df.csv')

    if os.path.isfile(output) and not isForce :
        print('Formating seed based data already done')
        return 
    else :
        # Load the seed-based metrics
        seeds = ['v_d_Ca_L', 'v_d_Ca_R', 'vm_dl_PU_L', 'vm_dl_PU_R']
        seed_metric = np.zeros([int(len(subjects)), len(seeds)])
        for i, sub in enumerate(subjects):
            folder_path = os.path.join(data_path, 'derivatives', '01_tracts', sub, 'ses-baseline', 'striat') 
            for j,s in enumerate(seeds):
                file_name=os.path.join(folder_path, sub + "_ses-baseline_"+ s +"_metric.csv")
                df = pd.read_csv(file_name)  
                seed_metric[i, j]=float(df.columns[0])

        dict_= {}
        for i,s in enumerate(seeds) : 
            dict_[s] = seed_metric[:,i]

        df = pd.DataFrame(dict_)
        df.to_csv(output, index=False)

def formate_behav(subjects:list, data_path:str, isForce:bool):
    # Load behavioral data
    file_name=os.path.join(data_path, 'behav_output', 'gain_blinded.csv')
    
    if os.path.isfile(file_name) : 
        df_behav = pd.read_csv(file_name) 

        output = os.path.join(data_path, 'derivatives', '01_analysis','behav.csv')
        if os.path.isfile(output) and not isForce :
            print('Formating behavioral data already done')
            return 
        else :
            # remove the subject if it's ont in the list
            for sub_id in df_behav['CODE'] : 
                sub_name = 'sub-' + sub_id
                if not sub_name in subjects : 
                    df_behav = df_behav.drop(df_behav[df_behav['CODE'] == sub_id].index)

            # Extract the gain and the condition
            df_gain_1 = df_behav[df_behav['CONDITION'] == 1]
            gain_1 = np.array(df_gain_1['gain'])


            df_gain_2 = df_behav[df_behav['CONDITION'] == 2]
            gain_2 = np.array(df_gain_2['gain'])


            gain = gain_2 - gain_1 
            ID_subj = df_behav['CODE'][:14]

            df_b = pd.DataFrame()
            df_b['gain'] = gain
            
            df_b.to_csv(output, index=False)
    else : 
        raise FileNotFoundError(file_name + 'not existing')
        return 
    
    
def formate_roi2roi_Pu(subjects:list, sess:str, datapath:str, isForce:bool):
    print('Formating data for Putamen tract analysis ...')
    output = os.path.join(data_path, 'derivatives', '01_analysis','Pu.csv')

    label_file = os.path.join(data_path,'sub-51T01', 'ses-baseline','roi2roi', 'fMRI_study','masks', 'sub-51T01_ses-baseline_global_mask.csv')
    
    if not os.path.isfile(label_file) : 
        print('Should provide a label file: ', label_file)
        return 

    labels_df = pd.read_csv(label_file)
    labels = np.array(labels_df['roi'])

    to_remove = ['v_d_Ca_L', 'v_d_Ca_R', 'vm_dl_PU_L']

    df_con = pd.DataFrame()
    for k, sub in enumerate(subjects): 
        folder_path = os.path.join(data_path,'derivatives', '01_tracts', sub, "ses-baseline",'roi2roi', 'fMRI_study')
        file_name=os.path.join(folder_path, sub + "_ses-baseline_connect_matrix.csv")
        if os.path.isfile(file_name) :
            if os.path.isfile(output) and not isForce :
                print('Formating data for putamen tracts analysis already done')
            else : 
                df = pd.read_csv(file_name)
                suj = {}
                for i, row in enumerate(df.to_numpy()):
                    row = row[0].split(' ')
                    for j, el in enumerate(row): 
                        if labels[i] == 'vm_dl_PU_R' or labels[j]== 'vm_dl_PU_R': 
                            # Remove inter striatum connection
                            if not (labels[i] in to_remove or labels[j] in to_remove) : 
                                suj[f'{labels[i]}-{labels[j]}'] = el
                            
                if k == 0: df_con = pd.DataFrame(columns = suj.keys())
                df_con.loc[k] = suj
        else : 
            raise FileNotFoundError(file_name + 'not existing')
            return 

    #remove 0 column : 
    for col in df_con.columns :
        if not np.any(df_con[col]) : 
            df_con = df_con.drop(columns=col)
        
    df_con.to_csv(output, index=False)

def formate_roi2roi_Ca(subjects:list, sess:str, datapath:str, isForce:bool):
    print('Formating data for Putamen tract analysis ...')
    output = os.path.join(data_path, 'derivatives', '01_analysis','Ca.csv')

    # Load the roi2roi metrics
    label_file = os.path.join(data_path,'sub-51T01', 'ses-baseline','roi2roi', 'fMRI_study','masks', 'sub-51T01_ses-baseline_global_mask.csv')
    
    if not os.path.isfile(label_file) : 
        print('Should provide a label file: ', label_file)
        return 

    labels_df = pd.read_csv(label_file)
    labels = np.array(labels_df['roi'])
    to_remove = ['v_d_Ca_L','vm_dl_PU_L','vm_dl_PU_R']

    df_con = pd.DataFrame()
    for k, sub in enumerate(subjects): 
        folder_path = os.path.join(data_path,'derivatives', '01_tracts', sub, "ses-baseline",'roi2roi', 'fMRI_study')
        file_name=os.path.join(folder_path, sub + "_ses-baseline_connect_matrix.csv")
        if os.path.isfile(file_name) :
            if os.path.isfile(output) and not isForce :
                print('Formating data for caudate tracts analysis already done')
            else : 
                df = pd.read_csv(file_name)
                suj = {}
                for i, row in enumerate(df.to_numpy()):
                    row = row[0].split(' ')
                    for j, el in enumerate(row): 
                        if labels[i] == 'v_d_Ca_R' or labels[j]== 'v_d_Ca_R': 
                            # Remove inter striatum connection
                            if not (labels[i] in to_remove or labels[j] in to_remove) : 
                                suj[f'{labels[i]}-{labels[j]}'] = el
                            
                if k == 0: df_con = pd.DataFrame(columns = suj.keys())
                df_con.loc[k] = suj
        else : 
            raise FileNotFoundError(file_name + 'not existing')
            return 

    #remove 0 column : 
    for col in df_con.columns :
        if not np.any(df_con[col]) : 
            df_con = df_con.drop(columns=col)
        

    df_con.to_csv(output, index=False)

def formate_roi2roi_Pu_net(subjects:list, sess:str, datapath:str, isForce:bool):
    print('Formating data for Putamen tract analysis ...')
    output = os.path.join(data_path, 'derivatives', '01_analysis','Pu_network.csv')
    
    # Load the roi2roi metrics
    label_file = os.path.join(data_path,'sub-51T01', 'ses-baseline','roi2roi', 'fMRI_study','masks', 'sub-51T01_ses-baseline_global_mask.csv')
    
    if not os.path.isfile(label_file) : 
        print('Should provide a label file: ', label_file)
        return 

    labels_df = pd.read_csv(label_file)
    labels = np.array(labels_df['roi'])
    to_remove = ['v_d_Ca_L', 'v_d_Ca_R', 'vm_dl_PU_L']

    df_con = pd.DataFrame()
    for k, sub in enumerate(subjects): 
        folder_path = os.path.join(data_path,'derivatives', '01_tracts', sub, "ses-baseline",'roi2roi', 'fMRI_study')
        file_name=os.path.join(folder_path, sub + "_ses-baseline_connect_matrix.csv")

        if os.path.isfile(file_name) :
            if os.path.isfile(output) and not isForce :
                print('Formating data for putamen network analysis already done')
                 
            else : 
                df = pd.read_csv(file_name)
                suj = {}
                for i, row in enumerate(df.to_numpy()):
                    row = row[0].split(' ')
                    for j, el in enumerate(row): 
                        # Remove inter striatum connection
                        if not (labels[i] in to_remove or labels[j] in to_remove) : 
                            suj[f'{labels[i]}-{labels[j]}'] = el
        
                if k == 0: df_con = pd.DataFrame(columns = suj.keys())
                df_con.loc[k] = suj
                
        else : 
            raise FileNotFoundError(file_name + 'not existing')
            return 

    #remove 0 column : 
    for col in df_con.columns :
        if not np.any(df_con[col]) : 
            df_con = df_con.drop(columns=col)
                
    df_con.to_csv(output, index=False)
            
            

if __name__ == "__main__":
    print('Formating data for statistical analysis ... ')

    parser = buildArgsParser()
    args = parser.parse_args()

    isForce = args.isForce

    if args.isVerbose:
        logging.basicConfig(level=logging.DEBUG)
    
    subj_list = [subj for subj in args.subj]

    data_path = args.data_path

    if "all" in subj_list:
        subjects = [s for s in os.listdir(data_path) if os.path.isdir(os.path.join(data_path,'derivatives', '01_dwi', s)) if "sub-51T" in s]
    else:
        subjects = ['sub-' + subj for subj in subj_list]
    
    
    sess_list = [sess for sess in args.sess]
    if "all" in sess_list:
       sessions = ["ses-T1", "ses-T2", "ses-T3", "ses-T4"]
    else:
        sessions = ['ses-' + sess for sess in sess_list]
    
    date = datetime.now()
    formatted_datetime = date.strftime("%Y-%m-%d-%H-%M-%S")
    fail_list_filename = f"fail_list_13_seed_based{formatted_datetime}.txt"
    for subj, sess in itertools.product(subjects, sessions):        
            try:
                formate_seed_based(subj_list, sess_list[0], data_path, isForce)
                formate_roi2roi_Pu(subj_list,sess_list[0], data_path, isForce)
                formate_roi2roi_Ca(subj_list,sess_list[0], data_path, isForce)
                formate_roi2roi_Pu_net(subj_list,sess_list[0], data_path, isForce)
                
            except Exception as e:
                with open(fail_list_filename, "+a") as f:
                    f.write(f"{subj} {sess} \n")
                    f.write(f"{str(e)} \n")
    