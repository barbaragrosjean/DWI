#!/bin/bash

SCRIPT=$(basename $0)


data_path='/data/PlasMA/wp_51T' 
bids_path='/mnt/Hummel-Data/TI/51T' #NOTUSED 
local_path='/data/PlasMA/wp_51T' # local = server uphummelsrv2 data/PlasMA/wp_51T
server_path='/home/bgrosjea/mnt/Hummel-Data/TI/mri/51T' # Data server_path from my home to right on it

#
lesion=false
MNI=false
save_to_server=true

#      0 1 2 3 4 5 6 7 8 9 10
works=(0 0 0 0 0 0 0 0 0 0 0)
# 0: 01_copy_file_locally
# 1: 02_dwi_preprocessing 
# 2: 03_lesionTransplantation_anat 
# 3: 04_freesurfer 
# 4: 05_anat_registration_dwi, 05_lesion_registration, 05_T1w2MNI_reg 
# 5: 06_dwi_processing, 06_compute_scalar_maps.py

## ROI Analysis 
# 6: roi_analysis/11_register_rois_MNI2B0.py
# 7: roi_analysis/12_create_parc.py
# 8: roi_analysis/13_dwi_extract_tracts_tckedit.py 
# 9: roi_analysis/13_seed_based.py
# 10: roi_analysis/GLM.R 


# Transform long options to short ones
for arg in "$@"; do
  shift
  case "$arg" in
    "--subj") set -- "$@" "-s" ;;
    "--sess") set -- "$@" "-t" ;;
    "--data_path")   set -- "$@" "-i" ;;
    "--bids_path")   set -- "$@" "-b" ;;
    "--server_path")   set -- "$@" "-o" ;;
    "--local_path")   set -- "$@" "-p" ;;
    "--lesion")   set -- "$@" "-l" ;;
    "--works")   set -- "$@" "-w" ;;
    *)        set -- "$@" "$arg"
  esac
done

while getopts ":s:t:i:b:o:p:l:w:" arg ; do
    case ${arg} in
        s ) set -f # disable glob
            IFS=' ' # split on space characters
            subj=($OPTARG) # use the split+glob operator
            ;;
        t ) set -f # disable glob
            IFS=' ' # split on space characters
            sess=($OPTARG)
            ;;
        i )
            data_path=${OPTARG}
            ;;
        b )
            bids_path=${OPTARG}
            ;;
        o )
            server_path=${OPTARG}
            ;;
        p )
            local_path=${OPTARG}
            ;;
        l )
            lesion=${OPTARG}
            ;;
        w )
            works=${OPTARG}
            ;;
        # Option entered by user is invalid entry
        \? )
            echo "${SCRIPT}: Invalid option --${OPTARG} ignored"
            ;;
        # Option entered by user is missing the option parameter
        : )
            echo "${SCRIPT}: must supply a parameter to --${OPTARG}"
            ;;
    esac
done

usage() {
    echo "Usage: $1"
    echo "subj subject id (i.e. "103")"
    echo "sess session (i.e. "T1")"
    echo "[data_path] folder where data are stored (i.e. /mnt/hummelarch/AVANCER/Data/raw/stroke/)"
    echo "[bids_path] folder where bids format need to be stored (i.e. /mnt/hummelarch/AVANCER/Data/analysis/stroke/group_analysis/mri/sourcedata/BIDS)"
    echo "[server_path] folder data need to be stored on the server (i.e. /mnt/hummelarch/AVANCER/Data/analysis/stroke/group_analysis/mri/)"
    echo "[local_path] folder where bids format need to be copied locally (i.e. /home/elena/data/subjects/patients/AVANCER)"
    echo "[works] Which steps need to be performed"
}

if [[ $# -eq 0 ]]; then
    usage ${SCRIPT}
    exit 1
fi


for sbj in "${subj[@]}"
do
    echo "$sbj"
    for ses in "${sess[@]}"

    do
        echo "$ses"
        # run DICOM 2 BIDS - add pat to Julia's script

        # If the session folder doesn't exist creat it on local 
        if [ ! -d "${local_path}/derivatives/01_dwi/sub-$sbj/ses-$ses" ]; then
            if [ ! "$sbj" == "all" ]; then
                echo "mkdir -p  $local_path/derivatives/01_dwi/sub-$sbj/ses-$ses"
                mkdir -p  $local_path/derivatives/01_dwi/sub-$sbj/ses-$ses 
                chmod -R 755 $local_path/derivatives/01_dwi/sub-$sbj/ses-$ses 
            fi
        fi

        # If the data are not on the running server copy them from the data server 
        if [ -d "$server_path/sub-$sbj/ses-$ses/dwi" ]; then
            if [ ! -d "${local_path}/sub-$sbj/ses-$ses" ]; then
                echo "mkdir -p  $local_path/sub-$sbj/ses-$ses"
                mkdir -p  $local_path/sub-$sbj/ses-$ses 
                chmod -R 755 $local_path/sub-$sbj/ses-$ses 
                echo "cp -r $server_path/sub-$sbj/ses-$ses/. $local_path/sub-$sbj/ses-$ses"
                cp -r $server_path/sub-$sbj/ses-$ses/. $local_path/sub-$sbj/ses-$ses 
            fi
            # If the lesion data are not on the runnning server copy them from the data server
            if $lesion ; then
                if [ ! -f "${local_path}/sub-AVC$sbj/ses-$ses/anat/sub-AVC${sbj}_ses-${ses}_acq-mprage_T1w_label-lesion_roi.nii.gz" ]; then
                    echo "Please draw lesion for subject sub-${sbj} session ${ses}"
                else
                    if [[ ${works[2]} -eq 1 ]]; then
                        python processing_scripts/diffusion/3_anatomical_lesion_transplantation.py --subjects $sbj --sessions $ses --data_path ${local_path} --output_path ${local_path}/derivatives/ -v
                    fi
                fi
            fi
        fi

        #------------------------------------------------------------#
        #### PREPROCESSING #### 
        #------------------------------------------------------------#
 
        # Step 4: free surfer 
        if [[ ${works[3]} -eq 1 ]]; then
            # If folder freesurfer doesn't existe create it on local 
            if [ ! -d "${local_path}/derivatives/01_freesurfer/sub-$subj-ses-$sess" ]; then
                echo "mkdir -p  ${local_path}/derivatives/01_freesurfer/sub-$subj-ses-$sess"
                mkdir -p  $local_path/derivatives/01_freesurfer/sub-$subj-ses-$sess
                chmod -R 755 $local_path/derivatives/01_freesurfer/sub-$subj-ses-$sess
            fi
            
            echo "python 04_freesurfer.py --subj $sbj --sess $ses --data_path ${local_path} -v" 
            python 04_freesurfer.py --subj $sbj --sess $ses --data_path ${local_path} -v
    
        fi
        
        # If it does't already existe in the data server copy the 01_freesurfer result to server
        if $save_to_server ; then 
            if [ ! -d "${server_path}/derivatives/01_freesurfer/sub-${sbj}_ses-${ses}" ]; then
                if [ -d "${local_path}/derivatives/01_freesurfer/sub-${sbj}_ses-${ses}" ]; then
                    mkdir -p  $server_path/derivatives/01_freesurfer/sub-${sbj}_ses-${ses}
                    cp -r $local_path/derivatives/01_freesurfer/sub-${sbj}_ses-${ses}. $server_path/01_freesurfer/sub-${sbj}_ses-${ses}
                fi
            fi
        fi

        # Check if the data are aviable on the running or data server under dwi folder
        if [ -d "$local_path/sub-$sbj/ses-$ses/dwi" ] || [ -d "$local_path/sub-$sbj/ses-$ses/dwi" ]; then
            echo "Data found"

            # Step 2: DWI preprossessing 
            if [[ ${works[1]} -eq 1 ]]; then
                # If folder  preproc doesn't exist create it on local 
                if [ ! -d "${local_path}/derivatives/01_dwi/sub-$sbj/ses-$ses/dwi/preproc" ]; then
                    echo "mkdir -p  $local_path/derivatives/01_dwi/sub-$sbj/ses-$ses/dwi/preproc"
                    mkdir -p  $local_path/derivatives/01_dwi/sub-$sbj/ses-$ses/dwi/preproc
                    chmod -R 755 $local_path/derivatives/01_dwi/sub-$sbj/ses-$ses/dwi/preproc
                fi
            
                echo "python 02_dwi_preprocessing.py --subj $sbj --sess $ses --data_path ${local_path}"
                python 02_dwi_preprocessing.py --subj $sbj --sess $ses --data_path ${local_path}
                
            fi

            # Step 5: DWI data registration  
            if [[ ${works[4]} -eq 1 ]]; then
                # Anat registration 
                echo "05_anat_registration_dwi.py --subj $sbj --sess $ses --data_path ${local_path}" 
                python 05_anat_registration_dwi.py --subj $sbj --sess $ses --data_path ${local_path} 
                    
                # Lesion registration 
                if $lesion ; then 
                    echo "python 05_lesion_registration.py --subjects $sbj --sessions $ses --data_path ${local_path} \
                    --output_path ${local_path}/derivatives/01_dwi/ -v" 

                    python 05_lesion_registration.py --subjects $sbj --sessions $ses --data_path ${local_path} \
                    --output_path ${local_path}/derivatives/01_dwi/ -v 
                        
                fi

                # T1 to MNI: universal template - to ignore
                if $MNI ; then 
                    echo "python 05_T1w2MNI.py --subjects $sbj --sessions $ses --data_path ${local_path} \
                    --output_path ${local_path}/derivatives/01_dwi/preproc -v" 

                    python 05_T1w2MNI.py --subjects $sbj --sessions $ses --data_path ${local_path} \
                    --output_path ${local_path}/derivatives/01_dwi/preproc -v 
                fi
            fi

            #------------------------------------------------------------#
            #### PROCESSING #### 
            #------------------------------------------------------------#

            # Step 6: Processing dwi : tractography algo
            if [[ ${works[5]} -eq 1 ]]; then

                # If folder proc doesn't exist create it on local 
                if [ ! -d "${local_path}/derivatives/01_dwi/sub-$sbj/ses-$ses/dwi/proc" ]; then
                    echo "mkdir -p  $local_path/derivatives/01_dwi/sub-$sbj/ses-$ses/dwi/proc"  
                    mkdir -p  $local_path/derivatives/01_dwi/sub-$sbj/ses-$ses/dwi/proc
                    chmod -R 755 $local_path/derivatives/01_dwi/sub-$sbj/ses-$ses/dwi/proc
                fi
     
                echo "python 06_dwi_processing.py --subj $sbj --sess $ses --data_path ${local_path}"
                python 06_dwi_processing.py --subj $sbj --sess $ses --data_path ${local_path} 

                echo "python 06_compute_scalar_maps.py --subj $sbj --sess $ses --data_path ${local_path}"
                python 06_compute_scalar_maps.py --subj $sbj --sess $ses --data_path ${local_path} 
            fi

            # If it does't already existe in the data server copy the 01_dwi result to server
            if $save_to_server ; then 
                if [ ! -d "${server_path}/derivatives/01_dwi/sub-${sbj}_ses-${ses}" ]; then
                    if [ -d "${local_path}/derivatives/01_dwi/sub-${sbj}_ses-${ses}" ]; then
                        mkdir -p  $server_path/derivatives/01_dwi/sub-${sbj}_ses-${ses}
                        cp -r $local_path/derivatives/01_dwi/sub-${sbj}_ses-${ses}. $server_path/01_dwi/sub-${sbj}_ses-${ses}
                    fi
                fi
            fi

        fi
    done
done

#--------------- ROI-2-ROI and SEED-BASED ANALYSIS -------------- 
# If all the previsous steps are done succefully. You can start the analysis
# this scrip provide you, 2 analysis 
# - Roi to roi analysis based voxel define by fMRI study from Wessel, Benaeto and al. (2023)
# - Seed based approach for seed placed in the striatum 
#------------------------------------------------------------#

for sbj in "${subj[@]}"
do
    echo "$sbj"
    for ses in "${sess[@]}"
    do
        echo "$ses"
        if [[ ${works[6]} -eq 1 ]]; then
            # Registration the roi to dwi space
            echo "python roi_analysis/11_register_rois_MNI2B0.py --subj $sbj --sess $ses --data_path ${local_path}"
            python roi_analysis/11_register_rois_MNI2B0.py --subj $sbj --sess $ses --data_path ${local_path} 
        fi

        if [[ ${works[7]} -eq 1 ]]; then
            # create mask with roi 
            echo "python roi_analysis/12_create_parc.py --subj $sbj --sess $ses --data_path ${local_path}"
            python roi_analysis/12_create_parc.py --subj $sbj --sess $ses --data_path ${local_path} 
            
        fi 

        if [[ ${works[8]} -eq 1 ]]; then
            # tract extraction for roi2roi analysis
            echo "python roi_analysis/13_dwi_extract_tracts_tckedit.py --subj $sbj --sess $ses --data_path ${local_path}"
            python roi_analysis/13_dwi_extract_tracts_tckedit.py --subj $sbj --sess $ses --data_path ${local_path}  
        fi

        if [[ ${works[9]} -eq 1 ]]; then
            # tract extraction for seed based analysis 
            echo "python roi_analysis/13_seed_based.py --subj $sbj --sess $ses --data_path ${local_path}"
            python roi_analysis/13_seed_based.py --subj $sbj --sess $ses --data_path ${local_path} 
        fi 

        if [[ ${works[10]} -eq 1 ]]; then
            # Before starting the statistic analysis need to formate well the data
            echo "python roi_analysis/14_formate_data.py --subj $sbj --sess $ses --data_path ${local_path}"
            python roi_analysis/14_formate_data.py --subj $sbj --sess $ses --data_path ${local_path}
             
            # Rune the R file to get the statistic
            echo "Rscript roi_analysis/14_GLM.R"
            Rscript roi_analysis/14_GLM.R "$sbj" "$ses" "$local_path"

        fi 

    # If it does't already existe in the data server copy the proc result from the running server
    if $save_to_server ; then 
        if [ ! -d "${server_path}/derivatives/01_tracts/sub-${sbj}_ses-${ses}" ]; then
            if [ -d "${local_path}/derivatives/01_tracts/sub-${sbj}_ses-${ses}" ]; then
                mkdir -p  $server_path/derivatives/01_tracts/sub-${sbj}_ses-${ses}
                cp -r $local_path/derivatives/01_tracts/sub-${sbj}_ses-${ses}. $server_path/derivatives/01_tckedit/sub-${sbj}_ses-${ses}
            fi
        fi

        if [ ! -d "${server_path}/derivatives/01_mni/sub-${sbj}_ses-${ses}" ]; then
            if [ -d "${local_path}/derivatives/01_mni/sub-${sbj}_ses-${ses}" ]; then
                mkdir -p  $server_path/derivatives/01_mni/sub-${sbj}_ses-${ses}
                cp -r $local_path/derivatives/01_mni/sub-${sbj}_ses-${ses}. $server_path/derivatives/01_mni/sub-${sbj}_ses-${ses}
            fi
        fi
    fi


    done
done