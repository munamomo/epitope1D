# *************************************************************************************
# *               University of Melbourne
# *   ----------------------------------------------
# * Bruna Moreira da Silva - bmoreiradasi@student.unimelb.edu.au
# *   ----------------------------------------------
# @article {da Silva2022,
# 	author = {da Silva, Bruna Moreira and Ascher, David B. and Pires, Douglas E. V.},
# 	title = {epitope1D: Accurate Taxonomy-Aware B-Cell Linear Epitope Prediction},
# 	journal={bioRxiv},
# 	year = {2022},
# }
# *************************************************************************************

import os
import sys
import re
import time
import argparse
import numpy as np
import pandas as pd
import pickle
import joblib
from epitope1d.codes.CTDC import *
from epitope1d.codes.pepCSM import *

def epitope_pred(info_pickle):

    basedir = os.path.abspath(os.path.dirname(__file__))
    CODE_FOLDER = os.path.join(basedir, 'epitope1d', 'codes')
    MODEL_FOLDER = os.path.join(basedir, 'epitope1d', 'models')

    infos = pickle.load(open(info_pickle, 'rb')) #Open the pickle
    fasta_file = infos['fasta_list'] #retrieve fasta in list format
    organism = [infos['organism']]  # retrieve organism taxonomy as list to subset the dataframe below
    window = infos['window'] #retrieve window size in integer format
    run_id = infos['run_id']
    output = infos['working_dir']

    print("\nepitope1D starting...Running id information: ",run_id)
    sliced_fasta = []

    for each in fasta_file:
        size_pep = len(each[1])
        full_seq = each[1].upper()  # pep seq vector will be from position 0 to (size_pep - 1)
        print("Full sequence: ", full_seq, " with size of: ", size_pep)
        start = 0

        ##Rules for window and peptide size:
        if size_pep < 6:
            print("Warning! Peptide size should be at minimum 6 aa length. Sequences smaller than that will be disregarded.")
            continue  # skip sequence.

        elif size_pep >= window:  ##Slice the protein or long peptides into sliding windows, keeping the same fasta header.
            total_pep = int(size_pep - window + 1)
            for i in range(total_pep):
                sliced_fasta.append([each[0], full_seq[start:start + window]])
                start += 1

        elif (size_pep < window):  # If the peptide size is smaller than the window, consider the entire sequence as it is.
            print("Warning! Peptide size smaller than the informed window length.")
            sliced_fasta.append([each[0], full_seq])

    data = pd.DataFrame(sliced_fasta, columns=["Fasta_header", "Peptide"])
    rows = len(data.index)

    # -----------------------START FEATURES CALCULATION---------------------#:
    print("\n Start Feature calculation...\n")

    ##(1) aat_antigenicity(inputfile) #The output is a dataframe
    output_file_aat = data.copy()
    aat_scales = CODE_FOLDER + '/aat_iedb_normal_may22.txt'
    aat_dict = {}

    with open(aat_scales) as f_aat:
        for line in f_aat:
            (key, val) = line.split()
            aat_dict[str(key)] = float(val)

    for index, row in data.iterrows():
        entry = str(row[1])
        size_entry_t = len(entry) - 2
        score_aat = 0
        score_aat_list = []

        for i in range(0, (len(entry) - 2)):
            try:
                slice_pep = str(entry[i: i + 3])
                score_aat += round(float(aat_dict[slice_pep]), 4)
                score_aat_list.append(score_aat)
            except KeyError:
                score_aat += float(-1)
                score_aat_list.append(score_aat)

        output_file_aat.loc[index, "AAT_score"] = round(float(score_aat / size_entry_t), 4)
        output_file_aat.loc[index, "AAT_max"] = round(max(score_aat_list), 4)
        output_file_aat.loc[index, "AAT_min"] = round(min(score_aat_list), 4)

    # pepCSM
    output_data_csm = data.copy()
    csm_feat = pepCSM(output_data_csm, '9', '1', '1', '4')  # dmax=9, dmin=1, cutoff_step=1, sign_type=4
    csm_feat_pd = pd.DataFrame(csm_feat)
    csm_feat_pd.columns = csm_feat_pd.iloc[0]
    csm_feat_pd = csm_feat_pd.iloc[1:].reset_index(drop=True)

    # CTDC
    ctdc_feat = CTDC(sliced_fasta)
    ctdc_feat_pd = pd.DataFrame(ctdc_feat)
    ctdc_feat_pd.columns = ctdc_feat_pd.iloc[0]
    ctdc_feat_pd = ctdc_feat_pd.iloc[1:].reset_index(drop=True)
    del ctdc_feat_pd["Peptides"]

    # Create Organims Data
    organism_list = ['Metamonada', 'Discoba', 'Sar', 'Viridiplantae',  'Opisthokonta', 'Terrabacteria group', 'Proteobacteria', 'PVC group',\
                     'Spirochaetes', 'FCB group', 'Thermodesulfobacteria', 'Fusobacteria', 'Riboviria', 'Duplodnaviria', 'Monodnaviria', \
                      'Varidnaviria', 'Ribozyviria', 'Anelloviridae', 'Naldaviricetes', 'Adnaviria']

    org_df = pd.DataFrame(np.zeros((rows, 20)), columns = organism_list)
    org_df[organism] = org_df[organism].replace(0, int(1))

    ###Merge features
    feat_result_final = pd.concat([output_file_aat, csm_feat_pd, ctdc_feat_pd, org_df], axis=1)

    # Replace infinite, nan, empty with 0.
    feat_result_id = feat_result_final.copy()
    feat_result_id = feat_result_id.replace([np.inf, -np.inf, np.nan], 0)
    feat_result_id.insert(0, "ID", value=np.arange(start=1, stop=(len(feat_result_id) + 1), step=1))

    # -----------------------START MODEL PREDICTION---------------------#:
    print("\n Start model prediction...\n")

    ## Prediction Model Random Forest
    model = '/epitope1d_RF_model.sav'
    feature_list_model = '/epitope1d_order_features.sav'

    loaded_model = joblib.load(MODEL_FOLDER + model)
    loaded_feature = joblib.load(MODEL_FOLDER + feature_list_model)
    final_temp = feat_result_id[loaded_feature]

    predicted_value = loaded_model.predict(final_temp) #This is the output array of prediction.
    prob_value = loaded_model.predict_proba(final_temp) #predict_proba() returns a numpy array: first column is the probability that target=0, second is target=1.

    #Output Dataframe
    output_pred = feat_result_id[['ID','Fasta_header','Peptide']]
    output_pred.insert(3, "Prediction", value=predicted_value) #Insert a new column named Prediction and place the predicted_value information.
    output_pred.insert(4, "Score_Epitope", value=prob_value[:,1])
    output_pred.insert(5, "Classification", value=output_pred['Prediction'].map({1: 'epitope', 0: 'non-epitope'}))

    output_pred.to_csv(os.path.join(output, "output-prediction.csv"), index=False)

    print("\n\n epitope1D finished. Result saved at: ",str(output))
    return True

def main(fasta,organism,window,output_dir):

    #Create random run_id identifier and folder to save results.
    run_id = str(time.time()).replace('.', '')
    cmd = "mkdir " + output_dir + "/" + str(run_id)
    os.system(cmd)
    folder = output_dir + "/" + str(run_id)

    if organism in ['Terrabacteria_group','PVC_group','FCB_group']:
        org_first,org_second = organism.split("_")
        organism = str(org_first) + " " + str(org_second)

    # Check FASTA format
    with open(fasta) as f:
        records = f.read()

    if re.search('>', records) == None:
        sys.exit('Input Error. Please provide input file in fasta format with corresponding headers '>'.')

    # Transform FASTA format into list
    records = records.split('>')[1:]
    myFasta = []

    for fastas in records:
        array = fastas.split('\n')
        name, sequence = array[0].split()[0], re.sub('[^ARNDCQEGHILKMFPSTWYV-]', '', ''.join(array[1:]).upper())
        myFasta.append([name, sequence])

    infos_pickle = os.path.join(folder, "infos.p")
    infos = {'fasta_list': myFasta, \
            'organism': organism, 'window':window, 'run_id': str(run_id), 'working_dir': folder}
    pickle.dump(infos, open(infos_pickle, 'wb'))

    epitope_pred(infos_pickle)

if __name__ == "__main__":
    output_dir = os.path.dirname(os.path.abspath(__file__))

    organism_list = ['Metamonada', 'Discoba', 'Sar', 'Viridiplantae', 'Opisthokonta', 'Terrabacteria_group',
                     'Proteobacteria', 'PVC_group', \
                     'Spirochaetes', 'FCB_group', 'Thermodesulfobacteria', 'Fusobacteria', 'Riboviria', 'Duplodnaviria',
                     'Monodnaviria', \
                     'Varidnaviria', 'Ribozyviria', 'Anelloviridae', 'Naldaviricetes', 'Adnaviria']

    # REQUIRED
    parser = argparse.ArgumentParser(description='usage: python run_epitope1d.py -fasta fastafile.fasta -organism Ribozyviria -window 25')
    parser.add_argument("-fasta", help="Path to the input fasta file. Sequences must be from the same Taxonomy.",required=True)
    parser.add_argument("-organism", help="Organism Taxonomy name. Only one per fasta file is accepted.", type = str, choices = organism_list,required=True)
    parser.add_argument("-window", help="Desired window size (amino acid length) to screen for linear epitopes.", default = 25, type = int, choices = range(6,26),required=True)

    args = parser.parse_args()
    fasta = args.fasta
    organism = args.organism
    window = args.window

    main(fasta,organism,window,output_dir)