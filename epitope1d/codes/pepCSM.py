#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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

##pep_CSM performs a scanning matrix into peptide sequence according to cutoff_limit size (window length) and cutoff_step (incremental step).
##This functions is cumulative, meaning that for each sequence, each step count is added to the next step of the same group.

# ____________________________________________________________________________________________________________________
# Residue classification 4 (5 classes)

# Apolar:         G,A,V,L,I,P,M (7)
# Aromatic:       F,Q,Y         (3)
# Polar neutral:  S,C,T,N,Q     (5)
# Basic:          R,H,K         (3)
# Acidic:         E,D           (2)

Group4_values = {}
Group4_values['Apolar'] = 1
Group4_values['Aromatic'] = 1
Group4_values['PolarNeutral'] = 1
Group4_values['Acidic'] = 1
Group4_values['Basic'] = 1

Group4 = {}
Group4['G'] = "Apolar"
Group4['A'] = "Apolar"
Group4['V'] = "Apolar"
Group4['L'] = "Apolar"
Group4['I'] = "Apolar"
Group4['P'] = "Apolar"
Group4['M'] = "Apolar"

Group4['F'] = "Aromatic"
Group4['W'] = "Aromatic"
Group4['Y'] = "Aromatic"

Group4['S'] = "PolarNeutral"
Group4['C'] = "PolarNeutral"
Group4['T'] = "PolarNeutral"
Group4['N'] = "PolarNeutral"
Group4['Q'] = "PolarNeutral"

Group4['R'] = "Basic"
Group4['H'] = "Basic"
Group4['K'] = "Basic"

Group4['E'] = "Acidic"
Group4['D'] = "Acidic"

#pep_CSM.py <input_data> <cuttof_limit> <cuttof_step> <sign_type {4}> \n
   # \n\tWhere:\
   # \n\t<input_data> Input file - should have 2 columns (csv file): \
   # \n\t<Fasta_header>,<Peptide>\n\
   # Where:\n\
   # \n\t\t<Fasta_header> - name of the sequence\
   # \n\t\t<Peptide> - peptide sequence\
   # \n\t<cutoff_limit> and <cuttof_step> are the (integer) parameters for the cutoff scanning;\
   # \n\t<sign_type> {4} \
   #         \n\t\t# 4:  5 residue classes {Apolar, Acidic, Aromatic, Basic, PolarNeutral,}\n\
#____________________________________________________________________________________________________________________
import pandas as pd
import warnings
import itertools as it
warnings.filterwarnings("ignore")

def pepCSM(input_data, dmax, dmin, cutoff_step, sign_type):

    in_file = input_data  # Input csv file
    data = pd.DataFrame(in_file)
    cutoff_limit = int(dmax) #maximum distance
    dmin = int(dmin)  # minimum distance to evaluate
    cutoff_step = int(cutoff_step)   # Cutoff Step distance
    sign_type = int(sign_type)

    # Signature type
    # Setting parameters for automatic file generation
    if (cutoff_step < 1):
        print ("\nError: Cutoff Step too small.\nChoose another please.\n\n")
        exit()
    if (cutoff_step >= cutoff_limit):
        print ("\nError: Cutoff Step should be smaller than cutoff_limit (window size).\nChoose another please.\n\n")
        exit()
    if(sign_type > 4 or sign_type < 1):
        print ("\nError: Please select an appropriate signature type {1,2,3,4}.\n\n")
        exit()

    file = []
    # Generating header for the output file
    header = []
    # Signature header
    for x in range(dmin, (cutoff_limit+1), cutoff_step):
        keys = {}
        if(sign_type == 4):
            keys = sorted(Group4_values)

    # Getting linear combination 2 to 2 with keys values
        combination_keys = [list(it.combinations_with_replacement(keys,2))]
        for a in range(0,len(combination_keys[0])):
            header.append(":".join(combination_keys[0][a]) + "-" + str(x))

    file.append(header)
    #--------------------------------------Classification Functions----------------------------------------------------------------
    for index, row in data.iterrows():
        features_list = []
        seqs = str(row[1])
        size_seqs = len(seqs)

        num_AcAc = 0
        num_AcA = 0
        num_AcAr = 0
        num_AcB = 0
        num_AcPn = 0
        num_AA = 0
        num_AAr = 0
        num_AB = 0
        num_APn = 0
        num_ArAr = 0
        num_ArB = 0
        num_ArPn = 0
        num_BB = 0
        num_BPn = 0
        num_PnPn = 0
        temp_feat = []

        for dist in range(dmin, (cutoff_limit+1), cutoff_step):
            for i in range(0, size_seqs):
                if (i + dist) < size_seqs:
                    slice_pep_1 = seqs[i]
                    slice_pep_2 = seqs[(i + dist)]

                    if sign_type == 4:
                        if Group4.get(slice_pep_1) == 'Acidic' and Group4.get(slice_pep_2) == 'Acidic':
                            num_AcAc +=1 #Acidic-Acidic
                        elif (Group4.get(slice_pep_1) == 'Acidic' and Group4.get(slice_pep_2) == 'Apolar') or (Group4.get(slice_pep_1) == 'Apolar' and Group4.get(slice_pep_2) == 'Acidic'):
                            num_AcA += 1 #Acidic-Apolar
                        elif (Group4.get(slice_pep_1) == 'Acidic' and Group4.get(slice_pep_2) == 'Aromatic') or (Group4.get(slice_pep_1) == 'Aromatic' and Group4.get(slice_pep_2) == 'Acidic'):
                            num_AcAr += 1 #Acidic-Aromatic
                        elif (Group4.get(slice_pep_1) == 'Acidic' and Group4.get(slice_pep_2) == 'Basic') or (Group4.get(slice_pep_1) == 'Basic' and Group4.get(slice_pep_2) == 'Acidic'):
                            num_AcB += 1 #Acidic-Basic
                        elif (Group4.get(slice_pep_1) == 'Acidic' and Group4.get(slice_pep_2) == 'PolarNeutral') or (Group4.get(slice_pep_1) == 'PolarNeutral' and Group4.get(slice_pep_2) == 'Acidic'):
                            num_AcPn += 1 #Acidic-Polar Neutral
                        elif Group4.get(slice_pep_1) == 'Apolar' and Group4.get(slice_pep_2) == 'Apolar':
                            num_AA += 1 #Apolar-Apolar
                        elif Group4.get(slice_pep_1) == 'Apolar' and Group4.get(slice_pep_2) == 'Aromatic' or (Group4.get(slice_pep_1) == 'Aromatic' and Group4.get(slice_pep_2) == 'Apolar'):
                            num_AAr += 1 #Apolar-Aromatic
                        elif (Group4.get(slice_pep_1) == 'Apolar' and Group4.get(slice_pep_2) == 'Basic') or (Group4.get(slice_pep_1) == 'Basic' and Group4.get(slice_pep_2) == 'Apolar'):
                            num_AB += 1 #Apolar-Basic
                        elif (Group4.get(slice_pep_1) == 'Apolar' and Group4.get(slice_pep_2) == 'PolarNeutral') or (Group4.get(slice_pep_1) == 'PolarNeutral' and Group4.get(slice_pep_2) == 'Apolar'):
                            num_APn += 1 #Apolar-Polar Neutral
                        elif Group4.get(slice_pep_1) == 'Aromatic' and Group4.get(slice_pep_2) == 'Aromatic':
                            num_ArAr += 1 #Aromatic-Aromatic
                        elif (Group4.get(slice_pep_1) == 'Aromatic' and Group4.get(slice_pep_2) == 'Basic') or (Group4.get(slice_pep_1) == 'Basic' and Group4.get(slice_pep_2) == 'Aromatic'):
                            num_ArB += 1 #Aromatic-Basic
                        elif (Group4.get(slice_pep_1) == 'Aromatic' and Group4.get(slice_pep_2) == 'PolarNeutral') or (Group4.get(slice_pep_1) == 'PolarNeutral' and Group4.get(slice_pep_2) == 'Aromatic'):
                            num_ArPn += 1 #Aromatic - Polar Neutral
                        elif Group4.get(slice_pep_1) == 'Basic' and Group4.get(slice_pep_2) == 'Basic':
                            num_BB += 1 #Basic-Basic
                        elif (Group4.get(slice_pep_1) == 'Basic' and Group4.get(slice_pep_2) == 'PolarNeutral') or (Group4.get(slice_pep_1) == 'PolarNeutral' and Group4.get(slice_pep_2) == 'Basic'):
                            num_BPn += 1 #Basic-Polar Neutral
                        elif Group4.get(slice_pep_1) == 'PolarNeutral' and Group4.get(slice_pep_2) == 'PolarNeutral':
                            num_PnPn += 1 #Polar Neutral-Polar Neutral
                        else:
                            print("## ERROR: Peptide Error! {} and {} \n".format(slice_pep_1, slice_pep_2))
                        temp_feat = [num_AcAc,num_AcA,num_AcAr,num_AcB,num_AcPn,num_AA,num_AAr,num_AB,num_APn,num_ArAr,num_ArB,num_ArPn,num_BB,num_BPn,num_PnPn]
                else:
                    continue
            features_list.append(temp_feat)

        features_str = str(features_list).replace('[','').replace(']','').strip(" ").replace(" ","")  ##removing characters
        features = list(map(int,features_str.split(','))) #transfom string of numbers into list of integers
        file.append(features)

    return file






