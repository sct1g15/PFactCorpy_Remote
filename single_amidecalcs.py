import argparse
import os
import pandas as pd
import numpy
import itertools
from kint import *
from SM8 import run_SM8
from Bio import SeqIO
from more_itertools import unique_everseen




def find_plus1(Protein, State, Exposure):
    '''
    :param Protein: Protein name in excel file, e.g. HLA A0201
    :param State: Experimental conditions label
    :param Exposure: Time exposure
    :return:
    '''
    ex_indata = pd.read_excel(args.excelin)

    ex_data = ex_indata[ex_indata['Protein'] == Protein]
    ex_data = ex_data[ex_data['State'] == State]
    ex_data = ex_data[ex_data['Exposure'] == Exposure]

    duplicatestart = ex_data[ex_data.duplicated(['Start'])]
    duplicatestart = unique_everseen(duplicatestart['Start'])
    duplicatestart = list(duplicatestart)

    sapos = pd.DataFrame()
    single_amides = pd.DataFrame()
    for i in duplicatestart:
        duplicateset = ex_data[ex_data['Start'] == i]
        for j in range(0, len(duplicateset)-1):
            #print(duplicateset.iloc[j,2]+1)
            #print(duplicateset.iloc[j+1,2])
            if duplicateset.iloc[j,2]+1 == duplicateset.iloc[j+1,2]:
                duplicateset2 = ex_data[ex_data['Start'] == i]
                for q in range(0, len(duplicateset2)-1):
                    if duplicateset2.iloc[q, 2]+1 == duplicateset2.iloc[q+1, 2]:
                        sapos = pd.DataFrame({
                            'Resi' : duplicateset2.iloc[q+1, 2],
                            'Resn' : str(duplicateset2.iloc[q+1, 3])[-1],
                            'Uptake' : float(duplicateset2.iloc[q+1, 12] - duplicateset2.iloc[q, 12]),
                            'Uptake SD' : duplicateset2.iloc[q+1, 13],
                            'Type' : '=-',
                            'State' : str(State),
                            'Exposure' : str(Exposure),
                            'Retention Time' : duplicateset2.iloc[q+1, 14],
                            'Retention Time SD': duplicateset2.iloc[q + 1, 15]
                        },
                        index = [0])
                    single_amides = single_amides.append(sapos)
    single_amides1 = single_amides.drop_duplicates()
    return single_amides1

def find_neg1(Protein, State, Exposure):
    '''

    :param Protein: Protein name in excel file, e.g. HLA A0201
    :param State: Experimental conditions label
    :param Exposure: Time exposure
    :return:
    '''
    ex_indata = pd.read_excel(args.excelin)

    ex_data = ex_indata[ex_indata['Protein'] == Protein]
    ex_data = ex_data[ex_data['State'] == State]
    ex_data = ex_data[ex_data['Exposure'] == Exposure]

    duplicateend = ex_data[ex_data.duplicated(['End'])]
    duplicateend = unique_everseen(duplicateend['End'])
    duplicateend = list(duplicateend)

    single_amides = pd.DataFrame()
    sapos=pd.DataFrame()
    for i in duplicateend:
        duplicateset = ex_data[ex_data['End'] == i]
        for j in range(0, len(duplicateset)-1):
            #print(duplicateset.iloc[j,2]+1)
            #print(duplicateset.iloc[j+1,2])
            if duplicateset.iloc[j,1]+1 == duplicateset.iloc[j+1,1]:
                duplicateset2 = ex_data[ex_data['End'] == i]
                for q in range(0, len(duplicateset2)-1):
                    if duplicateset2.iloc[q, 1]+1 == duplicateset2.iloc[q+1, 1]:
                        sapos = pd.DataFrame({
                            'Resi' : duplicateset2.iloc[q, 1]+1,
                            'Resn' : str(duplicateset2.iloc[q, 3])[1],
                            'Uptake' : float(duplicateset2.iloc[q, 12] - duplicateset2.iloc[q+1, 12]),
                            'Uptake SD': duplicateset2.iloc[q, 13],
                            'Type' : '-=',
                            'State': str(State),
                            'Exposure': str(Exposure),
                            'Retention Time': duplicateset2.iloc[q + 1, 14],
                            'Retention Time SD': duplicateset2.iloc[q + 1, 15]
                        },
                        index = [0])
                    single_amides = single_amides.append(sapos)
    single_amides2 = single_amides.drop_duplicates()
    return single_amides2

def collate_singleamide_data():
    Exposures = ex_indata.Exposure.dropna().unique()
    States = ex_indata.State.dropna().unique()
    Proteins = ex_indata.Protein.dropna().unique()

    output = pd.DataFrame()
    for i in Exposures:
        for j in States:
            for p in Proteins:
                output = output.append(find_plus1(p, j, i))
                output = output.append(find_neg1(p, j, i))
    print(output)
    output.to_csv(os.path.splitext(args.excelin)[0] + '_singleamide.csv')

# collate_singleamide_data()