import argparse
import os
import pandas as pd
import numpy as np
import itertools
from math import exp

# parser = argparse.ArgumentParser()
# parser.add_argument('-x', '--excelin', help='Excel Input File')
# args = parser.parse_args()

# def PepFrag_find_solo(Protein, State, Exposure):
#     ex_indata = pd.read_excel(args.excelin)
#
#     ex_data = ex_indata[ex_indata['Protein'] == Protein]
#     ex_data = ex_data[ex_data['State'] == State]
#     ex_data = ex_data[ex_data['Exposure'] == Exposure]
#
#     #standardise numbering, due to M starting residue in Mass Spec, not in MD
#     ex_data = ex_data[ex_data["Start"]!=1.0]
#     ex_data["Start"] -= 1
#     ex_data["End"] -= 1
#
#     print(ex_data)

# PepFrag_find_solo('1A02_HUMAN', 'UV exposed', 0.5)

def PepFrag_find(Protein, State, Exposure, excelin):
    print(Protein, State, Exposure)
    ex_indata = pd.read_excel(excelin)

    ex_data = ex_indata[ex_indata['Protein'] == Protein]
    ex_data = ex_data[ex_data['State'] == State]
    ex_data = ex_data[ex_data['Exposure'] == Exposure]

    #standardise numbering, due to M starting residue in Mass Spec, not in MD
    ex_data = ex_data[ex_data["Start"]!=1.0]
    ex_data["Start"] -= 1
    ex_data["End"] -= 1
    return(ex_data)

def PepFrag_uptake(Protein, State, Exposure, excelin):
    print(Protein, State, Exposure)
    ex_indata = pd.read_excel(excelin)

    ex_data = ex_indata[ex_indata['Protein'] == Protein]
    ex_data = ex_data[ex_data['State'] == State]
    ex_data = ex_data[ex_data['Exposure'] == Exposure]

    #standardise numbering, due to M starting residue in Mass Spec, not in MD
    ex_data = ex_data[ex_data["Start"]!=1.0]
    ex_data["Start"] -= 1
    ex_data["End"] -= 1
    return(np.divide(ex_data["Uptake"], ex_data["MaxUptake"]))

def PepFrag_predict(P, kint, time):
    pep_len = len(P)
    P = P[1:]
    kint = kint[1:]

    Frag_DU_boltz = 0
    if len(P) != len(kint):
        print('ERROR: Different number of P-values to Kint values passed to function')
    for i, k in zip(P, kint):
        Frag_DU_boltz += 1-exp((-k/i)*time)
    Frag_DU = (1/pep_len)*Frag_DU_boltz
    return(Frag_DU)

