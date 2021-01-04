# -*- coding: utf-8 -*-

import argparse
import numpy as np
from kint import *
from PepFrag_predict import PepFrag_find, PepFrag_predict, PepFrag_uptake
from SM8 import *
from metric_measure import *
#from more_itertools import unique_everseen
import matplotlib.pyplot as plt
import scipy.stats

time = 0.5

parser = argparse.ArgumentParser()
parser.add_argument('-x', '--excelin', help='Excel Input File')
parser.add_argument('-t', '--temp', help='Experimental Temp')
parser.add_argument('-p', '--pD', help='Experimental pD')
parser.add_argument('-f', '--trajin', help='Trajectory xtc file')
parser.add_argument('-s', '--topolin', help='Topology file')
parser.add_argument('-st', '--state', help='State to analyse (from excel file header)')
args = parser.parse_args()

Seq = 'GSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDGETRKVKAHSQTHRVDLGTLRGYYNQSEAGSHTVQRMYGCDVGSDWR' \
      'FLRGYHQYAYDGKDYIALKEDLRSWTAADMAAQTTKHKWEAAHVAEQLRAYLEGTCVEWLRRYLENGKETLQRTDAPKTHMTHHAVSDHEATLRCWALSFYPAEITLT' \
      'WQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGQEQRYTCHVQHEGLPKPLTLRWE'

#Get Intrinsic kint rates for sequence
kint = []
for counter, aminoacid in enumerate(Seq[0:len(Seq)]):
    kint.append(calc_kint_Alt(Seq[counter], Seq[counter-1], int(args.temp), int(args.pD), counter, len(Seq)))
print(kint)

# Get Single Amide P Values
# sasaavg, sasamax = get_sasa(args.trajin, args.topolin, step=15000) #approx 40 minute operation for 150,000 frames, 1 core, 40 threads
# distavg = get_amideH_polardistances_Parallel(args.trajin, args.topolin, step=200) #approx 25 hour operarion, same stats

sasamax = np.loadtxt('amidesasamax.out')
sasamax = sasamax*100
distavg = np.loadtxt('polardistance_mean.out')


P_values = []
for i, j in zip(sasamax, distavg):
    x = run_SM8(i, j)
    P_values.append(x)
#print(P_values)
np.savetxt('P_values_predictions.out', P_values)

#Deuterium uptake calculation per peptide
ex_indata = pd.read_excel(args.excelin)
Protein_name = ex_indata['Protein'].unique()

Pepfrags = PepFrag_find(Protein_name[0], args.state, time, args.excelin)

Du_values = []
for index, row in Pepfrags.iterrows():
    Pval = []
    Kintval = []
    if row['End'] > len(kint):
        continue
    for j in range(int(row['Start']), int(row['End'])):
        Pval.append(P_values[j])
        Kintval.append(kint[j])
    x = PepFrag_predict(Pval, Kintval, time)
    Du_values.append(x)

Pepfrags_uptake = PepFrag_uptake(Protein_name[0], args.state, time, args.excelin)
Pepfrags_uptake = Pepfrags_uptake.to_numpy()

#plot R2 graphs###
plt.scatter(Du_values, Pepfrags_uptake[0:len(Du_values)])
plt.ylabel('Actual Du values (D/amide/min)')
plt.xlabel('Predicted Du values (D/amide/min)')
print(scipy.stats.pearsonr(Du_values, Pepfrags_uptake[0:len(Du_values)]))
print(scipy.stats.kendalltau(Du_values, Pepfrags_uptake[0:len(Du_values)]))
plt.show()