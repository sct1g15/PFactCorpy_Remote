#!/scratch/sct1g15/mydocuments/conda/pymol/bin/python

#SBATCH --ntasks-per-node=40
#SBATCH --nodes=1
#SBATCH --time=02:00:00
#SBATCH --partition=scavenger

import os
import sys

os.chdir('/scratch/sct1g15/mydocuments/year1/PFactcorpy')
sys.path.append('/scratch/sct1g15/mydocuments/year1/PFactcorpy')

from math import tanh, exp, log
import numpy as np
import argparse
from scipy.optimize import curve_fit
import pandas as pd
from kint import calc_kint_Alt
from PepFrag_predict import PepFrag_uptake, PepFrag_find, PepFrag_predict
from multiprocessing import Pool, get_context
import matplotlib.pyplot as plt

time = 0.5

parser = argparse.ArgumentParser()
parser.add_argument('-x', '--excelin', help='Excel Input File', default='20180705 A0201 4 state.xlsx')
parser.add_argument('-t', '--temp', help='Experimental Temp', default=298)
parser.add_argument('-p', '--pD', help='Experimental pD', default=7)
parser.add_argument('-st', '--state', help='State to analyse (from excel file header)', default='UV EXP')
parser.add_argument('-np', '--procs', help='Number of processors', default=40)
args = parser.parse_args()

# Input SASA and DIST Histogram data
SASA = np.load('dist_sasa_20bins/SASA_hist.npy', allow_pickle=True)
print('SASA loaded.')
DIST = np.load('dist_sasa_20bins/DIST_hist.npy', allow_pickle=True)
print('DIST loaded.')

# Get kint data
Seq = 'MGSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDGETRKVKAHSQTHRVDLGTLRGYYNQSEAGSHTVQRMYGCDVGSDWR' \
      'FLRGYHQYAYDGKDYIALKEDLRSWTAADMAAQTTKHKWEAAHVAEQLRAYLEGTCVEWLRRYLENGKETLQRTDAPKTHMTHHAVSDHEATLRCWALSFYPAEITLT' \
      'WQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGQEQRYTCHVQHEGLPKPLTLRWE'

kint = []
for counter, aminoacid in enumerate(Seq[0:len(Seq)]):
    kint.append(calc_kint_Alt(Seq[counter], Seq[counter - 1], int(args.temp), int(args.pD), counter, len(Seq)))

# Input Peptide DU data
ex_indata = pd.read_excel(args.excelin)
Protein_name = ex_indata['Protein'].unique()
Pepfrags = PepFrag_find(Protein_name[0], args.state, time, args.excelin)
Pepfrags_uptake = PepFrag_uptake(Protein_name[0], args.state, time, args.excelin)
Pepfrags_uptake = Pepfrags_uptake.to_numpy()

count = 0
Pval = []
for index, row in Pepfrags.iterrows():
    print(row)
    kintval = []
    if row['End'] > len(kint):
        continue
    for j in range(int(row['Start']), int(row['End'])):
        kintval.append(kint[j])
    kintval_avg = np.array(kintval).mean()
    #P = ((-1 * kintval / log(Pepfrags_uptake[count])) * time) / (int(row['End']) - int(row['Start']))
    print(index,':kintval_avg:',kintval_avg)
    print(index,':kintval:',kintval)
    print(index,':DU Uptake:',Pepfrags_uptake[count])
    print(index,':N-count:',int((int(row['End']) - int(row['Start']))))
    print(index,':N-count_noPro:', int(row['MaxUptake']))
    #Calculate Average Protection Factor from DU and Kint
    v1 = (Pepfrags_uptake[count])
    v2 = -1*(v1-1)
    v3 = log(v2)/time
    print(v1)
    print(v2)
    print(v3)
    P_factor = -1*(kintval_avg)/v3
    lnPF = log(P_factor)
    # P = (-1 * (kintval/(int(row['End']) - int(row['Start'])))) / log((Pepfrags_uptake[count] / (int(row['End']) - int(row['Start'])))) * time
    print(index,':Protection factor:',P_factor)
    print(index,':lnPF:',lnPF)
    # if P > 10:
    #     P = 10
    Pval.append(lnPF)
    count += 1

np.savetxt('Pval_exp_fittingset.txt', Pval)

start = []
end = []
for index, row in Pepfrags.iterrows():
    start.append(int(row['Start']))
    end.append(int(row['End']))
start = np.array(start)
end = np.array(end)

# Generate comm split
procs = int(args.procs)
ave, res = divmod(count, procs)
counts = [ave + 1 if p < res else ave for p in range(procs)]
starts = [sum(counts[:p]) for p in range(procs)]
ends = [sum(counts[:p + 1]) for p in range(procs)]
data = [(starts[p], ends[p]) for p in range(procs)]

p0 = 100, 1, 2, 2, 2, 5, 4.5, 8
p0 = 0.0013, 2.64, 36.5, 1.27
print(data)

def HM8_model_multiprocessing(X, alpha, beta, gamma, delta, epsilon, xp, yp, zp):
    # print('alpha:', alpha, 'beta:', beta, 'gamma:', gamma, 'delta:', delta, 'epsilon:', epsilon, 'x,y,z:', xp, yp,
    #       zp)
    # print(X)
    Start, End = X
    pool = Pool(processes=args.procs)
    result = pool.starmap(HM8_model_method,
                               [(X, alpha, beta, gamma, delta, epsilon, xp, yp, zp) for X in zip(Start, End)])
    pool.close()
    pool.join()
    output = [y for x in result for y in x]
    return output

def HM8_model_multiprocessing_M8(X, betas, gammas, betap, gammap):
    # print('alpha:', alpha, 'beta:', beta, 'gamma:', gamma, 'delta:', delta, 'epsilon:', epsilon, 'x,y,z:', xp, yp,
    #       zp)
    # print(X)
    Start, End = X
    pool = Pool(processes=args.procs)
    result = pool.starmap(M8_model_method,
                               [(X, betas, gammas, betap, gammap) for X in zip(Start, End)])
    pool.close()
    pool.join()
    output = [y for x in result for y in x]
    return output


def M8_model_method(X, betas, gammas, betap, gammap):
    Start, End = X
    Pval_pred = []
    Pval_comb = []
    for j in range(int(Start), int(End)):
        # if kint[j] == -1:
        #     P = 0
        P = 0
        count = 0
        for q in range(0, len(SASA[j][0])):
            for p in range(0, len(DIST[j][0])):
                count += 1
                P += (((SASA[j][0][q] * DIST[j][0][p]) / (sum(SASA[j][0]) * sum(DIST[j][0]))) * ((betas*(SASA[j][1][q])**gammas)+(betap*(DIST[j][1][p])**gammap)))
        #print('Individual Amide ',j,' P-value:',P)
        Pval_comb.append(P)
    print(np.array(Pval_comb).mean())
    Pval_pred.append(np.array(Pval_comb).mean())
    #frag_DU.append(np.array(sa_DU).mean())
    return Pval_pred


def HM8_model_method(X, alpha, beta, gamma, delta, epsilon, xp, yp, zp):
    Start, End = X
    Pval_pred = []
    Pval_comb = []
    for j in range(int(Start), int(End)):
        # if kint[j] == -1:
        #     P = 0
        P = 0
        count = 0
        for q in range(0, len(SASA[j][0])):
            for p in range(0, len(DIST[j][0])):
                count += 1
                P += (((SASA[j][0][q] * DIST[j][0][p]) / (sum(SASA[j][0]) * sum(DIST[j][0]))) * (((epsilon * (
                    tanh((-alpha * SASA[j][1][q]) + xp)) * tanh((beta * DIST[j][1][p]) - yp)) - gamma * (tanh(
                    (-alpha * SASA[j][1][q]) + xp)) + delta * (tanh((beta * DIST[j][1][p]) - yp))) + zp))
        #print('Individual Amide ',j,' P-value:',P)
        Pval_comb.append(P)
        # if P >= 0:
        #     Pval_comb.append(P)
        #     sa_DU.append(1-exp((-kint[j]/exp(P))*time))
        # if P < 0:
        #     sa_DU.append(9999999)
    print(np.array(Pval_comb).mean())
    Pval_pred.append(np.array(Pval_comb).mean())
    #frag_DU.append(np.array(sa_DU).mean())
    return Pval_pred

start = start[0:count]
end = end[0:count]

if __name__ == "__main__":
    start = start[0:count]
    end = end[0:count]

    counter = 0

    popt, pcov = curve_fit(HM8_model_multiprocessing_M8, xdata=(start, end), ydata=Pval, p0=p0)

    #alpha, beta, gamma, delta, epsilon, xp, yp, zp = popt
    #alpha, beta, gamma, delta, epsilon, xp, yp, zp  = 100, 1, 2, 2, 2, 5, 4.5, 8
    betas, gammas, betap, gammap = popt
    print(popt)
    print(pcov)

    model='M8'
    if model == 'HM8':
        Pval_pred = []
        for X in zip(start, end):
            Start, End = X
            Pval_comb = []

            for j in range(int(Start), int(End)):
                # if kint[j] == -1:
                #     P = 0
                P = 0
                count = 0
                for q in range(0, len(SASA[j][0])):
                    for p in range(0, len(DIST[j][0])):
                        count += 1
                        P += (((SASA[j][0][q] * DIST[j][0][p]) / (sum(SASA[j][0]) * sum(DIST[j][0]))) * (((epsilon * (
                            tanh((-alpha * SASA[j][1][q]) + xp)) * tanh((beta * DIST[j][1][p]) - yp)) - gamma * (tanh(
                            (-alpha * SASA[j][1][q]) + xp)) + delta * (tanh((beta * DIST[j][1][p]) - yp))) + zp))
                if P >= 0:
                    Pval_comb.append(P)
                    #sa_DU.append(1-exp((-kint[j]/exp(P))*time))
                #if P < 0:
                    #sa_DU.append(9999999)
            Pval_pred.append(np.array(Pval_comb).mean())
            #frag_DU.append(np.array(sa_DU).mean())
    if model == 'M8':
        Pval_pred = []
        for X in zip(start, end):
            Start, End = X
            Pval_comb = []

            for j in range(int(Start), int(End)):
                # if kint[j] == -1:
                #     P = 0
                P = 0
                count = 0
                for q in range(0, len(SASA[j][0])):
                    for p in range(0, len(DIST[j][0])):
                        count += 1
                        P += (((SASA[j][0][q] * DIST[j][0][p]) / (sum(SASA[j][0]) * sum(DIST[j][0]))) * (
                                    (betas * (SASA[j][1][q]) ** gammas) + (betap * (DIST[j][1][p]) ** gammap)))
                        # P += (((SASA[j][0][q] * DIST[j][0][p]) / (sum(SASA[j][0]) * sum(DIST[j][0]))) * (((epsilon * (
                        #     tanh((-alpha * SASA[j][1][q]) + xp)) * tanh((beta * DIST[j][1][p]) - yp)) - gamma * (tanh(
                        #     (-alpha * SASA[j][1][q]) + xp)) + delta * (tanh((beta * DIST[j][1][p]) - yp))) + zp))
                if P >= 0:
                    Pval_comb.append(P)
                    # sa_DU.append(1-exp((-kint[j]/exp(P))*time))
                # if P < 0:
                # sa_DU.append(9999999)
            Pval_pred.append(np.array(Pval_comb).mean())
            # frag_DU.append(np.array(sa_DU).mean())

    # Pval_exp = []
    # count = 0
    # for index, row in Pepfrags.iterrows():
    #     kintval = []
    #     if row['End'] > len(kint):
    #         continue
    #     for j in range(int(row['Start']), int(row['End'])):
    #         kintval.append(kint[j])
    #     kintval = np.array(kintval).mean()
    #     DUval = Pepfrags_uptake[count]
    #     # P = (-1*kintval/log(DUval/(int(row['End']) - int(row['Start']))))*time
    #     v1 = (Pepfrags_uptake[count] / int(row['End']) - int(row['Start']))
    #     v2 = -1 * (v1 - 1)
    #     v3 = log(v2) / time
    #     P = kintval / v3
    #     P = log(P)
    #
    #     # P = (-1 * (kintval / (int(row['End']) - int(row['Start'])))) / log(
    #     #     (DUval / (int(row['End']) - int(row['Start'])))) * time
    #     Pval_exp.append(P)
    #     count += 1

    from scipy.stats import pearsonr

    np.savetxt('Pval_exp.txt', Pval)
    np.savetxt('Pval_pred.txt', Pval_pred)
    plt.scatter(Pval_pred, Pval)
    corr, _ = pearsonr(Pval, Pval_pred)
    plt.xlabel('Protection Factor Prediction')
    plt.ylabel('Protection Factor Experimental')
    plt.annotate("Pearsons R = {:.3f}".format(corr), (0, 1))
    plt.savefig('Protection Factors.png')

    # Pval_exp = []
    # sasaval_avg = []
    # distval_avg = []
    # kintval_avg = []
    # count = 0
    #
    # for index, row in Pepfrags.iterrows():
    #     sasaval_peptide=[]
    #     distval_peptide=[]
    #     if row['End'] > len(kint):
    #         continue
    #     for j in range(int(row['Start']), int(row['End'])):
    #         sasaval = []
    #         distval = []
    #         for i in range(0, len(SASA[j][0])): #across histogram bins
    #             sasaval.append(float(SASA[j][1][i])*float(SASA[j][0][i]))   #Slight innacuracy, each bin is averaged by its upper limit
    #             print(sasaval, SASA[j][0])
    #         sasaval_norm = sum(sasaval)/sum(SASA[j][0])
    #         for i in range(0, len(DIST[j][0])):
    #             distval.append(DIST[j][1][i]*DIST[j][0][i])
    #         distval_norm = sum(distval) / sum(DIST[j][0])
    #         sasaval_peptide.append(sasaval_norm)
    #         distval_peptide.append(distval_norm)
    #     sasaval_avg.append(np.array(sasaval_peptide).mean())
    #     distval_avg.append(np.array(distval_peptide).mean())
    #
    # np.savetxt('sasaval_avg.txt', sasaval_avg)
    # np.savetxt('distval_avg.txt', distval_avg)

