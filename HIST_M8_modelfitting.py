from math import tanh, exp, log
import numpy as np
import argparse
import scipy
from scipy.optimize import curve_fit
import pandas as pd
from kint import calc_kint_Alt
from PepFrag_predict import PepFrag_uptake, PepFrag_find, PepFrag_predict

time = 0.5

parser = argparse.ArgumentParser()
parser.add_argument('-x', '--excelin', help='Excel Input File')
parser.add_argument('-t', '--temp', help='Experimental Temp')
parser.add_argument('-p', '--pD', help='Experimental pD')
parser.add_argument('-st', '--state', help='State to analyse (from excel file header)')
args = parser.parse_args()

#Input SASA and DIST Histogram data
SASA = np.load('SASA_hist.npy', allow_pickle=True)
DIST = np.load('DIST_hist.npy', allow_pickle=True)

#Get kint data
Seq = 'MGSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDGETRKVKAHSQTHRVDLGTLRGYYNQSEAGSHTVQRMYGCDVGSDWR' \
      'FLRGYHQYAYDGKDYIALKEDLRSWTAADMAAQTTKHKWEAAHVAEQLRAYLEGTCVEWLRRYLENGKETLQRTDAPKTHMTHHAVSDHEATLRCWALSFYPAEITLT' \
      'WQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGQEQRYTCHVQHEGLPKPLTLRWE'

kint = []
for counter, aminoacid in enumerate(Seq[0:len(Seq)]):
    kint.append(calc_kint_Alt(Seq[counter], Seq[counter-1], int(args.temp), int(args.pD), counter, len(Seq)))


#Input Peptide DU data
ex_indata = pd.read_excel(args.excelin)
Protein_name = ex_indata['Protein'].unique()
Pepfrags = PepFrag_find(Protein_name[0], args.state, time, args.excelin)
Pepfrags_uptake = PepFrag_uptake(Protein_name[0], args.state, time, args.excelin)
Pepfrags_uptake = Pepfrags_uptake.to_numpy()

count=0
Pval = []
for index, row in Pepfrags.iterrows():
    kintval = []
    if row['End'] > len(kint):
        continue
    for j in range(int(row['Start']), int(row['End'])):
        kintval.append(kint[j])
    kintval = np.array(kintval).mean()
    P = ((-1 * kintval / log(Pepfrags_uptake[count])) * time) / (int(row['End']) - int(row['Start']))
    # if P > 7:
    #     P = 7
    Pval.append(P)
    count =+ 1

start = []
end = []
for index, row in Pepfrags.iterrows():
    start.append(int(row['Start']))
    end.append(int(row['End']))
start = np.array(start)
end = np.array(end)

p0 = 100,1,2,2,2,5,4.5,8

def HM8_model(X, alpha, beta, gamma, delta, epsilon, xp, yp, zp):
    #print('alpha:',alpha,'beta:',beta,'gamma:',gamma,'delta:',delta,'epsilon:', epsilon,'x,y,z:', xp, yp, zp)
    Start, End = X
    frag_DU = []
    Pval_pred = []
    for i in range(10, len(Start)):
        sa_DU = []
        Pval_comb = []
        if int(End[i]) > len(kint):
            continue
        for j in range(int(Start[i]), int(End[i])):
            # if kint[j] == -1:
            #     P = 0
            P=0
            count = 0
            for q in range(0, len(SASA[j][0])):
                for p in range(0, len(DIST[j][0])):
                    count += 1
                    P += (((SASA[j][0][q]*DIST[j][0][p])/(sum(SASA[j][0])*sum(DIST[j][0])))*(((epsilon * (tanh((-alpha * SASA[j][1][q])+xp))*tanh((beta * DIST[j][1][p])-yp)) - gamma*(tanh((-alpha * SASA[j][1][q])+xp)) + delta*(tanh((beta * DIST[j][1][p])-yp)))+zp))
                    # print('Number SASA:', SASA[j][0][q])
                    # print('Number DIST:', DIST[j][0][p])
                    # print('Relative Contribution:', (SASA[j][0][q]+DIST[j][0][p])/(sum(SASA[j][0])*sum(DIST[j][0])))
                    # print('P-value of selection:', (((epsilon * (tanh((-alpha * SASA[j][1][q])+xp))*tanh((beta * DIST[j][1][p])-yp)) - gamma*(tanh((-alpha * SASA[j][1][q])+xp)) + delta*(tanh((beta * DIST[j][1][p])-yp)))+zp))
                    # print('Final P value:', (((SASA[j][0][q]+DIST[j][0][p])/(sum(SASA[j][0])*sum(DIST[j][0])))*(((epsilon * (tanh((-alpha * SASA[j][1][q])+xp))*tanh((beta * DIST[j][1][p])-yp)) - gamma*(tanh((-alpha * SASA[j][1][q])+xp)) + delta*(tanh((beta * DIST[j][1][p])-yp)))+zp)))
                    # print('Cumulative P value:', P)
                    # print('count:', count)
                    # print(((epsilon * (tanh((-alpha * 0.045)+xp))*tanh((beta * 4.5)-yp)) - gamma*(tanh((-alpha * 0.045)+xp)) + delta*(tanh((beta * 4.5)-yp)))+zp)
            if P >= 0:
                Pval_comb.append(P)
                #print(kint[j], P)
                sa_DU.append(1-exp((-kint[j]/exp(P))*time))
                # print(Pval)
                # print(kint[j])
                # print(sa_DU)
            if P < 0:
                sa_DU.append(9999999)
                print('X')
        Pval_pred.append(np.array(Pval_comb).mean())
        frag_DU.append(np.array(sa_DU).mean())
        print('PVAL-Prediction:', Pval_pred[0+(i-10)])
        print('FRAG_DU-Prediction:', frag_DU[0])
        print('PVAL:', Pval[0+i])
    print('PVAL-Prediction:', Pval_pred[0])
    print('PVAL:',Pval[10])
    #print(len(Pval_pred))
    return Pval_pred

popt, pcov = scipy.optimize.curve_fit(HM8_model, xdata=(start, end), ydata=Pval, maxfev=2, p0=p0)

Pval_exp = []
sasaval_avg = []
distval_avg = []
kintval_avg = []
count = 0
for index, row in Pepfrags.iterrows():
    kintval = []
    sasaval = []
    distval = []
    if row['End'] > len(kint):
        continue
    for j in range(int(row['Start']), int(row['End'])):
        kintval.append(kint[j])
        sasaval.append(SASA[j])
        distval.append(DIST[j])
    kintval = np.array(kintval).mean()
    sasaval = np.array(sasaval).mean()
    distval = np.array(distval).mean()
    DUval = Pepfrags_uptake[count]
    P = ((-1*kintval/log(DUval))*time)/(int(row['End']) - int(row['Start']))
    Pval_exp.append(P)
    distval_avg.append(distval)
    sasaval_avg.append(sasaval)
    count = count + 1

np.savetxt('sasaval_avg.txt', sasaval_avg)
np.savetxt('distval_avg.txt', distval_avg)
np.savetxt('Pval_exp.txt', Pval_exp)