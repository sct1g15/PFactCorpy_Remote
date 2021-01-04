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


#Input SASA and DIST data
SASA = np.loadtxt('amidesasamax2_allframe.out')
SASA = SASA
DIST = np.loadtxt('polardistance_mean.out')

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

#Checking if large P-values are distorting model fitting, setting all deuterium uptake to effective P-value of 7
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

p0 = 1,1,1,1,1,8,8,4
bounds = ((0, 0, 0, 0, 0, 0, 0, 0), (6, 6, 6, 6, 6, 30, 30, 30))

def HM8_model(X, alpha, beta, gamma, delta, epsilon, xp, yp, zp):
    #print('alpha:',alpha,'beta:',beta,'gamma:',gamma,'delta:',delta,'epsilon:', epsilon,'x,y,z:', xp, yp, zp)
    Start, End = X
    frag_DU = []
    Pval_pred = []
    for i in range(0, len(Start)):
        sa_DU = []
        Pval = []
        if int(End[i]) > len(kint):
            continue
        for j in range(int(Start[i]), int(End[i])):
            # if kint[j] == -1:
            #     P = 0
            P = epsilon * (tanh((alpha * SASA[j])-xp)*tanh((beta * DIST[j])-yp) - gamma* tanh((alpha * SASA[j])-xp) + delta*tanh((beta * DIST[j])-yp))+zp
            if P >= 0:
                Pval.append(P)
                #print(kint[j], P)
                sa_DU.append(1-exp((-kint[j]/P)*time))
            if P < 0:
                sa_DU.append(9999999)
        Pval_pred.append(np.array(Pval).mean())
        frag_DU.append(np.array(sa_DU).mean())
    print(Pval_pred)
    #print(len(Pval_pred))
    return Pval_pred

#z = HM8_model((start, end), 1,1,1,1,1,1,1)

# popt, pcov = scipy.optimize.curve_fit(HM8_model, xdata=(start, end), ydata=Pval, bounds = bounds, p0=p0,
#                                           maxfev=5000)
popt, pcov = scipy.optimize.curve_fit(HM8_model, xdata=(start, end), ydata=Pval, maxfev=5000)

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

alpha, beta, gamma, delta, epsilon, xp, yp, zp = popt
print(epsilon, ' * ((tanh((', alpha, '* x) - ', xp, ')*tanh((', beta, ' * y) - ', yp, ')) - (', gamma,
      ' * tanh((', alpha, ' * x) - ', xp, ')) + (', delta, ' * tanh((', beta, ' * y) - ', yp, ')) + ', zp,
      ')', sep='')
print(popt)
print(pcov)



