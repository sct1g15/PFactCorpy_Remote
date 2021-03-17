import argparse
from PepFrag_predict import PepFrag_uptake, PepFrag_find, PepFrag_predict
from kint import calc_kint_Alt
from math import tanh, exp, log
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('-x', '--excelin', help='Excel Input File', default='20190225 A0201 automated vs manual 2 state.xlsx')
parser.add_argument('-t', '--temp', help='Experimental Temp', default=298)
parser.add_argument('-p', '--pD', help='Experimental pD', default=7)
parser.add_argument('-st', '--state', help='State to analyse (from excel file header)', default='UV_exposed_manual')
args = parser.parse_args()

class ModelTest:
    def __init__(self, inputdata):
        #Set custom parameters for specific model
        self.time = 0.5
        self.HM8_modelparams = [3.28979445e-08, 1.98441207e+00, 8.96847944e+00, -2.58808757e-01]

        #Create empty result objects
        self.Pval_exp = []
        self.Pval_pred = []

        #load simulation histogram data
        self.SASA = np.load('dist_sasa_20bins/SASA_hist.npy', allow_pickle=True)
        self.DIST = np.load('dist_sasa_20bins/DIST_hist.npy', allow_pickle=True)

        #Load experimental peptide fragment data
        ex_indata = pd.read_excel(inputdata)
        Protein_name = ex_indata['Protein'].unique()
        Pepfrags = PepFrag_find(Protein_name[0], args.state, self.time, args.excelin)
        Pepfrags_uptake = PepFrag_uptake(Protein_name[0], args.state, self.time, args.excelin)
        self.fitting_data = Pepfrags_uptake.to_numpy()

        #Calculate intrinsic protection factors
        self.Seq = 'MGSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDGETRKVKAHSQTHRVDLGTLRGYYNQSEAGSHTVQRMYGCDVGSDWR' \
              'FLRGYHQYAYDGKDYIALKEDLRSWTAADMAAQTTKHKWEAAHVAEQLRAYLEGTCVEWLRRYLENGKETLQRTDAPKTHMTHHAVSDHEATLRCWALSFYPAEITLT' \
              'WQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGQEQRYTCHVQHEGLPKPLTLRWE'
        self.kint = []
        for counter, aminoacid in enumerate(self.Seq[0:len(self.Seq)]):
            self.kint.append(calc_kint_Alt(self.Seq[counter], self.Seq[counter - 1], int(args.temp), int(args.pD), counter, len(self.Seq)))

        #Calculate P-factors for experimental peptide fragments
        count = 0
        Pval = []
        for index, row in Pepfrags.iterrows():
            print(row)
            kintval = []
            if row['End'] > len(self.kint):
                continue
            for j in range(int(row['Start']), int(row['End'])):
                kintval.append(self.kint[j])
            kintval_avg = np.array(kintval).mean()
            v1 = (self.fitting_data[count])
            v2 = -1 * (v1 - 1)
            v3 = log(v2) / self.time
            P_factor = -1 * (kintval_avg) / v3
            lnPF = log(P_factor)
            Pval.append(lnPF)
            count += 1
        self.count = count
        self.Pval_exp = Pval

        start = []
        end = []
        for index, row in Pepfrags.iterrows():
            start.append(int(row['Start']))
            end.append(int(row['End']))
        self.start = np.array(start)
        self.end = np.array(end)

    def HM8_method(self):
        betas, gammas, betap, gammap = self.HM8_modelparams
        Start = self.start[0:self.count]
        End = self.end[0:self.count]
        for X in zip(Start, End):
            Start, End = X
            Pval_comb = []
            for j in range(int(Start), int(End)):
                # if kint[j] == -1:
                #     P = 0
                P = 0
                count = 0
                for q in range(0, len(self.SASA[j][0])):
                    for p in range(0, len(self.DIST[j][0])):
                        count += 1
                        P += (((self.SASA[j][0][q] * self.DIST[j][0][p]) / (sum(self.SASA[j][0]) * sum(self.DIST[j][0]))) * (
                                    (betas * (self.SASA[j][1][q]) ** gammas) + (betap * (self.DIST[j][1][p]) ** gammap)))
                # print('Individual Amide ',j,' P-value:',P)
                Pval_comb.append(P)
            print(np.array(Pval_comb).mean())
            self.Pval_pred.append(np.array(Pval_comb).mean())

    def plotgraph(self):
        print(self.Pval_exp)
        print(self.Pval_pred)
        #np.savetxt('Pval_exp.txt', self.Pval_exp)
        #np.savetxt('Pval_pred.txt', self.Pval_pred)
        plt.scatter(self.Pval_pred, self.Pval_exp)
        #corr, _ = pearsonr(self.Pval_pred, self.Pval_exp)
        plt.xlabel('Protection Factor Prediction')
        plt.ylabel('Protection Factor Experimental')
        #plt.annotate("Pearsons R = {:.3f}".format(corr), (0, 1))
        plt.savefig('Protection_Factors_HM8.png')



Fitting = ModelTest(args.excelin)
Fitting.HM8_method()
Fitting.plotgraph()