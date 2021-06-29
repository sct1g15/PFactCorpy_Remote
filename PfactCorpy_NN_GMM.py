import numpy as np
from math import log
import matplotlib.pyplot as plt
from keras.models import Sequential
from keras.layers import Dense
from keras.callbacks import ModelCheckpoint
from PepFrag_predict import PepFrag_uptake, PepFrag_find, PepFrag_predict
from kint import calc_kint_Alt
import pandas as pd
import argparse
from sklearn.mixture import GaussianMixture
from scipy.stats import norm, mode

#User inputs
time=0.5
parser = argparse.ArgumentParser()
parser.add_argument('-x', '--excelin', help='Excel Input File')
parser.add_argument('-t', '--temp', help='Experimental Temp')
parser.add_argument('-p', '--pD', help='Experimental pD')
parser.add_argument('-st', '--state', help='State to analyse (from excel file header)')
args = parser.parse_args()

#parse DynamX cluster data
#Get kint data
Seq = 'MGSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDGETRKVKAHSQTHRVDLGTLRGYYNQSEAGSHTVQRMYGCDVGSDWR' \
      'FLRGYHQYAYDGKDYIALKEDLRSWTAADMAAQTTKHKWEAAHVAEQLRAYLEGTCVEWLRRYLENGKETLQRTDAPKTHMTHHAVSDHEATLRCWALSFYPAEITLT' \
      'WQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGQEQRYTCHVQHEGLPKPLTLRWE'

B2m_seq = 'MIQRTPKIQVYSRHPAENGKSNFLNCYVSGFHPSDIEVDLLKNGERIEKVEHSDLSFSKDWSFYLLYYTEFTPTEKDEYACRVNHVTLSQPKIVKWDRDM'



def get_kint(Seq, temp, pD):
    kint = []
    for counter, aminoacid in enumerate(Seq[0:len(Seq)]):
        kint.append(calc_kint_Alt(Seq[counter], Seq[counter - 1], int(temp), int(pD), counter, len(Seq)))
    return kint

def get_expPval(Pepfrags, Pepfrags_uptake, kint, prot_name):
    count=0
    Pval = []
    Pval_std = []
    for index, row in Pepfrags.iterrows():
        kintval = []
        if int(row['Start']) == 1 and prot_name == '1A02_HUMAN': #protein specific ommissions e.g. first amide
            continue
        if int(row['Start']) == 2 and prot_name == '1A02_HUMAN': #protein specific ommissions e.g. first amide
            continue
        if int(row['Start']) == 0 and prot_name == 'B2MG_HUMAN': #protein specific ommissions e.g. first amide
            continue
        if row['End'] > len(kint):
            continue
        for j in range(int(row['Start']), int(row['End'])):
            if row['Start'] == 15 and row['Protein'] == 'B2MG_HUMAN':
                print(kint[j])
                print(Pepfrags_uptake[count])
                print((log((Pepfrags_uptake[count] * time))))
            kintval.append(kint[j])
        kintval = np.array(kintval).mean()
        #P = ((-1 * kintval / log(Pepfrags_uptake[count])) * time) / (int(row['End']) - int(row['Start']))
        #P = ((np.mean(-log(kintval))/time)/(int(row['End']) - int(row['Start']))) * log(Pepfrags_uptake[count])
        P = []
        for j in range(int(row['Start']), int(row['End'])):
            if kint[j] == -1:
                continue
            if row['Start'] == 15 and row['Protein'] == 'B2MG_HUMAN':
                print(-(log(kint[j])) / (log((Pepfrags_uptake[count] * time))))
            P.append(-(log(kint[j])) / (log((Pepfrags_uptake[count] * time))))
            #P.append(-log(kint[j])/((log(Pepfrags_uptake[count]))*time))
        Pval.append(np.mean(P))
        # if P > 9:
        #     P = 9
        #Pval.append(P)
        count += 1
    return Pval

def get_pepstartend(Pepfrags):
    start = []
    end = []
    for index, row in Pepfrags.iterrows():
        start.append(int(row['Start']))
        end.append(int(row['End']))
    start = np.array(start)
    end = np.array(end)
    return start, end


#b2m test data
ex_indata = pd.read_excel(args.excelin)
Protein_name = ex_indata['Protein'].unique()
b2m_kint = get_kint(B2m_seq, args.temp, args.pD)
b2m_Pepfrags = PepFrag_find(Protein_name[2], args.state, time, args.excelin) #selects b2m as input name for subset, 1 is nan
b2m_Pepfrags_uptake = PepFrag_uptake(Protein_name[2], args.state, time, args.excelin)
b2m_Pepfrags_uptake = b2m_Pepfrags_uptake.to_numpy()
b2m_Pval = get_expPval(b2m_Pepfrags, b2m_Pepfrags_uptake, b2m_kint, Protein_name[2])
b2m_start, b2m_end = get_pepstartend(b2m_Pepfrags)

#Input training Peptide DU data
ex_indata = pd.read_excel(args.excelin)
Protein_name = ex_indata['Protein'].unique()
kint = get_kint(Seq, args.temp, args.pD)
Pepfrags = PepFrag_find(Protein_name[0], args.state, time, args.excelin)
Pepfrags_uptake = PepFrag_uptake(Protein_name[0], args.state, time, args.excelin)
Pepfrags_uptake = Pepfrags_uptake.to_numpy()
Pval = get_expPval(Pepfrags, Pepfrags_uptake, kint, Protein_name[0])
Start, End = get_pepstartend(Pepfrags)

#Input Histogram Data

DISTdat = np.load('DIST_histtrain.npy', allow_pickle=True)
DISTdat_test = np.load('DIST_histtest.npy', allow_pickle=True)
SASAdat = np.load('SASA_histtrain.npy', allow_pickle=True)
SASAdat_test = np.load('SASA_histtest.npy', allow_pickle=True)

#Combine single amide histograms into averaged peptide fragment histograms

def Pro_adjust(Seq, pos):
    output = pos
    for i in range(0, pos):
        if Seq[i] == 'P':
            pos = pos - 1
    return pos

def Pro_pos(Seq):
    pos = []
    for i in range(0, len(Seq)):
        if Seq[i] == 'P':
            pos.append(pos)
    return pos

def NormalizeData(data):
    return (data / sum(data))

def normalize_histogram(start, end, kint, DISTdat, SASAdat, prot_name, Pro_pos, Seq):
    Normalized_histograms_dist = []
    Normalized_histograms_sasa = []
    adjust_count_1 = 0
    adjust_count_2 = 0
    for x in zip(start, end):
        Start, End = x
        normdist = []
        normsasa = []
        counter = 0
        if Start == 1 and prot_name == '1A02_HUMAN': #protein specific ommissions e.g. first amide
            adjust_count_1 = 1
            continue
        if Start == 2 and prot_name == '1A02_HUMAN': #protein specific ommissions e.g. first amide
            adjust_count_2 = 1
            continue
        if Start == 0 and prot_name == 'B2MG_HUMAN': #protein specific ommissions e.g. first amide
            adjust_count_1 = 1
            continue
        if End >= len(kint): #c-terminal tail ommission
            continue
        for q in range(Start, End):
            if q in Pro_pos:
                continue
            else:
                normdist.append(DISTdat[Pro_adjust(Seq, q)-adjust_count_1-adjust_count_2][0])
                normsasa.append(SASAdat[Pro_adjust(Seq, q)-adjust_count_1-adjust_count_2][0])
                counter += 1
        normdistavg = sum(normdist)/counter
        normsasaavg = sum(normsasa)/counter
        Normalized_histograms_dist.append(NormalizeData(normdistavg))
        Normalized_histograms_sasa.append(NormalizeData(normsasaavg))
    return Normalized_histograms_dist, Normalized_histograms_sasa

print('Seq length (With prolines):',len(b2m_kint))
print('DIST training Data length (No prolines):',len(DISTdat))
print('SASA training Data length (No prolines):',len(SASAdat))
print('DIST test Data length (No prolines):',len(DISTdat_test))
print('SASA test Data length (No prolines):',len(SASAdat_test))

b2m_Pro_pos = Pro_pos(B2m_seq)
Pro_pos = Pro_pos(Seq)
Normalized_histograms_dist, Normalized_histograms_sasa = normalize_histogram(Start, End, kint, DISTdat, SASAdat, Protein_name[0], Pro_pos, Seq)
B2m_Normalized_histograms_dist, B2m_Normalized_histograms_sasa = normalize_histogram(b2m_start, b2m_end, b2m_kint, DISTdat_test, SASAdat_test, Protein_name[2], b2m_Pro_pos, B2m_seq)

Normalized_histograms_train = np.concatenate((Normalized_histograms_dist, Normalized_histograms_sasa), axis=1)
Normalized_histograms_train = Normalized_histograms_train/2 #renormalise to 1 after combination
Pval = np.asarray(Pval)
Normalized_histograms_test = np.concatenate((B2m_Normalized_histograms_dist, B2m_Normalized_histograms_sasa), axis=1)
Normalized_histograms_test = Normalized_histograms_test/2
b2m_Pval = np.asarray(b2m_Pval)

#GMM
#User defined variables
max_n_components = 10
n_repeat_BIC_calc = 10

#Run BIC Model Selection
# tensor_x = np.asarray(range(0, len(Normalized_histograms_train[0])))
# tensor_y = Normalized_histograms_train
tensor_x = np.asarray(range(0, len(Normalized_histograms_dist[0])))
tensor_y_dist = Normalized_histograms_dist
tensor_y_sasa = Normalized_histograms_sasa

def BIC_GMM_model_selection(tensor_y):
    amide_count = 0
    tensor_x = np.asarray(range(0, len(tensor_y[0])))
    all_models = []
    for amide in tensor_y:
        if amide_count == 2:
            break
        else:
            amide_count += 1
            print(amide_count)
            line_distribution   = np.random.choice(a = tensor_x, size = 100000, p = amide)
            number_points       = len(line_distribution)
            n_components = np.arange(1, max_n_components)

            component_n = []
            for q in range(0, n_repeat_BIC_calc):
                models = [GaussianMixture(n_components=n, covariance_type='full').fit(
                    np.reshape(line_distribution, (number_points, 1)))
                          for n in n_components]
                y = [m.bic(np.reshape(line_distribution, (number_points, 1))) for m in models]
                gradient = [y[n] - y[n - 1] for n in np.arange(1, max_n_components - 1)]
                component_n.append([gradient.index(min(gradient, key=abs))])
            all_models.append([models[round(np.amax(component_n))]])
    return all_models, max(component_n)

def plot_gmm(model, tensor_y, name):
    amide_count = 0
    for i, q in zip(model, tensor_y):
        gauss_mixt = []
        if amide_count == 2:
            break
        for x in zip(i[0].means_, i[0].covariances_, i[0].weights_):
            means, sd, weight = x
            gmm_x = np.linspace(0, len(Normalized_histograms_dist[0]), 1000)
            gauss_mixt.append(weight * norm.pdf(gmm_x, means[0], np.sqrt(sd[0][0])))
        gauss_mixt_t = np.sum(gauss_mixt, axis=0)
        fig, axis = plt.subplots(1, 1, figsize=(10, 12))
        axis.hist(tensor_x, weights=tensor_y[amide_count], bins=len(tensor_x))
        axis.set_xlabel('Input Tensor Histogram Datapoint')
        axis.set_ylabel('Normalised Probability')
        title = 'Peptide {} Gaussian Mixture Model Curve Fitting. {}'.format(str(amide_count), str(name))
        axis.set_title(title)
        plt.show()
        for i in range(len(gauss_mixt)):
            plt.plot(gmm_x, gauss_mixt[i], label='Gaussian ' + str(i + 1))
        axis.plot(gmm_x, gauss_mixt_t, label='GMM Model', color='black')
        axis.legend()
        amide_count += 1

    plt.show()

def combine_models(model1, model1_max, model2, model2_max):
    for i, q in zip(model1, model2):
        print(np.hstack(i[0].means_), np.hstack(i[0].covariances_), np.hstack(i[0].weights_))

# models_dist, max_n_dist = BIC_GMM_model_selection(tensor_y_dist)
# models_sasa, max_n_sasa = BIC_GMM_model_selection(tensor_y_sasa)
# combine_models(models_dist, max_n_dist, models_sasa, max_n_sasa)
# plot_gmm(models_dist, tensor_y_dist, 'Max N Components')
# plot_gmm(models_sasa, tensor_y_sasa, 'Max N Components')

#Kernel density estimation
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import LeaveOneOut


tensor_x = np.asarray(range(0, len(Normalized_histograms_dist[0])))
tensor_y_dist = Normalized_histograms_dist
tensor_y_sasa = Normalized_histograms_sasa

bandwidths = 10 ** np.linspace(-1, 1, 100)
grid = GridSearchCV(KernelDensity(kernel='gaussian'),
                    {'bandwidth': bandwidths},
                    cv=LeaveOneOut())
grid.fit([tensor_x, tensor_y_dist[0]]);
print(grid.best_params_)

fig, ax = plt.subplots()
kde = KernelDensity(kernel=kernel, bandwidth=0.5).fit([tensor_x, tensor_y_dist[0]])
log_dens = kde.score_samples(X_plot)
ax.plot(tensor_x, np.exp(log_dens), color=color, lw=lw,
        linestyle='-', label="kernel = '{0}'".format(kernel))