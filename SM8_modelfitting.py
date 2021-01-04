# -*- coding: utf-8 -*-

import sys
import numpy as np
import os
import subprocess
import re
import scipy
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import math
import mdtraj as md
import MDAnalysis as mda
import argparse
from matplotlib import animation
from multiprocessing import Pool
import time
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP

import pylab
from matplotlib import rc, rcParams
from mpl_toolkits.mplot3d import Axes3D

# parser = argparse.ArgumentParser()
# parser.add_argument('-f', '--trajin', help='Trajectory xtc file')
# parser.add_argument('-s', '--topolin', help='Topology file')
# parser.add_argument('-f2', '--trajin2', help='Trajectory xtc file')
# parser.add_argument('-s2', '--topolin2', help='Topology file')
# args = parser.parse_args()

def get_viable_clusters():
    """
    Only required for selection of residues present in HDX after clustering
    :return:
    """
    filenames = [f for f in os.listdir('/scratch/sct1g15/mydocuments/RPC/exPfact/A0201_20180903/default_state/output/singlep_contig/') if f.endswith('.mclust')]
    meanpos = np.empty((1,2))
    for contig in filenames:
        array = re.findall(r'[0-9]+', contig)
        data = np.loadtxt('/scratch/sct1g15/mydocuments/RPC/exPfact/A0201_20180903/default_state/output/singlep_contig/' + contig, comments='&')
        mean = []
        stddev = []
        position = []
        counter = 0
        for i in range(int(np.amin(data[:, 0])),(int(np.amax(data[:, 0]))+1)):
            b = data[data[:, 0] == i]
            mean = np.append(mean, np.mean(b[:, 1]))
            position = np.append(position, range(int(array[0]), int(array[1])+1)[counter])
            stddev = np.append(stddev, np.std(b[:, 1]))
            counter=counter+1
        meanposcontig = np.column_stack((position, mean))
        meanpos = np.concatenate((meanpos, meanposcontig))
    mean = np.delete(meanpos, (0), axis=0)
    mean = mean[np.argsort(mean[:, 0])]
    mean[:, 0] -= 1   #due to M starting residue in HDX data (not present in pdb structures)
    return(mean[:,0])

def get_hbonds_all():
    """
    Only required for selection of residues present in HDX after clustering
    :return:
    """
    data = np.loadtxt('/scratch/sct1g15/mydocuments/RPC/amber99ff_errorsims/3PWL_A0201/hbonds/hbonds.txt')
    return data

def get_sasa(trajectory, topology, step=100):
    """
    Calculates Solvent Accessible Surface area, currently step size of trajectory is reduced due to memory limits
    :param trajectory: string, path to trajectory location
    :param topology: string, path to topology location
    :return: list, mean SASA for each residue amide hydrogen
    """
    print('Calculating SASA')
    traj = md.load(trajectory, top=topology)

    topology = traj.topology

    sasa = md.shrake_rupley(traj[::step], probe_radius=0.014, n_sphere_points=960, mode='atom')

    amidesasa = sasa[:, topology.select('name H')]
    #np.savetxt('amidesasa.out', amidesasa)

    amidesasasum = amidesasa.sum(axis=0)
    #np.savetxt('amidesasasum.out', amidesasasum)
    amidesasamean = amidesasa.mean(axis=0)
    #np.savetxt('amidesasaverage.out', amidesasamean)
    amidesasamax = amidesasa.max(axis=0)

    return(amidesasamax)

    # plt.plot(traj[0:150002:1000].time, total_sasa)
    # plt.xlabel('Time [ps]', size=16)
    # plt.ylabel('Total SASA (nm)^2', size=16)
    # plt.savefig('sasa.png')

    # np.savetxt('sasa_total.out', sasa, delimiter=',')
    # np.savetxt('sasa_total.out', total_sasa, delimiter=',')

def get_amideH_polardistances():
    print('Generating contact maps:')
    universe = mda.Universe('/scratch/sct1g15/mydocuments/RPC/3PWL_A0201/500ns_con/md_0_1_protein.pdb',
                            '/scratch/sct1g15/mydocuments/RPC/3PWL_A0201/500ns_con/md_0_1_protein_molurcenter_uratom_urfit_1000skip.xtc')
    amides = universe.select_atoms('name H')
    amidespos = universe.select_atoms('name H').positions

    counter = 0
    all_dist = np.array([]).reshape(0, len(universe.trajectory))
    mean_dist = []
    for i in amidespos:
        counter2 = 0
        y = []
        for j in universe.trajectory:
            universe.trajectory[counter2]
            amidespos_update = universe.select_atoms('name H').positions[counter]
            counter2 += 1
            resid = int(amides[counter].resid)
            polaratomspos = universe.select_atoms(
                'not backbone and not type C and not type H and(resname ASP or resname GLU or resname ARG or resname LYS or resname HI[SDE] or '
                'resname ASP or resname GLU or resname SER or resname THR or resname TYR) and not resid '+str(resid-1)+'-'+str(resid+1)).positions
            dist = mda.lib.distances.distance_array(amidespos_update, polaratomspos)
            y.append(dist.min())
        counter += 1
        all_dist = np.vstack((all_dist, y))
        mean_dist = np.append(mean_dist, np.mean(y))


    np.savetxt('polardistance_mean.out', mean_dist)

def get_amideH_polardistances_Parallel(trajectory, topology, step=100):
    """
    :param trajectory: string, path to trajectory location
    :param topology: string, path to topology location
    :return: np.array, mean distance to nearest polar atom for each residue amide hydrogen
    """
    print('Generating contact maps:######################################')
    universe = mda.Universe(topology, trajectory)

    amidespos = universe.select_atoms('name H').positions
    mean_dist = []
    pool = Pool(processes=40)

    for c, pos in enumerate(amidespos):
        count = int(c)
        result = []

        result.append(pool.starmap(get_amideH_polardistances_method, [(ts, count, pos) for ts in range(0, len(universe.trajectory), step)]))
        mean_dist = np.append(mean_dist, np.max(result))

    #np.savetxt('polardistance_meanpara.out', mean_dist)
    return(mean_dist)

def get_amideH_polardistances_method(ts, amide, amidepos):
    """
    Methodology applied to get_amideH_polardistances_Parallel, not directly callable on its own.
    :param ts: timestep
    :param amide: mdanalysis position object
    :param amidepos: int, residue number
    :return: float, minimum distance for specified timestep and amide residue
    """
    universe = mda.Universe(args.topolin, args.trajin)
    universe.trajectory[ts]

    amides = universe.select_atoms('name H')
    resid = int(amides[amide].resid)
    amidespos_update = universe.select_atoms('name H').positions[amide]
    # polaratomspos = universe.select_atoms(
    #     'not backbone and not type C and not type H and(resname ASP or resname GLU or resname ARG or resname LYS or resname HI[SDE] or '
    #     'resname ASP or resname GLU or resname SER or resname THR or resname TYR) and not resid ' + str(resid - 1) + '-' + str(resid + 1)).positions
    polaratomspos = universe.select_atoms(
        'not type C and not type H and(resname ASP or resname GLU or resname ARG or resname LYS or resname HI[SDE] or '
        'resname ASP or resname GLU or resname SER or resname THR or resname TYR) and not resid ' + str(resid - 1) + '-' + str(resid + 1)).positions

    dist = mda.lib.distances.distance_array(amidespos_update, polaratomspos)

    return(dist.min())

def get_cluster_results(pos):
    filenames = [f for f in os.listdir('/scratch/sct1g15/mydocuments/RPC/exPfact/A0201_20180903/default_state/output/singlep_contig/') if f.endswith('.mclust')]
    meanpos = np.empty((1,2))
    for contig in filenames:
        array = re.findall(r'[0-9]+', contig)
        data = np.loadtxt('/scratch/sct1g15/mydocuments/RPC/exPfact/A0201_20180903/default_state/output/singlep_contig/' + contig, comments='&')
        mean = []
        stddev = []
        position = []
        counter = 0
        for i in range(int(np.amin(data[:, 0])),(int(np.amax(data[:, 0]))+1)):
            b = data[data[:, 0] == i]
            mean = np.append(mean, np.mean(b[:, 1]))
            position = np.append(position, range(int(array[0]), int(array[1])+1)[counter])
            stddev = np.append(stddev, np.std(b[:, 1]))
            counter=counter+1
        meanposcontig = np.column_stack((position, mean))
        meanpos = np.concatenate((meanpos, meanposcontig))
    mean = np.delete(meanpos, (0), axis=0)
    mean = mean[np.argsort(mean[:, 0])]
    mean[:, 0] -= 1 #due to M starting residue in HDX data (not present in pdb structures)
    return(mean[mean[:, 0] == pos][0][1])

def PF_function_M8(X , baseSASA, expSASA, baseDist, expDist):
    SASA, Dist = X
    return ((baseSASA)*(SASA)**(expSASA*-1))+((baseDist)*(Dist)**(expDist*-1))

def get_max(x, y):
    for i, j in zip(x, y):
        if i > j:
            return i
        if j > i:
            return j
        if i == j:
            return i

def PF_function_M8_sigmoid(X , baseSASA, expSASA, scaleSASA, baseDist, expDist, scaleDist):
    SASA, Dist = X
    sasa = scaleSASA / (1.0 + np.exp(baseSASA*(SASA-expSASA)))
    dist = scaleDist / (1.0 + np.exp(baseDist*(Dist-expDist)))
    return(np.maximum(sasa, dist))

def PF_function_M8_sigmoid_SASA(X , baseSASA, expSASA, scaleSASA):
    sasa = scaleSASA / (1.0 + np.exp(baseSASA*(X-expSASA)))
    return sasa

def PF_function_M8_sigmoid_DIST(X , baseDist, expDist, scaleDist):
    dist = scaleDist / (1.0 + np.exp(baseDist*(X-expDist)))
    return dist

def get_base_parameter_manualHDX_sigmoid(trajectory1, topology1, trajectory2, topology2, sasastep=100, diststep=100, preload=False):

    if preload == False:
        ts = time.time()
        print(ts)
        SASA = get_sasa(trajectory1, topology1, sasastep)
        ts = time.time()
        print(ts)
        Dist = get_amideH_polardistances_Parallel(trajectory1, topology1, diststep)
        ts = time.time()
        print(ts)
        SASA2 = get_sasa(trajectory2, topology2, sasastep)
        ts = time.time()
        print(ts)
        Dist2 = get_amideH_polardistances_Parallel(trajectory2, topology2, diststep)
        ts = time.time()
        print(ts)

        fit_index_full=[8, 9, 30, 61, 84, 156, 140, 168, 169, 231, 242, 8, 9, 30, 61, 84, 156, 140, 168, 169, 231, 242] #GSH indexing
        fit_index = [8, 9, 30, 61, 84, 156, 140, 168, 169, 231, 242]
        #default_30s_manual_PF = [11.51, 11.51, 6.56, 11.37, 3.97, 11.51, 7.05, 4.36, 8.34, 3.52, 8.69]
        #uv_exposed_30s_manual_PF = [11.51, 11.51, 6.58, 7.11, 3.74, 11.51, 11.51, 4.26, 11.51, 3.30, 7.19]

        inddataSASA = []
        inddataDIST = []

        for i in fit_index:
            inddataSASA.append(SASA[int(i) - 1])
            inddataDIST.append(Dist[int(i) - 1])

        for i in fit_index:
            inddataSASA.append(SASA2[int(i) - 1])
            inddataDIST.append(Dist2[int(i) - 1])

        inddataDIST = np.array(inddataDIST)
        inddataSASA = np.array(inddataSASA) * 100

        np.savetxt('SASAdata_preload.txt', inddataSASA)
        np.savetxt('Distdata_preload.txt', inddataDIST)
    if preload == True:
        inddataSASA = np.loadtxt('SASAdata_preload.txt')
        inddataDIST = np.loadtxt('Distdata_preload.txt')

    default_30s_manual_PF = [11.37, 11.37, 6.56, 11.37, 3.97, 11.37, 7.05, 4.36, 8.34, 3.52, 8.69, 11.37, 11.37, 6.58,
                             7.11, 3.74, 11.37, 11.37, 4.26, 11.37, 3.30, 7.19]
    fit_index = [8, 9, 30, 61, 84, 156, 140, 168, 169, 231, 242]
    fit_index_full = [8, 9, 30, 61, 84, 156, 140, 168, 169, 231, 242, 8, 9, 30, 61, 84, 156, 140, 168, 169, 231, 242]

    ###PLOT 3D Datapoints###
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    print(inddataDIST)
    print(inddataSASA)
    print(default_30s_manual_PF)
    ax.scatter(inddataSASA, inddataDIST, default_30s_manual_PF, alpha=1)
    ax.set_ylabel('Distance (Å)')
    ax.set_xlabel('SASA (Å²)')
    ax.set_zlabel('ln(PF)')
    plt.savefig('SASAmodel_scatter.png')

    # for i, txt in enumerate(fit_index_full):
    #     ax.annotate(txt, xy=(inddataSASA[i], default_30s_manual_PF[i]))

    def init():
        ax.scatter(inddataSASA, inddataDIST, default_30s_manual_PF)
        return fig,

    def animate(i):
        ax.view_init(elev=10., azim=i)
        return fig,

    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=360, interval=20, blit=True)
    anim.save('3dscatter_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
    ####################

    popt_sasa, pcov_sasa = scipy.optimize.curve_fit(PF_function_M8_sigmoid_SASA, xdata=inddataSASA, ydata=default_30s_manual_PF,
                                          maxfev=50000)
    popt_dist, pcov_dist = scipy.optimize.curve_fit(PF_function_M8_sigmoid_DIST, xdata=inddataDIST, ydata=default_30s_manual_PF,
                                          maxfev=50000)

    popts = []
    for i in range(0, 3):
        popts.append(popt_sasa[i])
    for i in range(0, 3):
        popts.append(popt_dist[i])
    np.savetxt('popt.output', popts)

    print('P-OPT_sasa:')
    print(popt_sasa)
    print('P-OPT_dist:')
    print(popt_dist)
    #perr = np.sqrt(np.diag(pcov))

    ###PLOT Model###
    x = np.linspace(4, 12 , 30)
    y = np.linspace(0, 8, 30)
    X, Y = np.meshgrid(x, y)

    def g(x, y):
        sasa = popt_sasa[2] / (1.0 + np.exp(popt_sasa[0]*(x-popt_sasa[1])))
        dist = popt_dist[2] / (1.0 + np.exp(popt_dist[0]*(y-popt_dist[1])))
        return (np.maximum(sasa, dist))
    Z2 = g(X, Y)
    ax.plot_surface(X, Y, Z2)

    def init():
        ax.plot_surface(X, Y, Z2)
        return fig,

    def animate(i):
        ax.view_init(elev=10., azim=i)
        return fig,

    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=360, interval=20, blit=True)
    anim.save('3dmodel_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
    plt.savefig('SASAmodel_3D.png')

    coords = []
    for a, b, c in zip(X, Y, Z2):
        for a1, b1, c1 in zip(a, b, c):
            coords.append((a1, b1, c1))
    np.savetxt('plane_coords.csv', coords, delimiter=',')

    coords2 = []
    for a, b, c, d in zip(inddataSASA, inddataDIST, default_30s_manual_PF, fit_index_full):
        coords2.append((a, b, c, d))
    np.savetxt('scatter_coords.csv', coords2, delimiter=',')

    plt.figure()
    fig, ax = plt.subplots()
    plt.plot(inddataDIST, default_30s_manual_PF, 'bo')
    ax.set_xlabel('Distance (Å)')
    ax.set_ylabel('ln(PF)')
    for i, txt in enumerate(fit_index):
        ax.annotate(txt, xy=(inddataDIST[i], default_30s_manual_PF[i]))
    def q(x):
        return popt_dist[2]*(11.51 / (1.0 + np.exp(popt_dist[0]*(x-popt_dist[1]))))
    Z3 = q(x)
    plt.plot(x, Z3)
    plt.savefig('DISTPLOT.png')


    plt.figure()
    fig, ax = plt.subplots()
    plt.plot(inddataSASA, default_30s_manual_PF, 'bo')
    ax.set_xlabel('SASA (Å²)')
    ax.set_ylabel('ln(PF)')
    for i, txt in enumerate(fit_index):
        ax.annotate(txt, xy=(inddataSASA[i], default_30s_manual_PF[i]))
    def q(x):
        return popt_sasa[2]*(11.51 / (1.0 + np.exp(popt_sasa[0]*(x-popt_sasa[1]))))
    Z3 = q(y)
    plt.plot(y, Z3)
    plt.savefig('SASAPLOT.png')



    ###Applying SM8 model prediction###

    SASA = np.array(SASA) * 100
    SASA2 = np.array(SASA2) * 100

    PFpred_def = []
    for i in range(0, 275):
        PFpred_def.append(g(SASA[i], Dist[i]))
    np.savetxt('PFpred_def.output', PFpred_def)
    PFpred_UV = []
    for i in range(0, 275):
        PFpred_UV.append(g(SASA2[i], Dist2[i]))
    np.savetxt('PFpred_UV.output', PFpred_UV)

#get_base_parameter_manualHDX_sigmoid(args.trajin, args.topolin, args.trajin2, args.topolin2, sasastep=10, diststep=50, preload=True)

def plot_secondary():
    default_30s_manual_PF = [11.37, 11.37, 6.56, 11.37, 3.97, 11.37, 7.05, 4.36, 8.34, 3.52, 8.69, 11.37, 11.37, 6.58,
                             7.11, 3.74, 11.37, 11.37, 4.26, 11.37, 3.30, 7.19]
    fit_index_full = [8, 9, 30, 61, 84, 156, 140, 168, 169, 231, 242, 8, 9, 30, 61, 84, 156, 140, 168, 169, 231, 242]

    inddataSASA = np.loadtxt('SASAdata_preload.txt')
    inddataDIST = np.loadtxt('Distdata_preload.txt')

    x = np.linspace(4, 12 , 100)
    y = np.linspace(0, 10, 100)
    X, Y = np.meshgrid(x, y)
    def g(x, y):
        sasa = 11.6408 / (1.0 + np.exp(1.2054*(x-10.0026)))
        dist = 9.5150 / (1.0 + np.exp(1.1609*(y-6.8153)))
        return (np.maximum(sasa, dist))
    Z2 = g(X, Y)

    p = PDBParser()
    structure = p.get_structure('3PWL', '3PWL_asym_rot_cap_htmd_check_amb_2.pdb')
    model = structure[0]
    dssp = DSSP(model, '3PWL_asym_rot_cap_htmd_check_amb_2.pdb')
    print(dssp.keys())

    dssp_output = []
    for i in range(0, 277):
        a_key = list(dssp.keys())[i]
        dssp_output.append(dssp[a_key][2])
    print(dssp_output)

    dssp_output_fit = []
    for i in fit_index_full:
        dssp_output_fit.append(dssp_output[i-1])
    dssp_output_fit = np.array(dssp_output_fit)
    print(dssp_output_fit)


    #np.savetxt('dsspoutput.txt', dssp_output)

    coords = []
    for a, b, c in zip(X, Y, Z2):
        for a1, b1, c1 in zip(a, b, c):
            coords.append((a1, b1, c1))
    np.savetxt('plane_coords.csv', coords, delimiter=',')

    coords2 = []
    for a, b, c, d, e in zip(inddataSASA, inddataDIST, default_30s_manual_PF, fit_index_full, dssp_output_fit):
        coords2.append((a, b, c, d, e))
    np.savetxt('scatter_coords.csv', coords2, delimiter=',', fmt='%s')

#plot_secondary()

def M8_model():
    inddataSASA = np.loadtxt('SASAdata_preload.txt')
    inddataDIST = np.loadtxt('Distdata_preload.txt')

    default_30s_manual_PF = [11.37, 11.37, 6.56, 11.37, 3.97, 11.37, 7.05, 4.36, 8.34, 3.52, 8.69, 11.37, 11.37, 6.58,
                             7.11, 3.74, 11.37, 11.37, 4.26, 11.37, 3.30, 7.19]

    popt, pcov = scipy.optimize.curve_fit(PF_function_M8, xdata=(inddataSASA, inddataDIST),
                                                    ydata=default_30s_manual_PF,
                                                    maxfev=50000)
    print(popt)

    x = np.linspace(5, 12 , 100)
    y = np.linspace(2, 10, 100)
    X, Y = np.meshgrid(x, y)
    def gM8(x, y):
        SASA = x
        Dist = y
        return ((popt[0])*(SASA)**(popt[1]*-1))+((popt[2])*(Dist)**(popt[3]*-1))

    Z2 = gM8(X, Y)

    coords=[]
    for a, b, c in zip(X, Y, Z2):
        for a1, b1, c1 in zip(a, b, c):
            coords.append((a1, b1, c1))
    np.savetxt('M8_plane_coords.csv', coords, delimiter=',')

    coords2=[]
    for a, b, c in zip(inddataSASA, inddataDIST, default_30s_manual_PF):
        coords2.append((a, b, c))
    np.savetxt('M8_scatter_coords.csv', coords2, delimiter=',', fmt='%s')

    def fM8(x, y):
        SASA = x
        Dist = y
        return ((1.3e-3)*(SASA)**(2.64*-1))+((3.65e1)*(Dist)**(1.27*-1))

    coords3 = []
    Z3 = fM8(X, Y)
    for a, b, c in zip(X, Y, Z3):
        for a1, b1, c1 in zip(a, b, c):
            coords3.append((a1, b1, c1))
    np.savetxt('M8_theirparams_plane_coords.csv', coords3, delimiter=',')

M8_model()