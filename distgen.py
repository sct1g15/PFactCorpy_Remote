import MDAnalysis as mda
import mdtraj as md
import numpy as np
import pandas as pd
import distributed
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--trajin', help='Trajectory xtc file')
parser.add_argument('-s', '--topolin', help='Topology file')
args = parser.parse_args()

def get_amideH_polardistances_method(ts, amide, amidepos, trajectory, topology):
    """
    Methodology applied to get_amideH_polardistances_Parallel, not directly callable on its own.
    :param ts: timestep
    :param amide: mdanalysis position object
    :param amidepos: int, residue number
    :return: float, minimum distance for specified timestep and amide residue
    """
    universe = mda.Universe(topology, trajectory)
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

    pairs, dist = mda.lib.distances.capped_distance(amidespos_update, polaratomspos, max_cutoff=12)

    return(dist.min())

def get_amideH_polardistances_Parallel(trajectory, topology, step=1000, procs=40):
    """
    :param trajectory: string, path to trajectory location
    :param topology: string, path to topology location
    :return: np.array, mean distance to nearest polar atom for each residue amide hydrogen
    """
    print('Generating contact maps:######################################')
    universe = mda.Universe(topology, trajectory)

    amidespos = universe.select_atoms('name H').positions
    mean_dist = []
    pool = Pool(processes=procs)

    for c, pos in enumerate(amidespos):
        count = int(c)
        result = []

        result.append(pool.starmap(get_amideH_polardistances_method, [(ts, count, pos, trajectory, topology) for ts in range(0, len(universe.trajectory), step)]))
        mean_dist = np.append(mean_dist, np.max(result))

    np.savetxt('polardistance_meanpara.out', mean_dist)
    return(mean_dist)

get_amideH_polardistances_Parallel(args.trajin, args.topolin, step=1000)

def plot_dist():
    print('Starting load:')
    distoutput_all = np.loadtxt('polardistance_all.out')
    print('Load complete.')

    plt.hist(distoutput_all[:,90], bins=200)
    plt.title('Distance to nearest polar atom residue 90')
    plt.xlabel('Distance to nearest polar atom (Â)')
    plt.savefig('amide_90_dist.png')
    plt.clf()

    plt.hist(distoutput_all[:,98], bins=200)
    plt.title('Distance to nearest polar atom residue 98')
    plt.xlabel('Distance to nearest polar atom (Â)')
    plt.savefig('amide_98_dist.png')
    plt.clf()

    plt.hist(distoutput_all[:,143], bins=200)
    plt.title('Distance to nearest polar atom residue 143')
    plt.xlabel('Distance to nearest polar atom (Â)')
    plt.savefig('amide_143_dist.png')
    plt.clf()

    plt.hist(distoutput_all[:,156], bins=200)
    plt.title('Distance to nearest polar atom residue 156')
    plt.xlabel('Distance to nearest polar atom (Â)')
    plt.savefig('amide_156_dist.png')
    plt.clf()

    plt.hist(distoutput_all[:,196], bins=200)
    plt.title('Distance to nearest polar atom residue 196')
    plt.xlabel('Distance to nearest polar atom (Â)')
    plt.savefig('amide_196_dist.png')
    plt.clf()

    plt.hist(distoutput_all[:,203], bins=200)
    plt.title('Distance to nearest polar atom residue 203')
    plt.xlabel('Distance to nearest polar atom (Â)')
    plt.savefig('amide_203_dist.png')
    plt.clf()

plot_dist()