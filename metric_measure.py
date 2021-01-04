import MDAnalysis as mda
import mdtraj as md
import numpy as np
import pandas as pd
import distributed
from multiprocessing import Pool



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

def get_amideH_polardistances_Parallel_V3(trajectory, topology, step=10000):
    """
    :param trajectory: string, path to trajectory location
    :param topology: string, path to topology location
    :return: np.array, mean distance to nearest polar atom for each residue amide hydrogen
    """

    client=distributed.Client()

    print('Generating contact maps:######################################')
    universe = mda.Universe(topology, trajectory)

    amidespos = universe.select_atoms('name H').positions
    mean_dist = []

    for ts in range(0, len(universe.trajectory), step):
        universe.trajectory[ts]

        for c, pos in enumerate(amidespos):
            count = int(c)
            result = []

            amides = universe.select_atoms('name H')
            resid = int(amides[count].resid)
            amidespos_update = universe.select_atoms('name H').positions[count]
            # polaratomspos = universe.select_atoms(
            #     'not backbone and not type C and not type H and(resname ASP or resname GLU or resname ARG or resname LYS or resname HI[SDE] or '
            #     'resname ASP or resname GLU or resname SER or resname THR or resname TYR) and not resid ' + str(resid - 1) + '-' + str(resid + 1)).positions
            polaratomspos = universe.select_atoms(
                'not type C and not type H and(resname ASP or resname GLU or resname ARG or resname LYS or resname HI[SDE] or '
                'resname ASP or resname GLU or resname SER or resname THR or resname TYR) and not resid ' + str(
                    resid - 1) + '-' + str(resid + 1)).positions

            dist = mda.lib.distances.distance_array(amidespos_update, polaratomspos, backend='OpenMP').run()
            result.append(dist)
        mean_dist = np.append(mean_dist, np.max(result))

    np.savetxt('polardistance_meanpara.out', mean_dist)
    return(mean_dist)

def get_amideH_polardistances_Parallel_Capped(trajectory, topology, step=10000):
    """
    :param trajectory: string, path to trajectory location
    :param topology: string, path to topology location
    :return: np.array, mean distance to nearest polar atom for each residue amide hydrogen
    """
    print('Generating contact maps:######################################')
    universe = mda.Universe(topology, trajectory)

    amidespos = universe.select_atoms('name H').positions
    mean_dist = np.empty((0, len(amidespos)), float)

    for ts in range(0, len(universe.trajectory), step):
        universe.trajectory[ts]
        result = []

        for c, pos in enumerate(amidespos):
            count = int(c)

            amides = universe.select_atoms('name H')
            resid = int(amides[count].resid)
            amidespos_update = universe.select_atoms('name H').positions[count]
            # polaratomspos = universe.select_atoms(
            #     'not backbone and not type C and not type H and(resname ASP or resname GLU or resname ARG or resname LYS or resname HI[SDE] or '
            #     'resname ASP or resname GLU or resname SER or resname THR or resname TYR) and not resid ' + str(resid - 1) + '-' + str(resid + 1)).positions
            polaratomspos = universe.select_atoms(
                'not type C and not type H and(resname ASP or resname GLU or resname ARG or resname LYS or resname HI[SDE] or '
                'resname ASP or resname GLU or resname SER or resname THR or resname TYR) and not resid ' + str(
                    resid - 1) + '-' + str(resid + 1)).positions

            contacts, dist = mda.lib.distances.capped_distance(amidespos_update, polaratomspos, max_cutoff=12)
            result.append(dist.min())
        mean_dist = np.vstack((mean_dist, result))
    outputmax = mean_dist.max(axis=0)
    outputavg = mean_dist.mean(axis=0)

    np.savetxt('polardistance_max.out', outputmax)
    np.savetxt('polardistance_mean.out', outputavg)
    return(outputmax, outputavg)

def get_amideH_polardistances_Parallel_V4(trajectory, topology, step=10000):
    """
    :param trajectory: string, path to trajectory location
    :param topology: string, path to topology location
    :return: np.array, mean distance to nearest polar atom for each residue amide hydrogen
    """
    print('Generating contact maps:######################################')
    universe = mda.Universe(topology, trajectory)

    amidespos = universe.select_atoms('name H').positions
    mean_dist = np.empty((0, len(amidespos)), float)

    for ts in range(0, len(universe.trajectory), step):
        universe.trajectory[ts]
        result = []

        for c, pos in enumerate(amidespos):
            count = int(c)

            amides = universe.select_atoms('name H')
            resid = int(amides[count].resid)
            amidespos_update = universe.select_atoms('name H').positions[count]
            # polaratomspos = universe.select_atoms(
            #     'not backbone and not type C and not type H and(resname ASP or resname GLU or resname ARG or resname LYS or resname HI[SDE] or '
            #     'resname ASP or resname GLU or resname SER or resname THR or resname TYR) and not resid ' + str(resid - 1) + '-' + str(resid + 1)).positions
            polaratomspos = universe.select_atoms(
                'not type C and not type H and(resname ASP or resname GLU or resname ARG or resname LYS or resname HI[SDE] or '
                'resname ASP or resname GLU or resname SER or resname THR or resname TYR) and not resid ' + str(
                    resid - 1) + '-' + str(resid + 1)).positions

            dist = mda.lib.distances.distance_array(amidespos_update, polaratomspos, backend='OpenMP')
            result.append(dist.min())
        print(result)
        mean_dist = np.vstack((mean_dist, result))
    print(mean_dist)
    outputmax = mean_dist.max(axis=0)
    outputavg = mean_dist.mean(axis=0)

    np.savetxt('polardistance_max.out', outputmax)
    np.savetxt('polardistance_mean.out', outputavg)
    return(outputmax, outputavg)

def get_sasa(trajectory, topology, step=15000):
    """
    Calculates Solvent Accessible Surface area, currently step size of trajectory is reduced due to memory limits
    :param trajectory: string, path to trajectory location
    :param topology: string, path to topology location
    :return: list, mean SASA for each residue amide hydrogen
    """
    print('Calculating SASA')
    print('step index')
    traj = md.load(trajectory, top=topology)

    topology = traj.topology

    frame_tot = 150000
    #frame_tot = traj.frame
    print(frame_tot)

    cores_tot = int(frame_tot/step)

    if float(cores_tot).is_integer() == False:
        print('ERROR: INCORRECT FRAME TO CORE RATIO, ADJUST STEP TO DIVISIBLE BY FRAME TOTAL i.e. 15000 for 150000 frames'
              'results in 10x seperated operations in serial')

    sasaoutputavg = pd.DataFrame()
    sasaoutputmax = pd.DataFrame()

    for i in range(0, cores_tot):
        sasa = md.shrake_rupley(traj[0+i::cores_tot], probe_radius=0.014, n_sphere_points=960, mode='atom')
        amidesasa = sasa[:, topology.select('name H')]

        amidesasaavg = pd.DataFrame(amidesasa.mean(axis=0), columns=[str(f"sasa{i}")])
        amidesasamax = pd.DataFrame(amidesasa.max(axis=0), columns=[str(f"sasa{i}")])

        print(sasaoutputavg)
        print(amidesasaavg)

        print(i)
        d = f"sasa{i}"
        print(d)
        pd.concat((sasaoutputavg, amidesasaavg), axis=1)
        pd.concat((sasaoutputmax, amidesasamax), axis=1)

    print(sasaoutputavg)
    sasaoutputavg = sasaoutputavg.mean(axis=1)
    sasaoutputmax = sasaoutputmax.mean(axis=1)

    np.savetxt('amidesasaavg.out', sasaoutputavg)
    np.savetxt('amidesasamax.out', sasaoutputmax)

    # amidesasa = sasa[:, topology.select('name H')]
    # amidesasaavg = amidesasa.mean(axis=0)
    # amidesasamax = amidesasa.max(axis=0)
    #np.savetxt('amidesasamax.out', amidesasamax)
    #np.savetxt('amidesasaverage.out', amidesasaavg)

    return(sasaoutputavg, sasaoutputmax)