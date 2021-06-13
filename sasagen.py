import MDAnalysis as mda
import mdtraj as md
import numpy as np
import pandas as pd
import distributed
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--trajin', help='Trajectory xtc file')
parser.add_argument('-s', '--topolin', help='Topology file')
args = parser.parse_args()

def get_sasa(trajectory, topology, coreruns = 10):
    """
    Calculates Solvent Accessible Surface area, currently step size of trajectory is reduced due to memory limits
    :param trajectory: string, path to trajectory location
    :param topology: string, path to topology location
    :param step: int, number of operations to split between serial runs to prevent memory overload
    :return: list, mean SASA for each residue amide hydrogen
    """
    print('Calculating SASA')
    print('step index')
    traj = md.load(trajectory, top=topology)

    topology = traj.topology

    max_frames = int(traj.n_frames)

    mem_max = 15000 #max number of frames per serial run, set due to memory limits
    stepsize = int(max_frames/mem_max)
    print('stepsize:',stepsize)

    if coreruns * mem_max > max_frames:
        print('ERROR')

    sasaoutputavg = pd.DataFrame()
    sasaoutputmax = pd.DataFrame()
    sasaoutput_all = pd.DataFrame()

    for i in range(0, coreruns):
        sasa = md.shrake_rupley(traj[0+i:max_frames:stepsize], probe_radius=0.014, n_sphere_points=960, mode='atom')
        amidesasa = sasa[:, topology.select('name H')]
        print(amidesasa)

        amidesasa_all = pd.DataFrame(amidesasa)
        np.savetxt('amidesasa_all', amidesasa_all)
        amidesasaavg = pd.DataFrame(amidesasa.mean(axis=0), columns=[str(f"sasa{i}_mean")])
        amidesasamax = pd.DataFrame(amidesasa.max(axis=0), columns=[str(f"sasa{i}_max")])

        d = f"sasa{i}"

        sasaoutput_all = pd.concat((sasaoutput_all, amidesasa_all), axis=0)
        sasaoutputavg = pd.concat((sasaoutputavg, amidesasaavg), axis=1)
        sasaoutputmax = pd.concat((sasaoutputmax, amidesasamax), axis=1)

        np.savetxt('amidesasa_all.out', sasaoutput_all)
        sasaoutputavg = sasaoutputavg.mean(axis=1)
        np.savetxt('amidesasaavg_allframe.out', sasaoutputavg)
        sasaoutputmax = sasaoutputmax.mean(axis=1)
        np.savetxt('amidesasamax_allframe.out', sasaoutputmax)

    return(sasaoutputavg, sasaoutputmax)

#get_sasa(args.trajin, args.topolin, coreruns=10)

def plot_sasa():
    print('Starting load:')
    sasaoutput_all = np.loadtxt('amidesasa_all.out')
    sasaoutput_all = sasaoutput_all[~np.isnan(sasaoutput_all)]
    sasaoutput_all = sasaoutput_all.reshape((150003, 360))
    print('Load complete.')

    Seq = 'GSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDGETRKVKAHSQTHRVDLGTLRGYYNQSEAGSHTVQRMYGCDVGSDWR' \
          'FLRGYHQYAYDGKDYIALKEDLRSWTAADMAAQTTKHKWEAAHVAEQLRAYLEGTCVEWLRRYLENGKETLQRTDAPKTHMTHHAVSDHEATLRCWALSFYPAEITLT' \
          'WQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGQEQRYTCHVQHEGLPKPLTLRWE'

    Pro_pos = []
    for count, i in enumerate(Seq):
        if i == 'P':
            Pro_pos.append(count)
        else:
            continue
    print(Pro_pos)

    for i in Pro_pos:
        sasaoutput_all = np.insert(sasaoutput_all, i, 9999, axis=1)

    bins = 25
    axis_range = 0, 0.12
    alpha = 0.4

    plt.hist(sasaoutput_all[:,89], bins=bins)
    plt.title('SASA residue 90')
    plt.xlabel('SASA (Â²)')
    plt.savefig('amide_90_sasa.png')
    plt.clf()

    plt.hist(sasaoutput_all[:,97], bins=bins)
    plt.title('SASA residue 98')
    plt.xlabel('SASA (Â²)')
    plt.savefig('amide_98_sasa.png')
    plt.clf()

    plt.hist(sasaoutput_all[:,142], bins=bins)
    plt.title('SASA residue 143')
    plt.xlabel('SASA (Â²)')
    plt.savefig('amide_143_sasa.png')
    plt.clf()

    plt.hist(sasaoutput_all[:,155], bins=bins)
    plt.title('SASA residue 156')
    plt.xlabel('SASA (Â²)')
    plt.savefig('amide_156_sasa.png')
    plt.clf()

    plt.hist(sasaoutput_all[:,195], bins=bins)
    plt.title('SASA residue 196')
    plt.xlabel('SASA (Â²)')
    plt.savefig('amide_196_sasa.png')
    plt.clf()

    plt.hist(sasaoutput_all[:,202], bins=bins)
    plt.title('SASA residue 203')
    plt.xlabel('SASA (Â²)')
    plt.savefig('amide_203_sasa.png')
    plt.clf()

    fig, ax = plt.subplots()
    ax.hist(sasaoutput_all[:,89], bins=bins, range=axis_range, label = '90', alpha=alpha)
    plt.title('Solvent Accessible Surface Area for HLA-A*02:01')
    plt.xlabel('Solvent Accessible Surface Area (Â²)')
    ax.hist(sasaoutput_all[:,97], bins=bins, range=axis_range, label = '98', alpha=alpha)
    ax.hist(sasaoutput_all[:,142], bins=bins, range=axis_range, label = '143', alpha=alpha)
    ax.hist(sasaoutput_all[:,155], bins=bins, range=axis_range, label = '156', alpha=alpha)
    ax.hist(sasaoutput_all[:,195], bins=bins, range=axis_range, label = '196', alpha=alpha)
    ax.hist(sasaoutput_all[:,202], bins=bins, range=axis_range, label = '203', alpha=alpha)
    leg = ax.legend(title='Residue Amide');
    plt.savefig('amide_collate_SASA.png')
    plt.clf


    hist_out = []
    print(len(Seq))
    print(np.histogram(sasaoutput_all[:, 0], bins=bins, range=axis_range))
    for x in range(0, 275):
        hist_out.append(np.histogram(sasaoutput_all[:, x], bins=bins))

    np.save('SASA_hist.npy', hist_out)


plot_sasa()

