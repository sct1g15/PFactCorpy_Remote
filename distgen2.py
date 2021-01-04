import MDAnalysis as mda
import numpy as np
from mpi4py import MPI
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--trajin', help='Trajectory xtc file')
parser.add_argument('-s', '--topolin', help='Topology file')
args = parser.parse_args()

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
print(size)
print(rank)

universe = mda.Universe(args.topolin, args.trajin)

if rank == 0:
    ave, res = divmod(len(universe.trajectory), size)
    counts = [ave + 1 if p < res else ave for p in range(size)]
    starts = [sum(counts[:p]) for p in range(size)]
    ends = [sum(counts[:p + 1]) for p in range(size)]
    data = [(starts[p], ends[p]) for p in range(size)]
else:
    data = None

d_loc = comm.scatter(data)

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
    polaratomspos = universe.select_atoms(
        'not type C and not type H and(resname ASP or resname GLU or resname ARG or resname LYS or resname HI[SDE] or '
        'resname ASP or resname GLU or resname SER or resname THR or resname TYR) and not resid ' + str(resid - 1) + '-' + str(resid + 1)).positions

    pairs, dist = mda.lib.distances.capped_distance(amidespos_update, polaratomspos, max_cutoff=12)

    return(dist.min())

def get_amideH_polardistances_Parallel(trajectory, topology, start, end, step=1):
    amidespos = universe.select_atoms('name H').positions

    result = []
    for c, pos in enumerate(amidespos):
        core_result = []
        count = int(c)
        for ts in range(start, end, step):
            core_result.append(get_amideH_polardistances_method(ts, count, pos, trajectory, topology))
        result.append(core_result)
    return(result)

print('rank results:',rank, d_loc, data)


print("After Scatter:")
NewData = get_amideH_polardistances_Parallel(args.trajin, args.topolin, d_loc[0], d_loc[1])


result = comm.gather(NewData)
comm.Barrier()

if rank == 0:
    result = np.array(result)
    result_collate = []
    print(result[0])
    print(result[0][0])
    for i in range(0, result.shape[1]): #360
        amide = []
        for j in range(0, result.shape[0]): #80
            amide.append(result[j][i])
        result_collate.append(amide)
    result_array = np.array(result_collate)
    np.save('dist_result_interim.txt', result_array)

def flatten_list(resid, result_array):
    flattened_list = []
    for x in result_array[resid]:
        for y in x:
            flattened_list.append(y)
    return flattened_list

def plot_dist():
    print('Starting load:')
    distoutput = np.load('dist_result_interim.txt.npy', allow_pickle=True)
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
        distoutput = np.insert(distoutput, i, np.nan, axis=0)

    bins = 500

    plt.hist(flatten_list(90, distoutput), bins=bins)
    plt.title('Distance to nearest polar atom residue 90')
    plt.xlabel('Distance to nearest polar atom (Â)')
    plt.savefig('amide_90_DIST.png')
    plt.clf()

    plt.hist(distoutput[98, :], bins=bins)
    plt.title('Distance to nearest polar atom residue 98')
    plt.xlabel('Distance to nearest polar atom (Â)')
    plt.savefig('amide_98_DIST.png')
    plt.clf()

    plt.hist(distoutput[143, :], bins=bins)
    plt.title('Distance to nearest polar atom residue 143')
    plt.xlabel('Distance to nearest polar atom (Â²)')
    plt.savefig('amide_143_DIST.png')
    plt.clf()

    plt.hist(distoutput[156, :], bins=bins)
    plt.title('Distance to nearest polar atom residue 156')
    plt.xlabel('Distance to nearest polar atom (Â)')
    plt.savefig('amide_156_DIST.png')
    plt.clf()

    plt.hist(distoutput[196, :], bins=bins)
    plt.title('Distance to nearest polar atom residue 196')
    plt.xlabel('Distance to nearest polar atom (Â)')
    plt.savefig('amide_196_DIST.png')
    plt.clf()

    plt.hist(distoutput[203, :], bins=bins)
    plt.title('Distance to nearest polar atom residue 203')
    plt.xlabel('Distance to nearest polar atom (Â)')
    plt.savefig('amide_203_DIST.png')
    plt.clf()

if rank == 0:
    plot_dist()

