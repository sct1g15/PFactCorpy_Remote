import MDAnalysis as mda
import numpy as np
from mpi4py import MPI
import matplotlib.pyplot as plt
import argparse

# result_array = np.load('dist_result_interim.txt.npy', allow_pickle=True)
# print(result_array.shape[0])
# print(result_array.shape[1])
#
# # print(result_array[94])
# flattened_list = []
# for x in result_array[149]:
#     for y in x:
#         flattened_list.append(y)
# print(flattened_list)
#
# bins=500
# plt.hist(flattened_list, bins=bins)
# plt.savefig('amide_149_DIST.png')

def flatten_list(resid, result_array):
    flattened_list = []
    for x in result_array[resid]:
        if x == 9999:
            flattened_list.append(9999)
        else:
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
        distoutput = np.insert(distoutput, i, 9999, axis=0)

    bins = 50
    axis_range = 0, 10
    alpha = 0.4

    plt.hist(flatten_list(89, distoutput), bins=bins, range=axis_range)
    plt.title('Distance to nearest polar atom residue 90')
    plt.xlabel('Distance to nearest polar atom (Â)')
    plt.savefig('amide_90_DIST.png')
    plt.clf()

    plt.hist(flatten_list(97, distoutput), bins=bins, range=axis_range)
    plt.title('Distance to nearest polar atom residue 98')
    plt.xlabel('Distance to nearest polar atom (Â)')
    plt.savefig('amide_98_DIST.png')
    plt.clf()

    plt.hist(flatten_list(142, distoutput), bins=bins, range=axis_range)
    plt.title('Distance to nearest polar atom residue 143')
    plt.xlabel('Distance to nearest polar atom (Â²)')
    plt.savefig('amide_143_DIST.png')
    plt.clf()

    plt.hist(flatten_list(155, distoutput), bins=bins, range=axis_range)
    plt.title('Distance to nearest polar atom residue 156')
    plt.xlabel('Distance to nearest polar atom (Â)')
    plt.savefig('amide_156_DIST.png')
    plt.clf()

    plt.hist(flatten_list(195, distoutput), bins=bins, range=axis_range)
    plt.title('Distance to nearest polar atom residue 196')
    plt.xlabel('Distance to nearest polar atom (Â)')
    plt.savefig('amide_196_DIST.png')
    plt.clf()

    plt.hist(flatten_list(202, distoutput), bins=bins, range=axis_range)
    plt.title('Distance to nearest polar atom residue 203')
    plt.xlabel('Distance to nearest polar atom (Â)')
    plt.savefig('amide_203_DIST.png')
    plt.clf()
    #
    # fig, ax = plt.subplots()
    # ax.hist(flatten_list(89, distoutput), bins=bins, range=axis_range, label = '90', alpha=alpha)
    # plt.title('Distance to nearest polar atom for HLA-A*02:01')
    # plt.xlabel('Distance to nearest polar atom (Â)')
    # ax.hist(flatten_list(97, distoutput), bins=bins, range=axis_range, label = '98', alpha=alpha)
    # ax.hist(flatten_list(142, distoutput), bins=bins, range=axis_range, label = '143', alpha=alpha)
    # ax.hist(flatten_list(155, distoutput), bins=bins, range=axis_range, label = '156', alpha=alpha)
    # ax.hist(flatten_list(195, distoutput), bins=bins, range=axis_range, label = '196', alpha=alpha)
    # ax.hist(flatten_list(202, distoutput), bins=bins, range=axis_range, label = '203', alpha=alpha)
    # leg = ax.legend(title='Residue Amide');
    # plt.savefig('amide_collate_DIST.png')

    hist_out = []
    for x in range(0, 275):
        hist_out.append(np.histogram(flatten_list(x, distoutput), bins=bins, range=axis_range))

    np.save('DIST_hist.npy', hist_out)



plot_dist()