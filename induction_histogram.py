#!/usr/bin/env python
import argparse
import collections
import itertools
import json
from color_palettes import palette
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np


def generate_plots(data_len_tuple_dict, axis, name, title):
    """makes pdf scatter plot"""
    # set up colors and names
    almost_gray = '#808080'
    almost_black = '#262626'
    # color_set = palette[len(data_len_tuple_dict)]
    color_set = {
        'CmeR': '#ab62c0',
        'DesT': '#50ab6d',
        'NalC': '#c95573',
        'PmeR': '#929d3d',
        'SmeT': '#648ace',
        'TtgR': '#ca743e'
    }
    log_axes = ['repressed', 'induced', 'free']
    normal_axes = ['fold_induction', 'fold_repression']

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    # max_val = -1
    # min_val = 1000000000
    # for i, length in enumerate(data_len_tuple_dict.keys()):
    #     temp_max = max(data_len_tuple_dict[length])
    #     temp_min = min(data_len_tuple_dict[length])
    #     if temp_max > max_val:
    #         max_val = temp_max
    #     if temp_min < min_val:
    #         min_val = temp_min
    # bins = np.logspace(np.log10(min_val), np.log10(max_val), num=20)
    # for i, length in enumerate(data_len_tuple_dict.keys()):
    #     if len(data_len_tuple_dict[length]) == 0:
    #         continue
    #     else:
    #         data_list = data_len_tuple_dict[length]
    #         ax1.hist(data_list, bins, alpha=0.5, color=color_set[i], label=str(length))
    data_list = []
    for single_list in data_len_tuple_dict.values():
        data_list.extend(single_list)
    bins = np.linspace(min(data_list), max(data_list))
    ax1.hist(data_list, bins, alpha=0.75, color=color_set[title], histtype='stepfilled',
             label='n = {0}'.format(len(data_list)))
    if axis in log_axes:
        plt.xscale('log')
    legend = ax1.legend(loc='best', scatterpoints=1, framealpha=0.5)
    rect = legend.get_frame()
    rect.set_linewidth(0.0)
    texts = legend.texts
    for t in texts:
        t.set_color(almost_black)

    spines_to_remove = ['top', 'right']
    for spine in spines_to_remove:
        ax1.spines[spine].set_visible(False)
    ax1.xaxis.set_ticks_position('none')
    ax1.yaxis.set_ticks_position('none')
    spines_to_keep = ['bottom', 'left']
    for spine in spines_to_keep:
        ax1.spines[spine].set_linewidth(0.5)
        ax1.spines[spine].set_color(almost_black)
    ax1.xaxis.label.set_color(almost_black)
    ax1.yaxis.label.set_color(almost_black)

    plt.title(title, fontsize=20, y=1.02)
    plt.xlabel(axis, fontsize=15)
    plt.ylabel('frequency', fontsize=15)
    plt.savefig('{0}'.format(name))
    plt.close()

    q25, q75 = np.percentile(data_list, [25, 75])
    iqr = q75 - q25
    top_whisker = q75 + 1.5 * iqr
    import sys
    sys.stdout.write('q25: {0}\tq75: {1}\tiqr: {2}\t tw: {3}\n'.format(q25, q75, iqr, top_whisker))
    # if axes[0] in log_axes:
    #     plt.xscale('log')
    #     ax1.set_xbound(lower=100.0)
    # else:
    #     ax1.set_xbound(lower=0.0)
    # if axes[1] in log_axes:
    #     plt.yscale('log')
    #     ax1.set_ybound(lower=100.0)
    # else:
    #     ax1.set_ybound(lower=0.0)
    # plt.xlabel('{0}'.format(axes[0]), fontsize=20)
    # plt.ylabel(axes[1], fontsize=20)
    #
    # pdf.savefig(fig)
    # plt.close()


def get_data_by_column(axis, data_file, seq_csi_score_dict):
    possible_lengths = [16, 17, 18, 19]
    indices = {}
    indices_found = False
    length_tuple_dict = collections.OrderedDict()
    for length in possible_lengths:
        length_tuple_dict[length] = []

    with open(data_file, 'r') as f:
        while not indices_found:
            line = f.readline().rstrip()
            if axis in line:
                split_line = line.split()
                x_index = split_line.index(axis) + 1
                indices[axis] = x_index
                indices_found = True
        for line in f:
            split_line = line.rstrip().split()
            try:
                x = float(split_line[indices[axis]])
                seq = split_line[0]
                seq_len = len(seq)
            except IndexError:
                continue
            if not math.isnan(x) and seq_len in possible_lengths:
                if seq_csi_score_dict:
                    if seq in seq_csi_score_dict.keys():
                        score = seq_csi_score_dict[seq]
                        length_tuple_dict[seq_len].append((x, score))
                    else:
                        continue
                else:
                    length_tuple_dict[seq_len].append(x)
    return length_tuple_dict


def induction_plotter(data_file, csi_score_json, name, x_axis):
    """create axes variable and calls functions"""
    if csi_score_json:
        with open(csi_score_json, 'r') as f:
            seq_csi_score_dict = json.load(f)
    else:
        seq_csi_score_dict = False
    data_len_tuple_dict = get_data_by_column(x_axis, data_file, seq_csi_score_dict)
    if not name:
        name = '{0}_{1}_hist.pdf'.format(data_file.split('/')[-1].split('.')[0], x_axis)
    generate_plots(data_len_tuple_dict, x_axis, name, data_file.split('/')[-1].split('.')[0])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="script generates scatter plot of induction data from "
                                                 "genotype_fluorescence.py")
    parser.add_argument("-x", "--x_axis", default='free',
                        help="criterion to be plotted on x-axis (default: free)")
    parser.add_argument("-n", "--name",
                        help='name of output pdf (default: <data_file>.pdf')

    required = parser.add_argument_group('required arguments')
    required.add_argument("-d", "--data_file", required=True,
                          help="columns of data separated by white space. columns used to generate plot must be "
                               "named the same as the x and y axes")
    parser.add_argument('-s', '--csi_scores',
                        help='json file from csi scoring script')
    args = parser.parse_args()
    induction_plotter(args.data_file, args.csi_scores, args.name, args.x_axis)
