#!/usr/bin/env python
import argparse
import collections
import itertools
import json
from color_palettes import palette
import numpy as np
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import pandas as pd
import seaborn as sns


def generate_plots(data_len_tuple_dict, x_axis, y_axis, name):
    """makes pdf scatter plot"""
    # set up colors and names
    almost_gray = '#808080'
    almost_black = '#262626'
    color_set = palette[len(data_len_tuple_dict)]
    log_axes = ['repressed', 'induced', 'free']
    normal_axes = ['fold_induction', 'fold_repression', 'csi_score']

    with PdfPages(name) as pdf:
        fig = plt.figure()
        ax1 = fig.add_subplot(111)

        for i, length in enumerate(data_len_tuple_dict.keys()):
            if len(data_len_tuple_dict[length]) == 0:
                continue
            else:
                x, y = zip(*data_len_tuple_dict[length])
                ax1.scatter(x, y, c=color_set[i-1], marker='o', alpha=0.5, edgecolor=almost_black, linewidth=0.15,
                            label='{0}'.format(length))

        if y_axis is 'csi_score':
            ax1.set_ybound(lower=0)
        legend = ax1.legend(loc='best', scatterpoints=1, framealpha=1)
        rect = legend.get_frame()
        rect.set_linewidth(0.25)
        texts = legend.texts
        for t in texts:
            t.set_color(almost_black)

        plt.xlabel(x_axis, fontsize=20)
        plt.ylabel(y_axis, fontsize=20)

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
        pdf.savefig(fig)
        plt.close()

        # fig = plt.figure()
        # for i, length in enumerate(data_len_tuple_dict.keys()):
        #     # make individual length plots
        #     ax1 = fig.add_subplot(2, 2, i + 1)
        #     if len(data_len_tuple_dict[length]) == 0:
        #         continue
        #     elif len(data_len_tuple_dict[length][0]) == 3:
        #         x, y, z = zip(*data_len_tuple_dict[length])
        #         ax1.scatter(x, y, c=z, marker='o', alpha=0.5, edgecolor=almost_black, linewidth=0.15, cmap='brg',
        #                     label='{0}'.format(length))
        #     else:
        #         x, y = zip(*data_len_tuple_dict[length])
        #         ax1.scatter(x, y, c=color_set[i-1], marker='o', alpha=0.5, edgecolor=almost_black, linewidth=0.15,
        #                     label='{0}'.format(length))
        #
        #     legend = ax1.legend(loc='best', scatterpoints=1, framealpha=1)
        #     rect = legend.get_frame()
        #     rect.set_linewidth(0.25)
        #     texts = legend.texts
        #     for t in texts:
        #         t.set_color(almost_black)
        #     if axes[0] in log_axes:
        #         plt.xscale('log')
        #         ax1.set_xbound(lower=100.0)
        #     else:
        #         ax1.set_xbound(lower=0.0)
        #     if axes[1] in log_axes:
        #         plt.yscale('log')
        #         ax1.set_ybound(lower=100.0)
        #     else:
        #         ax1.set_ybound(lower=0.0)
        #     ax1.set_xbound(upper=x_max)
        #     plt.xlabel('{0}'.format(axes[0]), fontsize=10)
        #     plt.ylabel(axes[1], fontsize=10)
        #
        #     spines_to_remove = ['top', 'right']
        #     for spine in spines_to_remove:
        #         ax1.spines[spine].set_visible(False)
        #     ax1.xaxis.set_ticks_position('none')
        #     ax1.yaxis.set_ticks_position('none')
        #     spines_to_keep = ['bottom', 'left']
        #     for spine in spines_to_keep:
        #         ax1.spines[spine].set_linewidth(0.5)
        #         ax1.spines[spine].set_color(almost_black)
        #     ax1.xaxis.label.set_color(almost_black)
        #     ax1.yaxis.label.set_color(almost_black)
        # pdf.savefig(fig)
        # plt.close()


def get_data_by_column(axis, data_file):
    possible_lengths = [16, 17, 18, 19]
    indices = {}
    indices_found = False
    seq_val_dict = collections.OrderedDict()
    with open(data_file, 'r') as f:
        while not indices_found:
            line = f.readline().rstrip()
            if axis in line:
                split_line = line.split()
                index = split_line.index(axis) + 1
                indices[axis] = index
                indices_found = True
        for line in f:
            split_line = line.rstrip().split()
            try:
                val = float(split_line[indices[axis]])
                seq = split_line[0]
                seq_len = len(seq)
            except IndexError:
                continue
            if not math.isnan(val) and seq_len in possible_lengths:
                seq_val_dict[seq] = val

    return seq_val_dict


def palindrome_score(dna_sequence):
    max_score = 0
    max_second_half = ''
    max_rev_comp = ''
    for i in range(1, len(dna_sequence) - 1):
        first = dna_sequence[:i]
        second = dna_sequence[i:]
        first_seq = Seq(first, generic_dna)
        rev_comp_first = str(first_seq.reverse_complement())
        score = sum(el1 == el2 for el1, el2 in zip(second, rev_comp_first))
        if score > max_score:
            max_score = score
            max_second_half = second
            max_rev_comp = rev_comp_first
    return max_score, max_second_half, max_rev_comp


def violin_plot(xy_dataframe, basename):
    pd.options.mode.chained_assignment = None
    # xy_dataframe = xy_dataframe[xy_dataframe['pssm score'] >= 0]
    # min_x = int(min(xy_dataframe['pssm score'])) / 5
    # max_x = int(max(xy_dataframe['pssm score'])) + 5
    # bins = range(min_x, max_x, 5)
    # labels = []
    # for i in range(0, len(bins) - 1):
    #     label = '{0}-{1}'.format(bins[i], bins[i + 1])
    #     labels.append(label)
    # print pd.cut(xy_dataframe.loc[:, 'pssm score'], bins=bins)
    # xy_dataframe['bins'] = pd.cut(xy_dataframe['pssm score'], bins=bins, labels=labels)
    sns.set_style("whitegrid")
    ax = sns.violinplot(x="palindrome_length", y="fold_induction", data=xy_dataframe,
                        cut=0, scale='width', palette='muted')
    fig = ax.get_figure()
    fig.savefig('{0}_violin.pdf'.format(basename))


def df_scatter_plot(xy_df, base_name):
    almost_gray = '#808080'
    almost_black = '#262626'
    color_set = {
        'CmeR': '#ab62c0',
        'DesT': '#50ab6d',
        'NalC': '#c95573',
        'PmeR': '#929d3d',
        'SmeT': '#648ace',
        'TtgR': '#ca743e'
    }
    log_axes = ['repressed', 'induced', 'free']
    normal_axes = ['fold_induction', 'fold_repression', 'csi_score']

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    max_pal = max(xy_df['palindrome_length'])
    max_fi = max(xy_df['fold_induction'])
    x = np.asarray(xy_df['palindrome_length'])
    y = np.asarray(xy_df['fold_induction'])
    q90 = np.percentile(y, 90)
    x_below_q90 = x[np.where(y < q90)]
    y_below_q90 = y[np.where(y < q90)]
    x_above_q90 = x[np.where(y >= q90)]
    y_above_q90 = y[np.where(y >= q90)]
    # ax1.scatter(xy_df['palindrome_length'], xy_df['fold_induction'], c=color_set[base_name], marker='o', alpha=0.75,
                # edgecolor=almost_black, linewidth=0.15)
    ax1.scatter(x_below_q90, y_below_q90, c=almost_gray, marker='o', alpha=0.5, edgecolor=almost_black, linewidth=0.15)
    ax1.scatter(x_above_q90, y_above_q90, c=color_set[base_name], marker='o', alpha=0.75,
                edgecolor=almost_black, linewidth=0.15)

    plt.xlabel('Number of Palindromic Nucleotides', fontsize=15)
    plt.ylabel('Fold Induction', fontsize=15)
    ax1.text(max_pal, max_fi, s=base_name, fontsize=15, horizontalalignment='right')
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
    plt.savefig('{0}_pal_fi.pdf'.format(base_name))
    plt.close()


def induction_plotter(data_file, csi_score_json, spacer_length, name):
    """create axes variable and calls functions"""
    if csi_score_json:
        with open(csi_score_json, 'r') as f:
            seq_val_dict = json.load(f)
        y_axis = 'csi_score'
    else:
        seq_val_dict = get_data_by_column('fold_induction', data_file)
        y_axis = 'fold_induction'
    data_len_tuple_dict = {16: [], 17: [], 18: [], 19: []}
    if spacer_length:
        x_axis = 'spacer'
        for sequence, val in seq_val_dict.iteritems():
            ps, second_half, rev_comp = palindrome_score(sequence)
            if ps >= spacer_length:
                i = 0
                matching = False
                while not matching:
                    if second_half[i] == rev_comp[i]:
                        matching = True
                    else:
                        i += 1
                data_len_tuple_dict[len(sequence)].append((i, val))
    else:
        xy_df = pd.DataFrame(data=None, columns=['palindrome_length', 'fold_induction'])
        x_axis = 'palindrome'
        for sequence, val in seq_val_dict.iteritems():
            ps, second_half, rev_comp = palindrome_score(sequence)
            data_len_tuple_dict[len(sequence)].append((ps, val))
            xy_df.loc[sequence] = {'palindrome_length': ps, 'fold_induction': val}
    if not name:
        name = '{0}_{1}_{2}.pdf'.format(data_file.split('/')[-1].split('.')[0], x_axis, y_axis)
    # generate_plots(data_len_tuple_dict, x_axis, y_axis, name)
    df_scatter_plot(xy_df, data_file.split('/')[-1].split('.')[0])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="script generates scatter plot of induction data from "
                                                 "genotype_fluorescence.py")
    parser.add_argument("-n", "--name", help='name of output pdf')
    parser.add_argument("-l", "--spacer_length", type=int,
                        help="plot length between palindromic half sites where more than "
                             "m minimum bases are palindromic. ")

    required = parser.add_argument_group('required arguments')
    required.add_argument("-d", "--data_file", required=True, help="file containing column of fold_induction data")
    parser.add_argument('-s', '--csi_scores', help='json file from csi scoring script')
    args = parser.parse_args()
    induction_plotter(args.data_file, args.csi_scores, args.spacer_length, args.name)
