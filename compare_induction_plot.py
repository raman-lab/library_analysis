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


def generate_plots(data_len_tuple_dict, axes, name):
    """makes pdf scatter plot"""
    # set up colors and names
    almost_gray = '#808080'
    almost_black = '#262626'
    color_set = palette[len(data_len_tuple_dict)]
    log_axes = ['repressed', 'induced', 'free', 'Repressed']
    normal_axes = ['fold_induction', 'fold_repression', 'csi_score', 'FoldChangeInducedToRepressed']

    with PdfPages(name) as pdf:
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        for i, length in enumerate(data_len_tuple_dict.keys()):
            if len(data_len_tuple_dict[length]) == 0:
                continue
            # elif len(data_len_tuple_dict[length][0]) == 3:
            #     x, y, z = zip(*data_len_tuple_dict[length])
            #     sc = plt.scatter(x, y, c=z, marker='o', alpha=0.5, edgecolor=almost_black, linewidth=0.15,
            #                      cmap='brg',
            #                 label='{0}'.format(length))
            else:
                x, y = zip(*data_len_tuple_dict[length])
                ax1.scatter(x, y, c=color_set[i-1], marker='o', alpha=0.5, edgecolor=almost_black, linewidth=0.15,
                            label='{0}'.format(length))

        legend = ax1.legend(loc='best', scatterpoints=1, framealpha=1)
        rect = legend.get_frame()
        rect.set_linewidth(0.25)
        texts = legend.texts
        for t in texts:
            t.set_color(almost_black)

        if axes[0] in log_axes:
            plt.xscale('log')
            ax1.set_xbound(lower=100.0)
        else:
            ax1.set_xbound(lower=0.0)
        if axes[1] in log_axes:
            plt.yscale('log')
            ax1.set_ybound(lower=100.0)
        else:
            ax1.set_ybound(lower=0.0)
        plt.xlabel('{0}'.format(axes[0]), fontsize=20)
        plt.ylabel(axes[1], fontsize=20)

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

        fig = plt.figure()
        for i, length in enumerate(data_len_tuple_dict.keys()):
            # make individual length plots
            ax1 = fig.add_subplot(2, 2, i + 1)
            if len(data_len_tuple_dict[length]) == 0:
                continue
            elif len(data_len_tuple_dict[length][0]) == 3:
                x, y, z = zip(*data_len_tuple_dict[length])
                ax1.scatter(x, y, c=z, marker='o', alpha=0.5, edgecolor=almost_black, linewidth=0.15, cmap='brg',
                            label='{0}'.format(length))
            else:
                x, y = zip(*data_len_tuple_dict[length])
                ax1.scatter(x, y, c=color_set[i-1], marker='o', alpha=0.5, edgecolor=almost_black, linewidth=0.15,
                            label='{0}'.format(length))

            legend = ax1.legend(loc='best', scatterpoints=1, framealpha=1)
            rect = legend.get_frame()
            rect.set_linewidth(0.25)
            texts = legend.texts
            for t in texts:
                t.set_color(almost_black)
            if axes[0] in log_axes:
                plt.xscale('log')
                ax1.set_xbound(lower=100.0)
            else:
                ax1.set_xbound(lower=0.0)
            if axes[1] in log_axes:
                plt.yscale('log')
                ax1.set_ybound(lower=100.0)
            else:
                ax1.set_ybound(lower=0.0)
            plt.xlabel('{0}'.format(axes[0]), fontsize=10)
            plt.ylabel(axes[1], fontsize=10)

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


def get_data_by_column(axes, txt_files, seq_list):
    possible_lengths = [16, 17, 18, 19]

    length_tuple_dict = collections.OrderedDict()
    for length in possible_lengths:
        length_tuple_dict[length] = []
    total_matched = 0
    total_seqs = 0
    seqs_with_vals = 0
    seqs_with_vals_matched = 0
    for data_file in txt_files:
        indices = {}
        indices_found = False
        with open(data_file, 'r') as f:
            while not indices_found:
                line = f.readline().rstrip()
                if axes[0] and axes[1] in line:
                    split_line = line.split()
                    x_index = split_line.index(axes[0])
                    y_index = split_line.index(axes[1])
                    indices[axes[0]] = x_index
                    indices[axes[1]] = y_index
                    indices_found = True
            for line in f:
                split_line = line.rstrip().split()
                try:
                    x = float(split_line[indices[axes[0]]])
                    y = float(split_line[indices[axes[1]]])
                    seq = split_line[0]
                    seq_len = len(seq)
                    total_seqs += 1
                    if seq in seq_list:
                        total_matched += 1
                except IndexError:
                    continue
                if x == 0.0 or y == 0.0:
                    continue
                if not math.isnan(x) and not math.isnan(y) and seq_len in possible_lengths:
                    seqs_with_vals += 1
                    length_tuple_dict[seq_len].append((x, y))
                    if seq in seq_list:
                        seqs_with_vals_matched += 1

    print "matched / total: {0} / {1} = {2}".format(float(total_matched), total_seqs, float(total_matched) / total_seqs)
    print "with values: {0} / {1} = {2}".format(
        float(seqs_with_vals_matched), seqs_with_vals, float(seqs_with_vals_matched) / seqs_with_vals)
    return length_tuple_dict


def get_seqs_from_data_file(data_file):
    seq_list = []
    with open(data_file, 'r') as f:
        f.next()
        f.next()
        for line in f:
            seq_list.append(line.split()[0])

    return seq_list


def induction_plotter(txt_files, data_file, name, x_axis, y_axis):
    """create axes variable and calls functions"""
    axes = [x_axis, y_axis]
    seq_list = get_seqs_from_data_file(data_file)
    data_len_tuple_dict = get_data_by_column(axes, txt_files, seq_list)

    if not name:
        name = '{0}_{1}_{2}.pdf'.format(data_file.split('/')[-1].split('.')[0], x_axis, y_axis)
    generate_plots(data_len_tuple_dict, axes, name)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="script generates scatter plot of induction data from "
                                                 "genotype_fluorescence.py")
    parser.add_argument("-x", "--x_axis", default='repressed',
                        help="criterion to be plotted on x-axis (default: repressed)")
    parser.add_argument("-y", "--y_axis", default='fold_induction',
                        help="criterion to be plotted on y-axis (default: fold_induction)")
    parser.add_argument("-n", "--name",
                        help='name of output pdf (default: <data_file>.pdf')

    required = parser.add_argument_group('required arguments')
    required.add_argument("-t", "--txt_file", required=True, nargs='*',
                          help="columns of data separated by white space. columns used to generate plot must be "
                               "named the same as the x and y axes")
    parser.add_argument('-d', '--data_file',
                        help='')
    args = parser.parse_args()
    induction_plotter(args.txt_file, args.data_file, args.name, args.x_axis, args.y_axis)
