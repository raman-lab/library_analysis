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
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import numpy as np
import plotly
import plotly.figure_factory as ff
import seaborn as sns
import pandas as pd


def generate_plots(data_len_tuple_dict, axes, name):
    """makes pdf scatter plot"""
    # set up colors and names
    almost_gray = '#808080'
    almost_black = '#262626'
    color_set = palette[len(data_len_tuple_dict)]
    log_axes = ['repressed', 'induced', 'free', 'csi_score']
    normal_axes = ['fold_induction', 'fold_repression']

    with PdfPages(name) as pdf:
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        x_max = -1
        x_min = 1e6
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
                if max(x) > x_max:
                    x_max = max(x)
                if min(x) < x_min:
                    x_min = min(x)

        legend = ax1.legend(loc='best', scatterpoints=1, framealpha=1)
        rect = legend.get_frame()
        rect.set_linewidth(0.25)
        texts = legend.texts
        for t in texts:
            t.set_color(almost_black)

        if axes[0] in log_axes:
            plt.xscale('log')
        else:
            ax1.set_xbound(lower=0.0)
        if axes[1] in log_axes:
            plt.yscale('log')
            ax1.set_ybound(lower=0.0)
        else:
            ax1.set_ybound(lower=0.0)
        ax1.set_xbound(lower=x_min, upper=x_max)
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
            else:
                ax1.set_xbound(lower=0.0)
            if axes[1] in log_axes:
                plt.yscale('log')
                ax1.set_ybound(lower=0.0)
            else:
                ax1.set_ybound(lower=0.0)
            ax1.set_xbound(lower=x_min, upper=x_max)
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


def get_data_by_column(axes, data_file, seq_csi_score_dict):
    possible_lengths = [16, 17, 18, 19]
    indices = {}
    indices_found = False
    length_tuple_dict = collections.OrderedDict()
    for length in possible_lengths:
        length_tuple_dict[length] = []
    number_missed = 0
    total = 0
    with open(data_file, 'r') as f:
        while not indices_found:
            line = f.readline().rstrip()
            if axes[0] and axes[1] in line:
                split_line = line.split()
                x_index = split_line.index(axes[0]) + 1
                y_index = split_line.index(axes[1]) + 1
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
            except IndexError:
                continue
            if not math.isnan(x) and not math.isnan(y) and seq_len in possible_lengths:
                total += 1
                rev_comp = str(Seq(seq, generic_dna).reverse_complement())
                if seq_csi_score_dict:
                    if seq in seq_csi_score_dict.keys():
                        score = seq_csi_score_dict[seq]
                        length_tuple_dict[seq_len].append((score, y))
                    elif rev_comp in seq_csi_score_dict.keys():
                        score = seq_csi_score_dict[rev_comp]
                        length_tuple_dict[seq_len].append((score, y))
                    else:
                        number_missed += 1
                        continue
                else:
                    length_tuple_dict[seq_len].append((x, y))
    # print float(number_missed)
    # print float(total)
    # print float(number_missed) / total
    return length_tuple_dict


def contour_plot(x_list, y_list):
    almost_gray = '#808080'
    x_list = np.log10(x_list)
    df = pd.DataFrame({'csi score': x_list, 'fold induction': y_list})
    sns_plot = sns.jointplot(x='csi score', y='fold induction', data=df, kind='kde', n_levels=25)
    sns_plot.plot_joint(plt.scatter, c=almost_gray, s=5, marker='o')
    # sns_plot.ax_joint.collections[0].set_alpha(0)
    # sns_plot.set_axis_labels("$X$", "$Y$")
    sns_plot.savefig("sns_test.png")


def plot_broken_yaxis(x, y, name):
    almost_black = '#262626'
    almost_gray = '#808080'
    color_set = {
        'CmeR': '#ab62c0',
        'DesT': '#50ab6d',
        'NalC': '#c95573',
        'PmeR': '#929d3d',
        'SmeT': '#648ace',
        'TtgR': '#ca743e'
    }
    break_set = {
        'CmeR': [(22, 26), (0, 13)],
        'DesT': [(280.0, 370.0), (0.0, 55.0)],
        'NalC': [(240.0, 250.0), (0.0, 115.0)],
        'PmeR': [(250, 260), (0, 170)],
        'SmeT': [(25.0, 35.0), (0.0, 20.0)],
        'TtgR': [(175.0, 225.0), (0.0, 125.0)]
    }
    color = color_set[name]

    wild_type_fi = {
        'CmeR': [2847, 0.906743941],
        'DesT': [3242, 1.340530537],
        'NalC': [2007, 11.94419532],
        'PmeR': [1416.5, 6.417578539],
        'SmeT': [351.5, 1.119487909],
        'TtgR': [1444, 5.823407202]
    }
    wt_rep, wt_fi = wild_type_fi[name]

    # For bold font:
    plt.rcParams["font.weight"] = "bold"
    plt.rcParams["axes.labelweight"] = "bold"
    # note original tick label size = 12, fontsize = 15
    if name != 'SmeT':
        ylim, ylim2 = break_set[name]

        top_yratio = float(ylim[1] - ylim[0]) / (ylim2[1] - ylim2[0] + ylim[1] - ylim[0])
        bot_yratio = float(ylim2[1] - ylim2[0]) / (ylim2[1] - ylim2[0] + ylim[1] - ylim[0])
        gs = gridspec.GridSpec(2, 1, height_ratios=[top_yratio, bot_yratio])

        fig = plt.figure()
        ax0 = fig.add_subplot(111)
        ax1 = fig.add_subplot(gs[0])
        ax2 = fig.add_subplot(gs[1])
        ax1.scatter(x, y, marker='o', alpha=0.5, edgecolors=almost_black, linewidth=0.15, c=color)
        ax2.scatter(x, y, marker='o', alpha=0.5, edgecolors=almost_black, linewidths=0.15, c=color)
        ax2.scatter(wt_rep, wt_fi, marker='o', alpha=1, edgecolors=almost_black, linewidths=0.5, c=almost_black, s=100)
        ax1.set_ylim(ylim[0], ylim[1])
        ax2.set_ylim(ylim2[0], ylim2[1])

        plt.subplots_adjust(hspace=0.15)

        ax1.set_xscale('log')
        ax2.set_xscale('log')
        spines_to_remove = ['top', 'right']
        for ax in [ax0, ax1, ax2]:
            for spine in spines_to_remove:
                ax.spines[spine].set_visible(False)
            ax.xaxis.set_ticks_position('none')
            ax.yaxis.set_ticks_position('none')
            spines_to_keep = ['bottom', 'left']
            for spine in spines_to_keep:
                ax.spines[spine].set_linewidth(0.5)
                ax.spines[spine].set_color(almost_black)
            ax.xaxis.label.set_color(almost_black)
            ax.yaxis.label.set_color(almost_black)
        ax1.spines['bottom'].set_visible(False)

        ax2.tick_params(axis='x', which='both', bottom='on', top='off')
        # ax2.tick_params(axis='y', which='both', right='off', left='off')
        ax1.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
        ax0.spines['bottom'].set_visible(False)
        ax0.spines['left'].set_visible(False)
        ax0.tick_params(top=0, bottom=0, right=0, left=0, labelcolor='w', labelsize=12)
        # ax1.tick_params(axis='y', which='both', right='off', left='off')

        # diagonal cut-out lines
        # d = 0.015
        # kwargs = dict(transform=ax1.transAxes, color=almost_black, clip_on=False)
        # ax1.plot((-d, +d), (-d, +d), **kwargs)
        # # ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)
        #
        # kwargs.update(transform=ax2.transAxes)
        # ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)
        # ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)

        # ax2.yaxis.set_label_coords(0.05, 0.5, transform=fig.transFigure)
        # kwargs = dict(color='k', clip_on=False)
        # xlim = ax1.get_xlim()
        # dx = .02 * (xlim[1] - xlim[0])
        # dy = .01 * (ylim[1] - ylim[0]) / top_yratio
        # ax1.plot((xlim[0] - dx, xlim[0] + dx), (ylim[0] - dy, ylim[0] + dy), **kwargs)
        # ax1.plot((xlim[1] - dx, xlim[1] + dx), (ylim[0] - dy, ylim[0] + dy), **kwargs)
        # dy = .01 * (ylim2[1] - ylim2[0]) / bot_yratio
        # ax2.plot((xlim[0] - dx, xlim[0] + dx), (ylim2[1] - dy, ylim2[1] + dy), **kwargs)
        # ax2.plot((xlim[1] - dx, xlim[1] + dx), (ylim2[1] - dy, ylim2[1] + dy), **kwargs)
        # ax1.set_xlim(xlim)
        # ax2.set_xlim(xlim)

        ax1.set_xlim(xmin=0.9*min(x), xmax=1.1*max(x))
        ax2.set_xlim(xmin=0.9*min(x), xmax=1.1*max(x))
        ax1.tick_params(labelsize=13, labelbottom=False)
        ax2.tick_params(labelsize=13)
        ax0.set_xlabel('Repression', fontsize=16)
        ax0.set_ylabel('Fold Induction', fontsize=16)
        ax1.text(max(x), ylim[1], s=name, fontsize=16, horizontalalignment='right')
        plt.savefig('{0}.pdf'.format(name))

        # formatting Vatsan likes
        # ax1.tick_params(labelsize=10, labelbottom=False)
        # ax2.tick_params(labelsize=10)
        # ax0.set_xlabel('Repression', fontsize=15, fontname='Times New Roman')
        # ax0.set_ylabel('Fold Induction', fontsize=15, fontname='Times New Roman')
        # ax2.xaxis.set_minor_formatter(matplotlib.ticker.LogFormatter(minor_thresholds=(2, 0.4)))
        # ax2.xaxis.set_major_formatter(matplotlib.ticker.LogFormatter())
        # plt.savefig('{0}.png'.format(name), dpi=1000)

        # pts = np.array([0.015, 0.166, 0.133, 0.159, 0.041, 0.024, 0.195,
        #                 0.039, 0.161, 0.018, 0.143, 0.056, 0.125, 0.096, 0.094, 0.051,
        #                 0.043, 0.021, 0.138, 0.075, 0.109, 0.195, 0.05, 0.074, 0.079,
        #                 0.155, 0.02, 0.01, 0.061, 0.008])
        # pts[[3, 14]] += .8
        #
        # ylim = [0.82, 1.0]
        # ylim2 = [0.0, 0.32]
        # ylimratio = (ylim[1] - ylim[0]) / (ylim2[1] - ylim2[0] + ylim[1] - ylim[0])
        # ylim2ratio = (ylim2[1] - ylim2[0]) / (ylim2[1] - ylim2[0] + ylim[1] - ylim[0])
        # gs = gridspec.GridSpec(2, 1, height_ratios=[ylimratio, ylim2ratio])
        # fig = plt.figure()
        # ax = fig.add_subplot(gs[0])
        # ax2 = fig.add_subplot(gs[1])
        # ax.plot(pts)
        # ax2.plot(pts)
        # ax.set_ylim(ylim)
        # ax2.set_ylim(ylim2)
        # plt.subplots_adjust(hspace=0.03)
        #
        # ax.spines['bottom'].set_visible(False)
        # ax2.spines['top'].set_visible(False)
        # ax.xaxis.tick_top()
        # ax.tick_params(labeltop='off')
        # ax2.xaxis.tick_bottom()
        #
        # ax2.set_xlabel('xlabel')
        # ax2.set_ylabel('ylabel')
        # ax2.yaxis.set_label_coords(0.05, 0.5, transform=fig.transFigure)
        #
        # kwargs = dict(color='k', clip_on=False)
        # xlim = ax.get_xlim()
        # dx = .02 * (xlim[1] - xlim[0])
        # dy = .01 * (ylim[1] - ylim[0]) / ylimratio
        # ax.plot((xlim[0] - dx, xlim[0] + dx), (ylim[0] - dy, ylim[0] + dy), **kwargs)
        # ax.plot((xlim[1] - dx, xlim[1] + dx), (ylim[0] - dy, ylim[0] + dy), **kwargs)
        # dy = .01 * (ylim2[1] - ylim2[0]) / ylim2ratio
        # ax2.plot((xlim[0] - dx, xlim[0] + dx), (ylim2[1] - dy, ylim2[1] + dy), **kwargs)
        # ax2.plot((xlim[1] - dx, xlim[1] + dx), (ylim2[1] - dy, ylim2[1] + dy), **kwargs)
        # ax.set_xlim(xlim)
        # ax2.set_xlim(xlim)
        #
        # plt.savefig('broken_axis-mod.png')
    else:
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.scatter(x, y, marker='o', alpha=0.5, edgecolors=almost_black, linewidth=0.15, c=color)
        ax1.scatter(wt_rep, wt_fi, marker='o', alpha=1, edgecolors=almost_black, linewidths=0.5, c=almost_black, s=100)
        ax1.set_xscale('log')
        spines_to_remove = ['top', 'right']
        for spine in spines_to_remove:
            ax1.spines[spine].set_visible(False)
        # ax.xaxis.set_ticks_position('none')
        ax1.yaxis.set_ticks_position('none')
        spines_to_keep = ['bottom', 'left']
        for spine in spines_to_keep:
            ax1.spines[spine].set_linewidth(0.5)
            ax1.spines[spine].set_color(almost_black)
        ax1.xaxis.label.set_color(almost_black)
        ax1.yaxis.label.set_color(almost_black)
        ax1.tick_params(axis='x', which='both', bottom='on', top='off')
        ax1.tick_params(labelsize=13)
        ax1.set_xlim(xmin=0.9*min(x), xmax=1.1*max(x))
        plt.xlabel('Repression', fontsize=16)
        plt.ylabel('Fold Induction', fontsize=16)
        ax1.text(max(x), max(y), s=name, fontsize=16, horizontalalignment='right')
        plt.savefig('{0}.pdf'.format(name))


def induction_plotter(data_file, csi_score_json, name, x_axis, y_axis):
    """create axes variable and calls functions"""
    axes = [x_axis, y_axis]
    if csi_score_json:
        with open(csi_score_json, 'r') as f:
            seq_csi_score_dict = json.load(f)
    else:
        seq_csi_score_dict = False
    data_len_tuple_dict = get_data_by_column(axes, data_file, seq_csi_score_dict)

    if seq_csi_score_dict:
        x_axis = 'csi_score'
        axes[0] = x_axis
    if not name:
        name = '{0}_{1}_{2}.pdf'.format(data_file.split('/')[-1].split('.')[0], x_axis, y_axis)
    # generate_plots(data_len_tuple_dict, axes, name)
    x_all = []
    y_all = []
    for length in data_len_tuple_dict.keys():
        x, y = zip(*data_len_tuple_dict[length])
        x_all.extend(x)
        y_all.extend(y)
    # contour_plot(x_all, y_all)

    plot_broken_yaxis(x_all, y_all, data_file.split('/')[-1].split('.')[0])

    # colorscale = ['#7A4579', '#D56073', 'rgb(236,158,105)', (1, 1, 0.2), (0.98,0.98,0.98)]
    # fig = ff.create_2d_density(
    #     np.log10(x_all), y_all, colorscale=colorscale,
    #     hist_color='rgb(255, 237, 222)', point_size='3', ncontours=5
    # )
    # plotly.offline.plot(fig, filename='test5.html', auto_open=False)


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
    required.add_argument("-d", "--data_file", required=True,
                           help="columns of data separated by white space. columns used to generate plot must be "
                                "named the same as the x and y axes")
    parser.add_argument('-s', '--csi_scores',
                        help='json file from csi scoring script')
    args = parser.parse_args()
    induction_plotter(args.data_file, args.csi_scores, args.name, args.x_axis, args.y_axis)
