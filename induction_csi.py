#!/usr/bin/env python
import argparse
import collections
import itertools
import json
from color_palettes import palette
import math
import matplotlib

from induction_plotly import dataframe_from_file

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import numpy as np
import plotly
import plotly.figure_factory as ff
import seaborn as sns
import pandas as pd


def contour_plot(xy_dataframe, name_prefix):
    almost_gray = '#808080'
    sns_plot = sns.jointplot(x='csi score', y='fold induction', data=xy_dataframe, kind='kde', n_levels=25)
    sns_plot.plot_joint(plt.scatter, c=almost_gray, s=5, marker='o')
    sns_plot.savefig("{0}.png".format(name_prefix))


# def median_from_binned_data(data_array, bins):
#     bin_indices = np.digitize(data_array, bins)
#     for index, bin_value in enumerate(bins):
#
#     return median_array, mad_array_positive, mad_array_negative


def marginal_plot(xy_dataframe, name_prefix):
    almost_gray = '#808080'
    almost_black = '#262626'
    x = xy_dataframe.ix[:, 0]
    y = xy_dataframe.ix[:, 1]

    nullfmt = matplotlib.ticker.NullFormatter()
    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left + width + 0.02

    rect_scatter = [left, bottom, width, height]
    rect_marginal_x = [left, bottom_h, width, 0.2]
    rect_marginal_y = [left_h, bottom, 0.2, height]

    # start with a rectangular Figure
    plt.figure(1, figsize=(8, 8))

    axScatter = plt.axes(rect_scatter)
    axHistx = plt.axes(rect_marginal_x)
    axHisty = plt.axes(rect_marginal_y)

    # no labels
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)

    # the scatter plot:
    axScatter.scatter(x, y, c=almost_gray, marker='o', alpha=0.5, edgecolor=almost_black)

    # now determine nice limits by hand:
    # binwidth = 0.25
    # xymax = np.max([np.max(np.fabs(x)), np.max(np.fabs(y))])
    # lim = (int(xymax / binwidth) + 1) * binwidth

    axScatter.set_xlim((min(x), max(x)))
    axScatter.set_ylim((min(y), max(y)))

    x_bins = np.linspace(min(x), max(x), 20)
    y_bins = np.linspace(min(y), max(y), 20)

    x_medians, x_mad_max, y_mad_min = median_from_binned_data(x, x_bins)

    print x_bins
    print y_bins
    axHistx.plot(x_bins, x_medians, color=almost_black, alpha=0.5)
    axHistx.fill_between(x_bins, x_mad_max, x_mad_min, color=almost_gray, alpha=0.3)

    axHisty.plot(y_bins, y_medians, color=almost_black, alpha=0.5, orientation='horizontal')
    axHisty.fill_between(ybins, y_mad_max, y_mad_min, color=almost_gray, alpha=0.3)

    axHistx.set_xlim(axScatter.get_xlim())
    axHisty.set_ylim(axScatter.get_ylim())

    plt.savefig('{0}_marginal.png'.format(name_prefix))


def mad_plot(xy_dataframe, name_prefix):
    almost_gray = '#808080'
    almost_black = '#262626'
    percentile = 80
    x = xy_dataframe.ix[:, 0]
    y = xy_dataframe.ix[:, 1]
    x_bins = np.linspace(min(x) - 1, max(x), 21)
    bin_indicies = np.digitize(x, x_bins)
    means = []
    stdevs_pos = []
    stdevs_neg = []
    sums = []

    for i in range(1, len(x_bins)):
        y_bin = y[bin_indicies == i]
        q = np.percentile(y_bin, percentile)

        stdev = np.std(y_bin)
        means.append(q)
        stdevs_pos.append(mu + stdev)
        stdevs_neg.append(mu - stdev)
        if len(y_bin) != 0:
            summed = np.sum(y_bin) / float(len(y_bin))
            sums.append(summed)
        else:
            sums.append('nan')

    plt.plot(x_bins[1:], means, color=almost_black, alpha=0.75)
    plt.fill_between(x_bins[1:], stdevs_pos, stdevs_neg, color=almost_gray, alpha=0.3)

    plt.savefig('{0}_stdev.png'.format(name_prefix))


def induction_csi_plotter(data_file, csi_json_files, name_prefix, filter_bool):
    # count dict
    count_var_dict = {'all_df_seqs': 0, 'induction_df_seqs': 0, 'post_filter_df_seqs': 0, 'csi_seqs': 0,
                      'csi_rev_comp_seqs': 0}

    if not name_prefix:
        name_prefix = data_file.split('/')[-1].split('.')[0]

    allowable_lengths = [16, 17, 18, 19]
    input_df = dataframe_from_file(data_file)
    count_var_dict['all_df_seqs'] = len(input_df.index)
    input_df.dropna(subset=['fold_induction'], how='any', inplace=True)
    count_var_dict['induction_df_seqs'] = len(input_df.index)
    csi_dict_list = []
    for csi_file in csi_json_files:
        with open(csi_file, 'r') as f:
            seq_csi_score_dict = json.load(f)
            csi_dict_list.append(seq_csi_score_dict)

    if filter_bool:
        input_df = input_df[
            (input_df.repressed < input_df.repressed.quantile(0.8)) & (input_df.fold_induction > 2)
        ]
        count_var_dict['post_filter_df_seqs'] = len(input_df.index)

    output_df = pd.DataFrame(columns=['csi score', 'fold induction'])
    for sequence in input_df.index:
        if sequence is not None and len(sequence) in allowable_lengths:
            csi_value = False
            rev_comp = str(Seq(sequence, generic_dna).reverse_complement())
            for seq_csi_score_dict in csi_dict_list:
                if sequence in seq_csi_score_dict.keys():
                    csi_value = seq_csi_score_dict[sequence]
                    count_var_dict['csi_seqs'] += 1
                    break
                elif rev_comp in seq_csi_score_dict.keys():
                    csi_value = seq_csi_score_dict[rev_comp]
                    count_var_dict['csi_rev_comp_seqs'] += 1
                    break

            if csi_value:
                fold_ind = input_df.loc[sequence]['fold_induction']
                output_df.loc[sequence] = {'csi score': np.log10(csi_value), 'fold induction': fold_ind}

    # contour_plot(output_df, name_prefix)
    #
    # colorscale = ['#7A4579', '#D56073', 'rgb(236,158,105)', (1, 1, 0.2), (0.98,0.98,0.98)]
    # fig = ff.create_2d_density(
    #     output_df['csi score'], output_df['fold induction'], colorscale=colorscale,
    #     hist_color='rgb(255, 237, 222)', point_size='3', ncontours=10
    # )
    # plotly.offline.plot(fig, filename='{0}.html'.format(name_prefix), auto_open=False)

    # marginal_plot(output_df, name_prefix)

    mad_plot(output_df, name_prefix)
    print count_var_dict


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="script generates scatter plot of induction data from "
                                                 "genotype_fluorescence.py")
    parser.add_argument("-n", "--name", help='prefix for output files')
    parser.add_argument('-f', '--filter', action='store_true',
                        help='enforce filter on induction sequences used. the filter requires sequences have a fold '
                             'induction greater than 2 and a repressed fluorescence value less than the 80th '
                             'percentile. in the future, the values used by the filter may be flexible')

    required = parser.add_argument_group('required arguments')
    required.add_argument("-d", "--data_file", required=True)
    parser.add_argument('-c', '--csi_scores', nargs='*', required=True,
                        help='json file(s) from csi scoring script')
    args = parser.parse_args()
    induction_csi_plotter(args.data_file, args.csi_scores, args.name, args.filter)
