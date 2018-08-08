#!/usr/bin/env python

# this file is part of the Raman Lab github repository: https://github.com/raman-lab/library_analysis
# author: nwhoppe
# created: 7/2/18

import argparse
import collections
import numpy as np
import pandas as pd
from color_palettes import palette


def plot_flow_hist(list_data_dict, medians_df, minimum, maximum, name, n_bins=100):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    color_list = palette[len(list_data_dict)]
    almost_black = '#262626'
    almost_gray = '#808080'
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    i = 0
    for file_name, data_list in list_data_dict.items():
        bins = np.logspace(np.log10(min(data_list)), np.log10(max(data_list)), n_bins)
        ax1.hist(data_list, bins, color=color_list[i], alpha=0.1)
        ax1.hist(data_list, bins, color=color_list[i], alpha=0.8, histtype='step')
        ax1.axvline(np.median(data_list), c=color_list[i])
        if not medians_df.empty:
            ax1.axvline(medians_df.loc[file_name]['gate_median'], c=color_list[i], linestyle='dashed')
        i += 1

    ax1.set_xscale('log')

    if minimum:
        ax1.set_xlim(left=minimum)
    if maximum:
        ax1.set_xlim(right=maximum)
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

    plt.savefig(name, dpi=300)
    plt.close()


def medians_from_csv(csv_list):
    medians = []
    for stat_csv in csv_list:
        temp_medians = []
        stats_df = pd.read_csv(stat_csv, index_col=False)
        median_list = stats_df[' Median'].tolist()
        for numeric_string in median_list:
            if not numeric_string.isspace():
                numeric_string.strip()
                temp_medians.append(float(numeric_string))
        temp_medians.pop(0)
        medians.extend(temp_medians)
    return medians


def main(input_csvs, column, medians_csv, name, minimum, maximum):

    column_data_dict = collections.OrderedDict()
    for csv in input_csvs:
        print csv
        df = pd.read_csv(csv, index_col=False)
        df = df[df[column] > 0]
        col_list = df[column].tolist()
        column_data_dict[csv] = col_list

    medians_df = pd.DataFrame()
    if medians_csv:
        medians_df = pd.read_csv(medians_csv, index_col=0)
        assert (len(medians_df) >= len(column_data_dict)), \
            'Error: each row of median_csv must correspond to an input csv'

    plot_flow_hist(column_data_dict, medians_df, minimum, maximum, name)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""plot histograms from input csvs. intended for use on csv files 
    from Sony FACS machine. may be amenable to other csvs""")
    parser.add_argument('-c', '--column', default='EGFP-A-Compensated',
                        help='header of column containing data to be plotted')
    parser.add_argument('-m', '--medians_csv',
                        help='csv containing gate-medians for each input csv. first column is input csv name. second '
                             'column is median for the gate. intended to compare gate-medians to reflow-medians')
    parser.add_argument('-n', '--name', default='flow_histogram.png', help='name of output png')
    parser.add_argument('-min', '--minimum', type=float, help='minimum x axis value')
    parser.add_argument('-max', '--maximum', type=float, help='maximum x axis value')
    required = parser.add_argument_group('required')
    required.add_argument('-i', '--input_csvs', nargs='*', required=True, help='input csv files')

    args = parser.parse_args()
    main(args.input_csvs, args.column, args.medians_csv, args.name, args.minimum, args.maximum)
