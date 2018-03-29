#!/usr/bin/env python
import argparse
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from plotly_csv_scatter import dataframe_from_file


def gc_content_of_seq(string):
    g_count = string.count('G')
    c_count = string.count('C')
    return float(g_count + c_count) / len(string)


def scatter_plot(x, y, name, y_axis):
    almost_black = '#262626'
    almost_gray = '#808080'
    log_axes = ['repressed', 'induced', 'free']
    normal_axes = ['fold_induction', 'fold_repression']

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.scatter(x, y, marker='o', alpha=0.5, edgecolors=almost_black, linewidth=0.2, c=almost_gray)

    plt.xlabel('GC Content', fontsize=20)
    plt.ylabel(y_axis, fontsize=20)

    if y_axis in log_axes:
        plt.yscale('log')
        ax1.set_ybound(upper=1e5)
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
    plt.savefig(name)
    plt.close()


def induction_gc_plot(input_file, y_axis):
    allowable_lengths = [16, 17, 18, 19]
    input_df = dataframe_from_file(input_file)
    input_df.dropna(subset=[y_axis], how='any', inplace=True)

    x = []
    y = []
    for sequence in input_df.index:
        if sequence and len(sequence) in allowable_lengths:
            y_val = input_df.loc[sequence][y_axis]
            if y_val:
                gc_content = gc_content_of_seq(sequence)
                x.append(gc_content)
                y.append(y_val)
    name = '{0}_gc_{1}.pdf'.format(input_file.split('/')[-1].split('.')[0], y_axis)
    scatter_plot(x, y, name, y_axis)


if __name__ == '__main__':
    parser = argparse.ArgumentParser("""script to plot gc content vs fold induction for promoter library data""")
    parser.add_argument("-y", "--y_axis", default='fold_induction',
                        help="criterion to be plotted on y-axis (default: fold_induction)")
    required = parser.add_argument_group('required')
    required.add_argument('-i', '--input_file', required=True, help='file from genotype_fluorescence.py')

    args = parser.parse_args()
    induction_gc_plot(args.input_file, args.y_axis)
