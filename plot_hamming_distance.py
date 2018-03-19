#!/usr/bin/env python
import argparse
import collections
import itertools

import plotly

from color_palettes import palette
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def generate_plots(data_tuple_dict, axes, name):
    """makes pdf scatter plot"""
    # set up colors and names
    almost_gray = '#808080'
    almost_black = '#262626'
    if len(data_tuple_dict) == 1:
        color_set = [almost_gray]
    else:
        color_set = palette[len(data_tuple_dict) - 1]

    with PdfPages(name) as pdf:
        for i, data_file in enumerate(data_tuple_dict.keys()):
            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            x, y = zip(*data_tuple_dict[data_file])
            ax1.scatter(x, y, c=color_set[i-1], marker='o', alpha=0.5, edgecolor=almost_black, linewidth=0.15,
                        label='{0}'.format(data_file.split('.')[0]))

            legend = ax1.legend(loc='best', scatterpoints=1, framealpha=1)
            rect = legend.get_frame()
            rect.set_linewidth(0.25)
            texts = legend.texts
            for t in texts:
                t.set_color(almost_black)
            # plt.xscale('log')
            ax1.set_xbound(lower=0.0)
            ax1.set_ybound(lower=0.0)
            plt.xlabel('log {0}'.format(axes[0]), fontsize=20)
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


def get_data_by_column(axes, data_file):
    indices = {}
    indices_found = False
    data_tuple_list = []
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
            except IndexError:
                continue
            if x != float('nan') and y != float('nan'):
                data_tuple_list.append((x, y))
    return data_tuple_list


def parse_data_from_genotype_fluorescence(data_file):
    length_seq_induction_dict = collections.defaultdict(dict)
    with open(data_file, 'r') as f:
        for line in f:
            if not line[0].isspace():
                split_line = line.rstrip().split()
                fold_induction = float(split_line[4])
                if not math.isnan(fold_induction):
                    sequence = split_line[0]
                    length_seq = len(sequence)
                    length_seq_induction_dict[length_seq][sequence] = fold_induction
    return length_seq_induction_dict


def compute_hamming_distance(string1, string2):
    """Return the Hamming distance between equal-length sequences"""
    if len(string1) != len(string2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(el1 != el2 for el1, el2 in zip(string1, string2))


def plot_score_files(data_files, name, x_axis, y_axis):
    """create axes variable and calls functions"""
    # axes = [x_axis, y_axis]
    # data_tuple_lists = collections.defaultdict(list)

    for data_file in data_files:
        plotting_dict = {
            16: {
                'hamming_distance': [],
                'fold_fold_induction': [],
                'seqs': []
            },
            17: {
                'hamming_distance': [],
                'fold_fold_induction': [],
                'seqs': []
            },
            18: {
                'hamming_distance': [],
                'fold_fold_induction': [],
                'seqs': []
            },
            19: {
                'hamming_distance': [],
                'fold_fold_induction': [],
                'seqs': []
            },
        }

        print "processing {0}".format(data_file)
        comparison_count = 0
        tuple_list = []
        length_seq_induction_dict = parse_data_from_genotype_fluorescence(data_file)
        for length in [16, 17, 18, 19]:
            seq_induction_dict = dict(length_seq_induction_dict[length].items())
            for seq1, seq2 in itertools.combinations(seq_induction_dict.keys(), 2):
                comparison_count += 1
                hamming_distance = compute_hamming_distance(seq1, seq2)
                induction_vals = [seq_induction_dict[seq1], seq_induction_dict[seq2]]
                fold_fold_ind = (max(induction_vals) - min(induction_vals)) / min(induction_vals)
                plotting_dict[length]['hamming_distance'].append(hamming_distance)
                plotting_dict[length]['fold_fold_induction'].append(fold_fold_ind)

                plotting_dict[length]['seqs'].append('{0}<br>{1}'.format(seq1, seq2))
                # data_tuple_lists[data_file].append((hamming_distance, fold_fold_ind))

        # sorted_list = sorted(data_tuple_lists[data_file], key=lambda x: (x[0], -x[1]))

        print 'comparisons made: {0}'.format(comparison_count)
        # data_tuple_lists[data_file] = tuple_list

        fig = {
            'data': [
                {
                    'x': plotting_dict[length]['hamming_distance'],
                    'y': plotting_dict[length]['fold_fold_induction'],
                    'text': plotting_dict[length]['seqs'],
                    'name': length,
                    'mode': 'markers',
                } for length in plotting_dict.keys()
                ],
            'layout': {
                'xaxis': {'title': 'hamming distance'},
                'yaxis': {'title': 'fold fold induction'},
                'title': data_file.split('.')[0]
            }
        }

        plotly.offline.plot(fig, filename='{0}_hamming_distance.html'.format(data_file.split('.')[0]), auto_open=False)

    # generate_plots(data_tuple_lists, axes, name)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="script generates plot of fold fold induction and hamming distance "
                                                 "for data in output files from genotype_fluorescence.py")
    parser.add_argument("-x", "--x_axis", default='hamming distance',
                        help="criterion to be plotted on x-axis (default: hamming distance)")
    parser.add_argument("-y", "--y_axis", default='fold fold induction',
                        help="criterion to be plotted on y-axis (default: fold fold induction)")
    parser.add_argument("-n", "--name", default='hamming_distance_plot.pdf',
                        help='name of output pdf (default: hamming_distance_plot.pdf')

    requiredO = parser.add_argument_group('required arguments')
    requiredO.add_argument("-d", "--data_files", nargs='*', required=True,
                           help="columns of data separated by white space")
    args = parser.parse_args()
    plot_score_files(args.data_files, args.name, args.x_axis, args.y_axis)
