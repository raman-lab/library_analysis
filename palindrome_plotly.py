#!/usr/bin/env python
import argparse
import collections
import json
import math
import plotly
import sys
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import plotly.plotly as py
import plotly.graph_objs as go
import pandas as pd


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
    max_first_half = ''
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
            max_first_half = first
            max_second_half = second
            max_rev_comp = rev_comp_first
    return max_score, max_first_half, max_second_half, max_rev_comp


def induction_plotter(data_file, csi_score_json, spacer_length, name):
    """create axes variable and calls functions"""
    # if csi_score_json:
    pssm_df = pd.read_json(csi_score_json)
    seq_val_dict = get_data_by_column('fold_induction', data_file)
    plotting_dict = {
        16: {
            'x': [],
            'y': [],
            'z': [],
            'text': []
        },
        17: {
            'x': [],
            'y': [],
            'z': [],
            'text': []
        },
        18: {
            'x': [],
            'y': [],
            'z': [],
            'text': []
        },
        19: {
            'x': [],
            'y': [],
            'z': [],
            'text': []
        }
    }
    sys.stdout.write('sequence\tfold_induction\tpssm_score\tpalindrome_length\tgap\n')
    for sequence, val in seq_val_dict.iteritems():
        seq_len = len(sequence)
        ps, first_half, second_half, rev_comp = palindrome_score(sequence)
        if ps >= spacer_length:
            i = 0
            matching = False
            while not matching:
                if second_half[i] == rev_comp[i]:
                    matching = True
                else:
                    i += 1
            plotting_dict[seq_len]['x'].append(ps)
            plotting_dict[seq_len]['y'].append(i)
            plotting_dict[seq_len]['z'].append(val)
            plotting_dict[seq_len]['text'].append('{0}|{1}'.format(first_half, second_half))
            try:
                pssm_val = pssm_df.loc[sequence]['pssm score']
            except KeyError:
                continue
            sys.stdout.write('{0}|{1}\t'.format(first_half, second_half))
            sys.stdout.write('{0}\t'.format(val))
            sys.stdout.write('{0}\t'.format(pssm_val))
            sys.stdout.write('{0}\t'.format(ps))
            sys.stdout.write('{0}\n'.format(i))

    if not name:
        name = '{0}'.format(data_file.split('/')[-1].split('.')[0])
    data = []
    for length in plotting_dict.keys():
        trace = go.Scatter3d(
            x=plotting_dict[length]['x'],
            y=plotting_dict[length]['y'],
            # z=plotting_dict[length]['z'],
            text=plotting_dict[length]['text'],
            name=length,
            mode='markers',
            opacity=0.8
        )
        data.append(trace)

    layout = go.Layout(
        title=name,
        scene=dict(
            xaxis=dict(
                title='number of palindromic bases'
            ),
            yaxis=dict(
                title='gap between palindromic bases'
            ),
            zaxis=dict(
                title='fold induction'
            ),
        ),
    )
    # fig = {
    #     'data': [
    #         {
    #             'x': plotting_dict[length]['x'],
    #             'y': plotting_dict[length]['y'],
    #             'z': plotting_dict[length]['z'],
    #             'text': plotting_dict[length]['text'],
    #             'name': length,
    #             'mode': 'markers',
    #         } for length in plotting_dict.keys()
    #         ],
    #     'layout': {
    #         'xaxis': {'title': 'palindrome score'},
    #         'yaxis': {'title': 'palindrome gap'},
    #         'zaxis': {'title': 'fold induction'},
    #         'title': data_file.split('/')[-1].split('.')[0]}
    # }
    fig = go.Figure(data=data, layout=layout)
    # plotly.offline.plot(fig, filename='{0}.html'.format(name), auto_open=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="script generates scatter plot of induction data from "
                                                 "genotype_fluorescence.py")
    parser.add_argument("-n", "--name", help='name of output pdf')
    parser.add_argument("-l", "--spacer_length", type=int,
                        help="plot length between palindromic half sites where more than "
                             "m minimum bases are palindromic. ")

    required = parser.add_argument_group('required arguments')
    required.add_argument("-d", "--data_file", required=True, help="file containing column of fold_induction data")
    parser.add_argument('-s', '--pssm_scores', required=True, help='json file from pssm scoring script')
    args = parser.parse_args()
    induction_plotter(args.data_file, args.pssm_scores, args.spacer_length, args.name)
