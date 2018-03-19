#!/usr/bin/env python
import argparse

import plotly

from induction_plotly import dataframe_from_file
import itertools


def parse_fasta_to_tag_seq_dict(fasta_file):
    tag_sequence_dict = {}
    with open(fasta_file, 'r') as f:
        for tag, sequence in itertools.izip_longest(f, f, fillvalue=None):
            tag = tag.rstrip().split('>', 1)[-1].split(' ')[0]
            sequence = sequence.rstrip()
            tag_sequence_dict[tag] = sequence
    return tag_sequence_dict


def parse_cdhit_to_cluster_tag_dict(cd_hit_cluster_file):
    cluster_dict = {}
    with open(cd_hit_cluster_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                cluster_number = int(line.rstrip().split()[-1])
                cluster_dict[cluster_number] = []
            else:
                pre_identifier = line.rstrip().split()[2]
                identifier_with_ellipsis = pre_identifier.split('>')[-1]
                identifier = identifier_with_ellipsis.split('...')[0]
                cluster_dict[cluster_number].append(identifier)
    return cluster_dict


def induction_cluster_plotly(data_file, cluster, fasta, name):
    dataframe = dataframe_from_file(data_file)
    tag_seq_dict = parse_fasta_to_tag_seq_dict(fasta)
    cluster_tag_list_dict = parse_cdhit_to_cluster_tag_dict(cluster)
    plotting_dict = {}
    base_pair_pairs = [''.join(x) for x in itertools.product('ACGT', 'ACGT') if x[0] is not x[1]]
    for pair in base_pair_pairs:
        plotting_dict[pair] = {'position': [], 'fold_fold': [], 'text': [],
                               'name': '{0} &#8594; {1}'.format(pair[0], pair[1])}

    for cluster_number, tag_list in cluster_tag_list_dict.iteritems():
        len_seq_dict = {16: [], 17: [], 18: [], 19: []}
        if len(tag_list) > 1:
            for tag in tag_list:
                seq = tag_seq_dict[tag]
                len_seq_dict[len(seq)].append(seq)

        for length in len_seq_dict.keys():
            for seq1, seq2 in itertools.combinations(len_seq_dict[length], 2):
                induction_vals = [dataframe.loc[seq1]['fold_induction'], dataframe.loc[seq2]['fold_induction']]
                repressed_vals = [dataframe.loc[seq1]['repressed'], dataframe.loc[seq2]['repressed']]
                induced_vals = [dataframe.loc[seq1]['induced'], dataframe.loc[seq2]['induced']]
                seq_ind_tuples = [(seq1, induction_vals[0]), (seq2, induction_vals[1])]
                seq_ind_tuples = sorted(seq_ind_tuples, key=lambda x: x[1], reverse=True)
                fold_fold_ind = seq_ind_tuples[0][1] / seq_ind_tuples[1][1]
                base_list = []
                changed_position = False
                for p, position in enumerate(seq_ind_tuples[1][0]):
                    if position != seq_ind_tuples[0][0][p]:
                        changed_position = p
                        low_base = position
                        high_base = seq_ind_tuples[0][0][p]
                        base_list.extend([low_base, high_base])
                if len(base_list) is not 2:
                    continue
                    # right now i do not want to handle double mutants, may add later
                    # raise Exception("there must be one and only one base changed in each compared pair of sequences\n"
                    #                 "this is not true for {0} and {1}".format(seq1, seq2))
                plotting_dict[''.join(base_list)]['position'].append(changed_position)
                plotting_dict[''.join(base_list)]['fold_fold'].append(fold_fold_ind)
                induction_vals = [round(x, 1) for x in induction_vals]
                repressed_vals = [round(x, 1) for x in repressed_vals]
                induced_vals = [round(x, 1) for x in induced_vals]

                plotting_dict[''.join(base_list)]['text'].append('{0} ({1}, {2}, {3})<br>{4} ({5}, {6}, {7})'.format(
                    seq1, repressed_vals[0], induced_vals[0], induction_vals[0],
                    seq2, repressed_vals[1], induced_vals[1], induction_vals[1]))

    protein = data_file.split('/')[-1].split('.')[0]
    fig = {
        'data': [
            {
                'x': plotting_dict[pair]['position'],
                'y': plotting_dict[pair]['fold_fold'],
                'text': plotting_dict[pair]['text'],
                'name': plotting_dict[pair]['name'],
                'mode': 'markers',
            } for pair in base_pair_pairs
            ],
        'layout': {
            'xaxis': {'title': 'position number'},
            'yaxis': {'title': 'fold fold induction'},
            'title': protein}
    }

    plotly.offline.plot(fig, filename='{0}'.format(name), auto_open=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="script generates scatter plot of induction data from "
                                                 "genotype_fluorescence.py")
    parser.add_argument("-n", "--name", default='cluster_comparison.html',
                        help='name of output html')

    required = parser.add_argument_group('required arguments')
    required.add_argument("-d", "--data_file", required=True,
                          help="columns of data separated by white space. columns used to generate plot must be "
                               "named the same as the x and y axes")
    parser.add_argument('-c', '--cluster',
                        help='a sorted cluster file from cd-hit')
    parser.add_argument('-f', '--fasta',
                        help='the fasta file containing sequences appearing in the sorted cluster file')
    args = parser.parse_args()
    induction_cluster_plotly(args.data_file, args.cluster, args.fasta, args.name)
