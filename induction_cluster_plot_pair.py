#!/usr/bin/env python
import argparse

import plotly
import sys

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


def induction_cluster_plotly(data_file1, data_file2, cluster, fasta1, fasta2, name, exact):
    dataframe1 = dataframe_from_file(data_file1)
    dataframe2 = dataframe_from_file(data_file2)
    tag_seq_dict1 = parse_fasta_to_tag_seq_dict(fasta1)
    tag_seq_dict2 = parse_fasta_to_tag_seq_dict(fasta2)
    cluster_tag_list_dict = parse_cdhit_to_cluster_tag_dict(cluster)
    plotting_dict = {
        'cluster': [],
        'ratio': [],
        'text': []
    }
    sys.stdout.write('seq\trepressed\tinduced\tfi\tffi\n')
    for cluster_number, tag_list in cluster_tag_list_dict.iteritems():
        if len(tag_list) > 1:
            # len_seq_dict = {16: [], 17: [], 18: [], 19: []}
            # for tag in tag_list:
            #     seq = tag_seq_dict[tag]
            #     len_seq_dict[len(seq)].append(seq)
            # for length in len_seq_dict.keys():
            for tag1, tag2 in itertools.combinations(tag_list, 2):
                pro1 = tag1.rsplit('_', 1)[0]
                if pro1 in tag2:
                    continue

                if tag1 in tag_seq_dict1.keys():
                    seq1 = tag_seq_dict1[tag1]
                elif tag1 in tag_seq_dict2.keys():
                    seq2 = tag_seq_dict2[tag1]
                else:
                    raise Exception('tag {0} not found in either fasta file'.format(tag1))
                if tag2 in tag_seq_dict1.keys():
                    seq1 = tag_seq_dict1[tag2]
                elif tag2 in tag_seq_dict2.keys():
                    seq2 = tag_seq_dict2[tag2]
                else:
                    raise Exception('tag {0} not found in either fasta file'.format(tag2))

                if exact:
                    if seq1 != seq2:
                        continue

                # for seq1, seq2 in itertools.combinations(len_seq_dict[length], 2):
                induction_vals = [dataframe1.loc[seq1]['fold_induction'], dataframe2.loc[seq2]['fold_induction']]
                repressed_vals = [dataframe1.loc[seq1]['repressed'], dataframe2.loc[seq2]['repressed']]
                induced_vals = [dataframe1.loc[seq1]['induced'], dataframe2.loc[seq2]['induced']]
                fold_fold_ind = induction_vals[0] / induction_vals[1]
                plotting_dict['cluster'].append(cluster_number)
                plotting_dict['ratio'].append(fold_fold_ind)
                induction_vals = [round(x, 1) for x in induction_vals]
                repressed_vals = [round(x, 1) for x in repressed_vals]
                induced_vals = [round(x, 1) for x in induced_vals]

                plotting_dict['text'].append('{0} ({1}, {2}, {3})<br>{4} ({5}, {6}, {7})'.format(
                    seq1, repressed_vals[0], induced_vals[0], induction_vals[0],
                    seq2, repressed_vals[1], induced_vals[1], induction_vals[1]))
                sys.stdout.write('{0}\t{1}\t{2}\t{3}\t{4}\n'.format(
                    seq1, repressed_vals[0], induced_vals[0], induction_vals[0], fold_fold_ind))
                sys.stdout.write('{0}\t{1}\t{2}\t{3}\t{4}\n'.format(
                    seq2, repressed_vals[1], induced_vals[1], induction_vals[1], fold_fold_ind))
                sys.stdout.write('\n')
    protein1 = data_file1.split('/')[-1].split('.')[0]
    protein2 = data_file2.split('/')[-1].split('.')[0]
    fig = {
        'data': [
            {
                'x': plotting_dict['cluster'],
                'y': plotting_dict['ratio'],
                'text': plotting_dict['text'],
                'mode': 'markers',
            }
            ],
        'layout': {
            'xaxis': {'title': 'cluster number'},
            'yaxis': {'title': 'fold fold induction'},
            'title': '{0} over {1}'.format(protein1, protein2)}
    }

    plotly.offline.plot(fig, filename='{0}'.format(name), auto_open=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="script generates scatter plot of induction data from "
                                                 "genotype_fluorescence.py")
    parser.add_argument("-n", "--name", default='cluster_comparison.html',
                        help='name of output html')
    parser.add_argument('-e', '--exact', action='store_true',
                        help='force sequences to be exactly the same')

    required = parser.add_argument_group('required arguments')
    required.add_argument("-d1", "--data_file1", required=True,
                          help="columns of data separated by white space. columns used to generate plot must be "
                               "named the same as the x and y axes")
    required.add_argument("-d2", "--data_file2", required=True,
                          help="columns of data separated by white space. columns used to generate plot must be "
                               "named the same as the x and y axes")

    parser.add_argument('-c', '--cluster',
                        help='a sorted cluster file from cd-hit')
    parser.add_argument('-f1', '--fasta1',
                        help='the fasta file containing sequences appearing in the sorted cluster file')
    parser.add_argument('-f2', '--fasta2',
                        help='the fasta file containing sequences appearing in the sorted cluster file')
    args = parser.parse_args()
    induction_cluster_plotly(args.data_file1, args.data_file2, args.cluster, args.fasta1, args.fasta2, args.name,
                             args.exact)
