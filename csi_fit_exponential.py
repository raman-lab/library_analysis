#!/usr/bin/env python
import argparse
import json

import sys

from induction_cluster_plotly import parse_fasta_to_tag_seq_dict, parse_cdhit_to_cluster_tag_dict


def score_csi_from_cdhit(cluster_files, fasta_files, name, minimum_size):
    seq_cluster_size_dict = {}
    for cluster_file, fasta_file in zip(cluster_files, fasta_files):
        temp_dict = {}
        tag_seq_dict = parse_fasta_to_tag_seq_dict(fasta_file)
        cluster_tag_list_dict = parse_cdhit_to_cluster_tag_dict(cluster_file)
        total_count = 0
        for cluster_number, tag_list in cluster_tag_list_dict.iteritems():
            cluster_size = len(tag_list)
            if cluster_size >= minimum_size:
                total_count += cluster_size
                for tag in tag_list:
                    temp_dict[tag_seq_dict[tag]] = float(cluster_size)

        for seq, count in temp_dict.iteritems():
            temp_dict[seq] = count / total_count
        seq_cluster_size_dict.update(temp_dict)

    with open(name, 'w') as o:
        json.dump(seq_cluster_size_dict, o)

    fasta_prefix = name.rstrip('.json')
    for s, seq in enumerate(seq_cluster_size_dict.iterkeys()):
        sys.stdout.write('>{0}_{1}\n'.format(fasta_prefix, s))
        sys.stdout.write('{0}\n'.format(seq))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-n", "--name", default='csi_score.json',
                        help='name of output file. default is csi_score.json')

    required = parser.add_argument_group('required arguments')
    parser.add_argument('-m', '--minimum_size', type=int, default=5,
                        help='minimum cluster size to consider')
    parser.add_argument('-c', '--cluster', nargs='*',
                        help='a sorted cluster file from cd-hit')
    parser.add_argument('-f', '--fasta', nargs='*',
                        help='the fasta file containing all sequences appearing in the sorted cluster file')
    args = parser.parse_args()
    score_csi_from_cdhit(args.cluster, args.fasta, args.name, args.minimum_size)
