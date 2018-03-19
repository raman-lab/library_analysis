#!/usr/bin/env python
import argparse
import collections
import itertools
import pandas as pd


def parse_fasta_to_tag_seq_dict(fasta_file):
    tag_sequence_dict = {}
    with open(fasta_file, 'r') as f:
        for tag, sequence in itertools.izip_longest(f, f, fillvalue=None):
            tag = tag.rstrip().split('>', 1)[-1].split(' ')[0]
            sequence = sequence.rstrip()
            tag_sequence_dict[tag] = sequence
    return tag_sequence_dict


def pull_top_blast_hit(input_file, fasta, expected_length, read_count_minimum):
    expected_lengths = {
        'TetR_seg1': 115,
        'TetR_seg2': 115,
        'TetR_seg3': 115,
        'TetR_seg4': 115,
        'TetR_seg5': 114,
        'TetR_seg6': 48,
        'TtgR_seg1': 114,
        'TtgR_seg2': 114,
        'TtgR_seg3': 115,
        'TtgR_seg4': 114,
        'TtgR_seg5': 122,
        'TtgR_seg6': 62,
        'TetR_B0': 129,
        'TetR_B1': 127,
        'TtgR_B0': 124,
        'TtgR_B1': 129
    }
    id_seq_dict = parse_fasta_to_tag_seq_dict(fasta)
    match_df = pd.DataFrame(data=None, columns=['library_match', 'percent_identity', 'length_difference',
                                                'read_count'])
    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith('Query='):
                query = line.rstrip().split('Query= ')[-1]
                f.next()
                length = int(f.next().rstrip().split('Length=')[-1])
                gathering = True
                while gathering:
                    line = f.next().rstrip()
                    if line.startswith('***** No hits found *****'):
                        gathering = False
                    elif line.startswith('Sequences producing significant alignments:'):
                        f.next()
                        top_hit = f.next().rstrip().split()[0]
                    elif line.startswith(' Identities ='):
                        identities_split = line.rstrip().split('/')
                        first_identity_count = identities_split[0].split(' ')[-1]
                        second_identity_count = identities_split[1].split(' ')[0]
                        length_diff = abs(expected_length - length)

                        if first_identity_count == second_identity_count and length_diff == 0:
                            fraction = 1.0
                        else:
                            fraction = float(first_identity_count) / length

                        read_count = int(query.split('_')[-1])
                        if fraction > 0.95 and length_diff < 10 and read_count >= read_count_minimum:
                            s = pd.Series(
                                name=id_seq_dict[query],
                                data={'library_match': top_hit, 'percent_identity': fraction,
                                      'length_difference': length_diff, 'read_count': read_count}
                            )
                            match_df = match_df.append(s)
                        gathering = False
    match_df.sort_values(['percent_identity', 'length_difference', 'read_count'], ascending=[False, True, False],
                         inplace=True)

    match_df.to_csv('{0}.csv'.format(input_file.split('.')[0]))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""script for retrieving the top sequence alignment from blast search. output is formatted in
        text file for import into excel. assumes naming convention for query sequences"""
    )
    parser.add_argument('-r', '--read_count', default=1, type=int, help='minimum read count to be included in output')
    required = parser.add_argument_group('required arguments')
    required.add_argument('-i', '--input', required=True, help='file from blast search')
    required.add_argument('-l', '--expected_length', type=int, required=True,
                          help='expected length of query sequences. queries within a 5 nt length of this will be '
                               'considered in output')
    required.add_argument('-f', '--fasta',
                          help='fasta file used to create the query sequences for blast search')
    args = parser.parse_args()
    pull_top_blast_hit(args.input, args.fasta, args.expected_length, args.read_count)
