#!/usr/bin/env python
import argparse
import collections
import pandas as pd
import sys

from ngs_facs import update_nested_counter_with_fastq


def seq_counts_from_nested_counter(name_state_seq_barcode_counter, barcode_df, minimum_read_count):
    drop_count = 0
    name_matrix_dict = {}
    for name in name_state_seq_barcode_counter.iterkeys():
        count_matrix_dict = {}
        for state in name_state_seq_barcode_counter[name].iterkeys():
            for seq in name_state_seq_barcode_counter[name][state].keys():
                seq_count = sum(name_state_seq_barcode_counter[name][state][seq].itervalues())
                if seq_count < minimum_read_count:
                    del name_state_seq_barcode_counter[name][state][seq]
                    drop_count += 1

            count_matrix = pd.DataFrame.from_dict(name_state_seq_barcode_counter[name][state], orient='index')
            count_matrix_dict[state] = count_matrix

        overall_count_matrix = pd.concat(count_matrix_dict, axis='columns')
        overall_count_matrix['sum'] = overall_count_matrix.sum(axis='columns')
        overall_count_matrix.sort_values(by='sum', ascending=False, inplace=True)
        name_matrix_dict[name] = overall_count_matrix
    return name_matrix_dict, drop_count


def ngs_demultiplex(forward_fastq, barcode_library_mapping_csv, reverse_fastq, minimum_read_count, indexed):
    barcode_df = pd.DataFrame.from_csv(barcode_library_mapping_csv)
    barcode_df.index = barcode_df.index.str.upper()
    count_series = pd.Series(0, index=barcode_df.index, name='Count')
    barcode_df = pd.concat([barcode_df, count_series], axis='columns')
    counters = {
        'forward': {'total': 0, 'with barcode': 0, 'with library': 0},
        'reverse': {'total': 0, 'with barcode': 0, 'with library': 0}
    }

    name_state_seq_barcode_counter = collections.defaultdict(
        lambda: collections.defaultdict(
            lambda: collections.defaultdict(
                lambda: collections.Counter()
            )
        )
    )

    name_state_seq_barcode_counter, barcode_df, counters = update_nested_counter_with_fastq(
        forward_fastq, barcode_df, name_state_seq_barcode_counter, counters, indexed=indexed)

    if reverse_fastq:
        name_state_seq_barcode_counter, barcode_df, counters = update_nested_counter_with_fastq(
            reverse_fastq, barcode_df, name_state_seq_barcode_counter, counters, reverse=True, indexed=indexed)

    name_matrix_dict, drop_count = seq_counts_from_nested_counter(
        name_state_seq_barcode_counter, barcode_df, minimum_read_count)

    for variable, count in counters['forward'].items():
        sys.stdout.write('seqs {0}: {1}\n'.format(variable, count))
    if reverse_fastq:
        for variable, count in counters['reverse'].items():
            sys.stdout.write('seqs {0}: {1}\n'.format(variable, count))
    sys.stdout.write('seqs with less than min read count: {0}\n'.format(drop_count))

    for name, count_df in name_matrix_dict.items():
        count_df.to_csv('{0}_counts.csv'.format(name))

    barcode_df.to_csv('{0}_counts.csv'.format(barcode_library_mapping_csv.split('/')[-1].split('.csv')[0]))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""
        script to extract nucleotide sequences from a fastq file using constraints supplied in csv format
        required input:
            forward fastq file from ngs run
            a csv file made from the template barcode_library_mapping.csv
        optional input:
            reverse fastq file from paired end ngs run
            minimum number of reads per sequence
            indicate if ngs run was indexed and reads have already been assigned a barcode
        """
    )
    required = parser.add_argument_group('required arguments')
    required.add_argument('-f', '--forward_fastq', required=True, help='forward read fastq.gz file')
    required.add_argument('-b', '--barcode_library_mapping', required=True,
                          help='csv file with the columns specified in barcode_library_mapping.csv')
    parser.add_argument('-r', '--reverse_fastq', help='reverse read fastq.gz file')
    parser.add_argument('-m', '--minimum_read_count', default=5, type=int,
                        help='minimum number of reads a sequence must have in a state to be included in fluorescence '
                             'calculations')
    parser.add_argument('-i', '--indexed', action='store_true',
                        help='use if ngs run was indexed and reads have been assigned a barcode in their sequence '
                             'identifier line')
    args = parser.parse_args()
    ngs_demultiplex(args.forward_fastq, args.barcode_library_mapping, args.reverse_fastq, args.minimum_read_count,
                    args.indexed)
