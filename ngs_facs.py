#!/usr/bin/env python
import argparse
import collections
import gzip
import itertools
import numpy as np
import pandas as pd
import re
import sys
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


def update_nested_counter_with_fastq(fastq_file, barcode_df, name_state_seq_barcode_counter, counters,
                                     reverse=False, indexed=False):
    if reverse:
        direction = 'reverse'
        barcode_const_array = barcode_df['Barcode Constant Region Reverse'].unique()
    else:
        direction = 'forward'
        barcode_const_array = barcode_df['Barcode Constant Region Forward'].unique()
    library_const_array = barcode_df['Library Constant Region'].unique()
    with gzip.open(fastq_file, 'r') as f:
        for identifier, sequence, spacer, quality_str in itertools.izip_longest(*[f] * 4, fillvalue=None):
            if reverse:
                sequence = str(Seq(sequence, generic_dna).reverse_complement())
            counters[direction]['total'] += 1

            barcode_match_list = []
            if indexed:
                barcode_match_list.append(identifier.rstrip().split(':')[-1])
                counters[direction]['with barcode'] += 1
            else:
                for barcode_const_region in barcode_const_array:
                    barcode_before, barcode_after = re.split('X+', barcode_const_region.upper())
                    barcode_matches = re.findall('{0}(.*){1}'.format(barcode_before, barcode_after), sequence)

                    for barcode_match in barcode_matches:
                        if barcode_match in barcode_df.index:
                            barcode_match_list.append(barcode_match)
                            counters[direction]['with barcode'] += 1

            for barcode_match in barcode_match_list:
                for library_const_region in library_const_array:
                    library_before, library_after = re.split('N+', library_const_region.upper())
                    library_matches = re.findall('{0}(.*){1}'.format(
                        library_before, library_after), sequence)

                    for library_match in library_matches:
                        counters[direction]['with library'] += 1
                        partial_df = barcode_df.loc[barcode_df['Library Constant Region'] == library_const_region]
                        try:
                            expected_length = partial_df.loc[barcode_match]['Expected Length(s)']
                        except KeyError:
                            sys.stdout.write('sequence with inconsistent barcode and library:\n')
                            sys.stdout.write(sequence)
                            sys.stdout.write('barcode:\n{0}\n'.format(barcode_match))
                            sys.stdout.write(partial_df.to_string())
                            sys.stdout.write('\n')
                            continue
                        name = partial_df.loc[barcode_match]['Name']
                        state = partial_df.loc[barcode_match]['State']
                        index = np.where(
                            (barcode_df.index == barcode_match) &
                            (barcode_df['Name'] == name) &
                            (barcode_df['State'] == state)
                        )[0]

                        if expected_length != expected_length:
                            barcode_df.iloc[index, -1] += 1
                            name_state_seq_barcode_counter[name][state][library_match][barcode_match] += 1
                        else:
                            expected_lengths = map(int, expected_length.split(','))
                            if len(library_match) in expected_lengths:
                                barcode_df.iloc[index, -1] += 1
                                name_state_seq_barcode_counter[name][state][library_match][barcode_match] += 1
    return name_state_seq_barcode_counter, barcode_df, counters


def calculate_fluorescence_from_nested_counter(name_state_seq_barcode_counter, barcode_df, minimum_read_count):
    drop_count = 0
    name_matrices_list_dict = collections.defaultdict(dict)
    for name in name_state_seq_barcode_counter.iterkeys():
        fluor_vec_list = []
        count_matrix_dict = {}
        for state in name_state_seq_barcode_counter[name].iterkeys():
            for seq in name_state_seq_barcode_counter[name][state].keys():
                seq_count = sum(name_state_seq_barcode_counter[name][state][seq].itervalues())
                if seq_count < minimum_read_count:
                    del name_state_seq_barcode_counter[name][state][seq]
                    drop_count += 1

            count_matrix = pd.DataFrame.from_dict(name_state_seq_barcode_counter[name][state], orient='index')
            fractions = barcode_df.loc[
                (barcode_df['Name'] == name) & (barcode_df['State'] == state)]['Population Percent']
            fractions = fractions.div(fractions.sum())
            median_fluors = barcode_df.loc[
                (barcode_df['Name'] == name) & (barcode_df['State'] == state)]['Median Fluorescence']
            fraction_matrix = fractions * count_matrix.div(count_matrix.sum(axis='index'), axis='columns')
            norm_fraction_matrix = fraction_matrix.div(fraction_matrix.sum(axis='columns'), axis='index')
            log_fluor_matrix = norm_fraction_matrix * np.log10(median_fluors)
            calculated_fluor_vec = np.power(10, log_fluor_matrix.sum(axis='columns'))
            calculated_fluor_vec.rename(state, inplace=True)
            fluor_vec_list.append(calculated_fluor_vec)
            count_matrix_dict[state] = count_matrix

        fluor_matrix = pd.concat(fluor_vec_list, axis='columns')
        overall_count_matrix = pd.concat(count_matrix_dict, axis='columns')
        overall_count_matrix['sum'] = overall_count_matrix.sum(axis='columns')
        fluor_matrix['count'] = overall_count_matrix['sum']
        overall_count_matrix.sort_values(by='sum', ascending=False, inplace=True)
        fluor_matrix.sort_values(by='count', ascending=False, inplace=True)
        name_matrices_list_dict[name]['fluor'] = fluor_matrix
        name_matrices_list_dict[name]['count'] = overall_count_matrix
    return name_matrices_list_dict, drop_count


def ngs_facs(forward_fastq, barcode_library_mapping_csv, reverse_fastq, minimum_read_count, indexed):
    barcode_df = pd.read_csv(barcode_library_mapping_csv, index_col=0)
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

    name_matrices_dict, drop_count = calculate_fluorescence_from_nested_counter(
        name_state_seq_barcode_counter, barcode_df, minimum_read_count)

    for variable, count in counters['forward'].items():
        sys.stdout.write('seqs {0}: {1}\n'.format(variable, count))
    if reverse_fastq:
        for variable, count in counters['reverse'].items():
            sys.stdout.write('seqs {0}: {1}\n'.format(variable, count))
    sys.stdout.write('seqs with less than min read count: {0}\n'.format(drop_count))

    for name in name_matrices_dict.keys():
        fluor_df = name_matrices_dict[name]['fluor']
        fluor_df.to_csv('{0}_fluorescence.csv'.format(name))
        count_df = name_matrices_dict[name]['count']
        count_df.to_csv('{0}_counts.csv'.format(name))

    barcode_df.to_csv('{0}_counts.csv'.format(barcode_library_mapping_csv.split('/')[-1].split('.csv')[0]))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""
        script to calculate the fluorescence value for nucleotide sequences present in sequencing data 
        from a FACS'ed library.
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
    parser.add_argument('-r', '--reverse_fastq',
                        help='reverse read fastq.gz file. Reverse reads are converted to their reverse complement '
                             '(ie the forward read) before barcode regions and library regions are matched')
    parser.add_argument('-m', '--minimum_read_count', default=5, type=int,
                        help='minimum number of reads a sequence must have in a state to be included in fluorescence '
                             'calculations')
    parser.add_argument('-i', '--indexed', action='store_true',
                        help='use if ngs run was indexed and reads have been assigned a barcode in their sequence '
                             'identifier line')
    args = parser.parse_args()
    ngs_facs(args.forward_fastq, args.barcode_library_mapping, args.reverse_fastq, args.minimum_read_count,
             args.indexed)
