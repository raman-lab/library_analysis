#!/usr/bin/env python

# this file is part of the Raman Lab github repository: https://github.com/raman-lab/library_analysis
# author: nwhoppe
# created: 8/8/18

import argparse
import collections
import gzip

import itertools
import numpy as np
import pandas as pd
import re
import sys
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq

from barcode_hamming_distance import compute_hamming_distance


def probability_incorrect_from_ascii(ascii_char, base=33):
    phred_score = ord(ascii_char) - base
    probability_incorrect = pow(10, -float(phred_score) / 10)
    return probability_incorrect


def expected_incorrect_from_quality(quality_str):
    """Calculate expected number of incorrect nucleotides from ASCII quality string"""
    expected_incorrect_num = 0
    for ascii_char in quality_str.rstrip():
        expected_incorrect_num += probability_incorrect_from_ascii(ascii_char)
    return expected_incorrect_num


def str_between_const_characters(search_str, const_str_array, variable_character='X'):
    match_const_chars_dict = {}
    for const_str_region in const_str_array:
        const_before_str, const_after_str = re.split('{0}+'.format(variable_character), const_str_region.upper())
        str_matches = re.findall('{0}(.*){1}'.format(const_before_str, const_after_str), search_str)
        match_const_chars_dict[const_str_region] = str_matches
    return match_const_chars_dict


def multiple_match_error(identifier, match_list1, match_list2, description_str):
    sys.stderr.write('Error: more than one {0} for read {1}\n'.format(description_str, identifier))
    sys.stderr.write('Forward matches found: {0}\n'.format(' '.join(match_list1)))
    sys.stderr.write('Reverse matches found: {0}\n'.format(' '.join(match_list2)))
    sys.stderr.write('For the read to be assigned, only one match can be present\n\n')


def correct_barcode(barcode_pair_match, barcode_list, max_correction_distance):
    hamming_distances = []
    for barcode_pair in barcode_list:
        if len(barcode_pair_match) == len(barcode_pair):
            hamming_distances.append(compute_hamming_distance(barcode_pair_match, barcode_pair))
        else:
            hamming_distances.append(1e6)
    min_indices = [i for i, x in enumerate(hamming_distances) if x == min(hamming_distances)]
    min_hamming_distance = hamming_distances[min_indices[0]]
    corrected_barcode = barcode_list[min_indices[0]]
    if len(min_indices) > 1 or min_hamming_distance >= max_correction_distance:
        return None, None
    else:
        return corrected_barcode, min(hamming_distances)


def ascii_from_probability_incorrect(probability_incorrect, base=33):
    phred_score = int(-10 * np.log10(probability_incorrect))
    ascii_char = chr(phred_score + base)
    return ascii_char


def correct_library(library_fwd_match, library_rev_match, quality_fwd_match, quality_rev_match):
    corrected_library_list = []
    corrected_quality_list = []
    changed_positions_counter = 0
    for index, fwd_rev_nt in enumerate(itertools.izip(library_fwd_match, library_rev_match)):
        if fwd_rev_nt[0] == fwd_rev_nt[1]:
            corrected_library_list.append(fwd_rev_nt[0])
            p_incorrect_fwd = probability_incorrect_from_ascii(quality_fwd_match[index])
            p_incorrect_rev = probability_incorrect_from_ascii(quality_rev_match[index])
            p_adjusted = p_incorrect_fwd * p_incorrect_rev
            corrected_quality_list.append(ascii_from_probability_incorrect(p_adjusted))
        else:
            if quality_fwd_match[index] == quality_rev_match[index]:
                return None, None, None
            elif quality_fwd_match[index] > quality_rev_match[index]:
                corrected_library_list.append(library_fwd_match[index])
                corrected_quality_list.append(quality_fwd_match[index])
                changed_positions_counter += 1
            else:
                corrected_library_list.append(library_rev_match[index])
                corrected_quality_list.append(quality_rev_match[index])
                changed_positions_counter += 1
    return ''.join(corrected_library_list), ''.join(corrected_quality_list), changed_positions_counter


def group_seq_counter_from_paired_fastqs(fastq_fwd, fastq_rev, barcode_df, group_seq_counter, read_counters,
                                         expected_error_max, max_correction_distance, indexed=False):
    """Extract barcodes and library sequences from reads and count occurrences by group and subgroup"""
    expected_errors = np.array([])
    hamming_distances = np.array([])
    base_changes = np.array([])
    const_barcode_fwd_array = barcode_df['Forward Constant Region'].unique()
    const_barcode_rev_array = barcode_df['Reverse Constant Region'].unique()
    const_library_array = barcode_df['Library Constant Region'].unique()

    barcode_list = ['_'.join(tup) for tup in zip(
        barcode_df['Forward Barcode'].str.upper(), barcode_df['Reverse Barcode'].str.upper())]

    with gzip.open(fastq_fwd, 'r') as r1, gzip.open(fastq_rev, 'r') as r2:
        for identifier_fwd, identifier_rev in itertools.izip(r1, r2):
            read_counters['total'] += 1

            read_fwd = r1.next().rstrip()
            r1.next()
            quality_fwd = r1.next().rstrip()

            read_rev = str(Seq(r2.next().rstrip(), generic_dna).reverse_complement())
            r2.next()
            quality_rev = r2.next().rstrip()[::-1]

            if identifier_fwd.split()[0] != identifier_rev.split()[0]:
                sys.stderr.write('Error: forward and reverse identifiers do not match\n')
                sys.stderr.write('forward identifier: {0}\n'.format(identifier_fwd))
                sys.stderr.write('reverse identifier: {0}\n'.format(identifier_rev))
                raise NameError('if forward and reverse reads do not appear in the same order in each file, '
                                'this script will be unable to match reads\n\n')
            else:
                common_id = identifier_fwd.split()[0]
            # PUT ERROR CHECK AFTER MERGE
            # expected_incorrect = quality_check_fastq(quality_str)
            # np.append(expected_errors, expected_incorrect)
            # expected_errors.append(expected_incorrect)
            # if expected_incorrect > expected_error_max:
            #     read_counters['failing qc'] += 1
            #     continue
            barcode_pair_match = None
            if indexed:
                barcode_fwd = identifier_fwd.rstrip().split(':')[-1]
                barcode_rev = identifier_rev.rstrip().split(':')[-1]
                barcode_pair_match = '_'.join([barcode_fwd, barcode_rev])
                read_counters['perfect barcode'] += 1
            else:
                barcode_fwd_match_dict = str_between_const_characters(read_fwd, const_barcode_fwd_array)
                barcode_rev_match_dict = str_between_const_characters(read_rev, const_barcode_rev_array)

                barcode_fwd_match_dict = {k: v for (k, v) in barcode_fwd_match_dict.items() if v}
                barcode_rev_match_dict = {k: v for (k, v) in barcode_rev_match_dict.items() if v}

                if len(barcode_fwd_match_dict.keys()) > 1 or len(barcode_rev_match_dict.keys()) > 1:
                    multiple_match_error(common_id, barcode_fwd_match_dict.keys(),
                                         barcode_rev_match_dict.keys(), 'barcode constant region')
                    read_counters['ambiguous barcode'] += 1
                elif len(barcode_fwd_match_dict.keys()) == 1 and len(barcode_rev_match_dict.keys()) == 1:
                    barcode_fwd_matches = barcode_fwd_match_dict.values()[0]
                    barcode_rev_matches = barcode_rev_match_dict.values()[0]
                    if len(barcode_fwd_matches) != 1 or len(barcode_rev_matches) != 1:
                        multiple_match_error(common_id, barcode_fwd_matches, barcode_rev_matches, 'barcode')
                        read_counters['ambiguous barcode'] += 1
                    else:
                        barcode_pair_match = '_'.join([barcode_fwd_matches[0], barcode_rev_matches[0]])
                        if barcode_pair_match in barcode_list:
                            read_counters['perfect barcode'] += 1
                        else:
                            barcode_pair_match, hamming_distance = correct_barcode(barcode_pair_match, barcode_list,
                                                                                   max_correction_distance)
                            if barcode_pair_match:
                                hamming_distances = np.append(hamming_distances, hamming_distance)
                                read_counters['corrected barcode'] += 1
                            else:
                                read_counters['ambiguous barcode'] += 1
                else:
                    read_counters['no barcode'] += 1

            library_match = None
            if barcode_pair_match:
                library_fwd_match_dict = str_between_const_characters(read_fwd, const_library_array)
                library_rev_match_dict = str_between_const_characters(read_rev, const_library_array)

                library_fwd_match_dict = {k: v for (k, v) in library_fwd_match_dict.items() if v}
                library_rev_match_dict = {k: v for (k, v) in library_rev_match_dict.items() if v}

                if len(library_fwd_match_dict.keys()) > 1 or len(library_rev_match_dict.keys()) > 1:
                    multiple_match_error(common_id, library_fwd_match_dict.keys(), library_rev_match_dict.keys(),
                                         'library constant region')
                    read_counters['ambiguous library'] += 1
                elif len(library_fwd_match_dict.keys()) == 1 and len(library_rev_match_dict.keys()) == 1 and \
                        library_fwd_match_dict.keys()[0] == library_rev_match_dict.keys()[0]:
                    const_library_match = library_fwd_match_dict.keys()[0]
                    library_fwd_matches = library_fwd_match_dict.values()[0]
                    library_rev_matches = library_rev_match_dict.values()[0]
                    if len(library_fwd_matches) != 1 or len(library_rev_matches) != 1:
                        multiple_match_error(common_id, library_fwd_matches, library_rev_matches, 'library sequence')
                        read_counters['ambiguous library'] += 1
                    else:
                        library_fwd_match = library_fwd_matches[0]
                        library_rev_match = library_rev_matches[0]
                        if library_fwd_match == library_rev_match:
                            library_match = library_fwd_match
                            quality_match = quality_fwd[read_fwd.index(library_fwd_match):
                                                        read_fwd.index(library_fwd_match) + len(library_fwd_match)]
                            read_counters['perfect library'] += 1
                        elif len(library_fwd_match) == len(library_rev_match):
                            quality_fwd_match = quality_fwd[read_fwd.index(library_fwd_match):
                                                            read_fwd.index(library_fwd_match) + len(library_fwd_match)]
                            quality_rev_match = quality_rev[read_rev.index(library_rev_match):
                                                            read_rev.index(library_rev_match) + len(library_rev_match)]
                            library_match, quality_match, change_count = correct_library(
                                library_fwd_match, library_rev_match, quality_fwd_match, quality_rev_match)
                            if library_match:
                                base_changes = np.append(base_changes, change_count)
                                read_counters['corrected library'] += 1
                            else:
                                read_counters['ambiguous library'] += 1
                        else:
                            read_counters['ambiguous library'] += 1
                else:
                    read_counters['no library'] += 1

                if library_match:
                    expected_incorrect = expected_incorrect_from_quality(quality_match)
                    expected_errors = np.append(expected_errors, expected_incorrect)
                    if expected_incorrect > expected_error_max:
                        read_counters['library failing qc'] += 1
                    else:
                        barcode_fwd, barcode_rev = barcode_pair_match.split('_')
                        barcode_series = barcode_df.loc[
                            (barcode_df['Library Constant Region'] == const_library_match) &
                            (barcode_df['Forward Barcode'] == barcode_fwd) &
                            (barcode_df['Reverse Barcode'] == barcode_rev)
                        ].iloc[0]
                        group = barcode_series['Group']
                        subgroup = barcode_series['Subgroup']
                        barcode_df.loc[barcode_series.name, 'Count'] += 1
                        group_seq_counter[group][subgroup][library_match] += 1

    sys.stdout.write('Expected errors per read: {0} +/- {1} (n = {2})\n'.format(
        np.average(expected_errors), np.std(expected_errors), len(expected_errors)))
    sys.stdout.write('Hamming distance per corrected barcode: {0} +/- {1} (n = {2})\n'.format(
        np.average(hamming_distances), np.std(hamming_distances), len(hamming_distances)))
    sys.stdout.write('Base changes per corrected library: {0} +/- {1} (n = {2})\n'.format(
        np.average(base_changes), np.std(base_changes), len(base_changes)))
    return group_seq_counter, barcode_df, read_counters


def ngs_paired_end_counts(forward_fastq, reverse_fastq, barcode_library_mapping, expected_error_max, indexed,
                          max_correction_distance):
    barcode_df = pd.read_csv(barcode_library_mapping, index_col=False)
    barcode_df.dropna(inplace=True)
    columns_to_capitalize = ['Forward Barcode', 'Reverse Barcode']
    for column in columns_to_capitalize:
        barcode_df[column] = barcode_df[column].str.upper()
    count_series = pd.Series(0, index=barcode_df.index, name='Count')
    barcode_df = pd.concat([barcode_df, count_series], axis='columns')

    counter_list = ['total', 'perfect barcode', 'corrected barcode', 'ambiguous barcode', 'no barcode',
                    'perfect library', 'library failing qc', 'no library', 'ambiguous library', 'corrected library']
    read_counters = {key: 0 for key in counter_list}

    group_seq_counter = collections.defaultdict(dict)
    for group, subgroup in zip(*[barcode_df['Group'], barcode_df['Subgroup']]):
        group_seq_counter[group][subgroup] = collections.Counter()

    group_seq_counter, barcode_df, counters = group_seq_counter_from_paired_fastqs(
        forward_fastq, reverse_fastq, barcode_df, group_seq_counter, read_counters, expected_error_max,
        max_correction_distance, indexed=indexed
    )

    for group, subgroup_dict in group_seq_counter.items():
        count_df = pd.DataFrame.from_dict(subgroup_dict, orient='index')
        if not count_df.empty:
            count_df.to_csv('{0}_counts.csv'.format(group))

    barcode_df.to_csv('{0}_counts.csv'.format(barcode_library_mapping.split('/')[-1].split('.csv')[0]))

    for read_type, count in read_counters.items():
        sys.stdout.write('Reads {0}: {1}\n'.format(read_type, count))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""
    script to extract nucleotide sequences from a paired fastq files using constraints supplied in csv format and 
    output counts per sequence
        required input:
            forward and reverse fastq file (typically have a R1 and R2 in filename, respectively)
            a csv file made from the template barcode_library_mapping.csv
        optional input:
            expected error limit of read
            indicate if ngs run was indexed and reads have already been assigned a barcode
    """)
    required = parser.add_argument_group('required')
    required.add_argument('-f', '--forward_fastq', required=True, help='forward read fastq.gz file')
    required.add_argument('-r', '--reverse_fastq', required=True,
                          help='reverse read fastq.gz file. Reverse reads are converted to their reverse complement '
                               '(ie the forward read) before barcode regions and library regions are matched')
    required.add_argument('-b', '--barcode_library_mapping', required=True,
                          help='csv file with the columns specified in barcode_library_mapping.csv')
    parser.add_argument('-e', '--expected_error_max', default=1, type=int,
                        help='integer number of the number of expected wrong bases in extracted library sequences '
                             'based on phred quality score')
    parser.add_argument('-i', '--indexed', action='store_true',
                        help='use if ngs run was indexed and reads have been assigned a barcode in their sequence '
                             'identifier line')
    parser.add_argument('-mcd', '--max_correction_distance', default=4, type=int,
                        help='maximum hamming distance allowed during barcode correction. this should be equal to '
                             'the minimum hamming distance between the pairs of barcodes')
    args = parser.parse_args()
    ngs_paired_end_counts(args.forward_fastq, args.reverse_fastq, args.barcode_library_mapping, args.expected_error_max,
                          args.indexed, args.max_correction_distance)
