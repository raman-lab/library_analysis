#!/usr/bin/env python
import argparse
import collections
import gzip
import itertools
import numpy as np
import pandas as pd
import re


def parse_barcode_descriptor_file(barcode_descriptor_file):
    barcode_desc_dict = {}
    with open(barcode_descriptor_file, 'r') as f:
        f.next()
        for line in f:
            barcode, desc = line.rstrip().split('\t')
            barcode_desc_dict[barcode.upper()] = desc
    print barcode_desc_dict.keys()
    return barcode_desc_dict


def quality_check_fastq(quality_str):
    expected_incorrect_num = 0
    for ascii_char in quality_str.rstrip():
        phred_score = ord(ascii_char) - 33
        prob_incorrect = pow(10, - float(phred_score) / 10)
        expected_incorrect_num += prob_incorrect
    return expected_incorrect_num


def parse_from_library_architecture(library_architecture):
    with open(library_architecture, 'r') as f:
        library_read = f.readline().rstrip().upper()
    library_start = library_read.find('N')
    library_end = library_read.rfind('N')
    barcode_start = library_read.find('X')
    barcode_end = library_read.rfind('X')

    lib_start_str = library_read[library_start - 6:library_start]
    lib_end_str = library_read[library_end + 1: library_end + 7]
    barcode_start_str = library_read[barcode_start - 6: barcode_start]
    barcode_end_str = library_read[barcode_end + 1: barcode_end + 7]

    return (lib_start_str, lib_end_str), (barcode_start_str, barcode_end_str)


def parse_fastqs_to_dict(fastq_files, max_expected_incorrect, const_library_tup, const_barcode_tup, barcode_desc_dict):
    expected_errors = []
    total_sequence_count = 0
    sequences_passing_qc = 0
    sequences_with_matching_const_region = 0
    sequence_dict = collections.defaultdict(collections.Counter)
    seq_size_dict = collections.defaultdict(collections.Counter)
    for fastq_file in fastq_files:
        with gzip.open(fastq_file, 'r') as f:
            for identifier, sequence, spacer, quality_str in itertools.izip_longest(*[f] * 4, fillvalue=None):
                total_sequence_count += 1
                expected_incorrect = quality_check_fastq(quality_str)
                expected_errors.append(expected_incorrect)
                if expected_incorrect <= max_expected_incorrect:
                    library_seq = re.search('{0}(.*){1}'.format(const_library_tup[0], const_library_tup[1]), sequence)
                    barcode_seq = re.search('{0}(.*){1}'.format(const_barcode_tup[0], const_barcode_tup[1]), sequence)
                    sequences_passing_qc += 1
                    if library_seq and barcode_seq: #and barcode_seq in barcode_desc_dict.keys():
                        library_seq = library_seq.group(1)
                        barcode_seq = barcode_seq.group(1)
                        sequence_dict[barcode_seq][library_seq] += 1
                        sequences_with_matching_const_region += 1
                    if barcode_seq in barcode_desc_dict.keys():
                        seq_size_dict[barcode_desc_dict[barcode_seq]][len(library_seq)] += 1

    barcode_read_df = pd.DataFrame.from_dict(seq_size_dict, orient='index')
    barcode_read_df['total'] = barcode_read_df.sum(axis='columns')
    df_sum = barcode_read_df.sum(axis='index')
    df_min = barcode_read_df.min(axis='index')
    df_max = barcode_read_df.max(axis='index')
    barcode_read_df.loc['sum'] = df_sum
    barcode_read_df.loc['max'] = df_max
    barcode_read_df.loc['min'] = df_min
    errors = np.array(expected_errors)
    med = np.median(errors)
    mad = np.median(np.abs(errors - med))
    matched_percent = float(sequences_with_matching_const_region) / total_sequence_count
    with open('barcode_read_counts.txt', 'w') as o:
        o.write('{0}\n'.format(barcode_read_df.to_string()))
        o.write('\n')
        o.write('total sequences: {0}\n'.format(total_sequence_count))
        o.write('median expected number of wrong bases per read: {0} +/- {1}\n'.format(med, mad))
        o.write('seqs passing qc: {0}\n'.format(sequences_passing_qc))
        o.write('seqs matching const areas: {0}\n'.format(sequences_with_matching_const_region))
        o.write('percentage of matched seqs: {0}\n'.format(matched_percent))

    return sequence_dict


def main(fastq_files, barcode_descriptor_file, library_architecture, expected_error_max, minimum):
    barcode_desc_dict = parse_barcode_descriptor_file(barcode_descriptor_file)
    const_library_tup, const_barcode_tup = parse_from_library_architecture(library_architecture)
    barcode_seq_dict = parse_fastqs_to_dict(fastq_files, expected_error_max, const_library_tup, const_barcode_tup,
                                            barcode_desc_dict)
    for barcode, desc in barcode_desc_dict.iteritems():
        output_lines = []
        for s, seq in enumerate(barcode_seq_dict[barcode].keys()):
            output_lines.append('>seq_{0}\n'.format(s))
            output_lines.append('{0}\n'.format(seq))
        with open('{0}.fasta'.format(desc), 'w') as o:
            o.writelines(output_lines)

        output_df = pd.DataFrame.from_dict(barcode_seq_dict[barcode], orient='index')
        output_df.to_json('{0}.json'.format(desc))
        output_df.to_csv('{0}.csv'.format(desc))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""script to parse library sequences from fastq file.
            input is fastq(s) of sequences and a tab delimited text file of barcodes and descriptors.
            for each barcode, a fasta file and pandas dataframe is output in json and csv format"""
    )
    required = parser.add_argument_group('required arguments')
    required.add_argument('-f', '--fastq', nargs='*', required=True, help='fastq.gz file(s) of sequences')
    required.add_argument('-b', '--barcodes', required=True,
                          help='tab delimited text file with one column for barcodes and one with one word descriptors')
    required.add_argument('-l', '--library_architecture', required=True,
                          help='text file with an example one-line library read used to specify sequence architecture.'
                               'library sequence should be denoted with "N"s, and barcodes should be denoted '
                               'with "X"s. currently this script is takes six residues before and after the library '
                               'and barcode sequences, assumed constant, to extract library and barcode sequences.'
                               'more functionality will be added as needed.')
    parser.add_argument('-e', '--expected_error_max', default=1, type=int,
                        help='integer number of the number of expected wrong bases per read based on phred quality '
                             'score')
    parser.add_argument('-m', '--minimum', default=5, type=int,
                        help='minimum number of reads a sequence must have to be included in fluorescence calculations')
    args = parser.parse_args()
    main(args.fastq, args.barcodes, args.library_architecture, args.expected_error_max, args.minimum)