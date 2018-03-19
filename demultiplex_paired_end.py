#!/usr/bin/env python
import argparse
import collections
import gzip
import itertools
import numpy as np
import pandas as pd
import re


class ArchitectureMatch(object):
    def __init__(self, start_string, end_string, start_offset_int, end_offset_int, name_string):
        self._start_str = None
        self._end_str = None
        self._start_offset = None
        self._end_offset = None
        self._name = None
        self.start_str = start_string
        self.end_str = end_string
        self.start_offset = start_offset_int
        self.end_offset = end_offset_int
        self.name = name_string

    @property
    def start_str(self):
        return self._start_str

    @start_str.setter
    def start_str(self, start_string):
        self._start_str = str(start_string).upper()

    @property
    def end_str(self):
        return self._end_str

    @end_str.setter
    def end_str(self, end_string):
        self._end_str = str(end_string.upper())

    @property
    def start_offset(self):
        return self._start_offset

    @start_offset.setter
    def start_offset(self, start_offset_int):
        self._start_offset = int(start_offset_int)

    @property
    def end_offset(self):
        return self._end_offset

    @end_offset.setter
    def end_offset(self, end_offset_int):
        self._end_offset = int(end_offset_int)

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name_str):
        self._name = str(name_str)


def parse_barcode_descriptor_file(barcode_descriptor_file):
    label_dict = collections.OrderedDict([('barcode', None), ('name', None), ('sub_arc', None)])
    barcode_desc_dict = {}

    with open(barcode_descriptor_file, 'r') as f:
        first_line_list = f.next().rstrip().lower().split('\t')
        labels = label_dict.keys()
        for index, item in enumerate(first_line_list):
            for label in labels:
                if label in item:
                    label_dict[label] = index
                    labels.remove(label)
                    break

        for line in f:
            if not line.isspace() and not line.startswith('#'):
                line_list = line.rstrip().split('\t')
                barcode = line_list[label_dict['barcode']].upper()
                desc = line_list[label_dict['name']]
                if label_dict['sub_arc']:
                    sub_arc = line_list[label_dict['sub_arc']]
                else:
                    sub_arc = None
                barcode_desc_dict[barcode.upper()] = {'name': desc, 'sub_arc': sub_arc}
    return barcode_desc_dict


def quality_check_fastq(quality_str):
    expected_incorrect_num = 0
    for ascii_char in quality_str.rstrip():
        phred_score = ord(ascii_char) - 33
        prob_incorrect = pow(10, - float(phred_score) / 10)
        expected_incorrect_num += prob_incorrect
    return expected_incorrect_num


def parse_library_architecture(library_architecture_file):
    with open(library_architecture_file, 'r') as f:
        library_architecture = f.readline().rstrip().upper()
    # library_start = library_read.find('N')
    # library_end = library_read.rfind('N')
    # barcode_start = library_read.find('X')
    # barcode_end = library_read.rfind('X')
    #
    # lib_start_str = library_read[library_start - 6:library_start]
    # lib_end_str = library_read[library_end + 1: library_end + 7]
    # barcode_start_str = library_read[barcode_start - 6: barcode_start]
    # barcode_end_str = library_read[barcode_end + 1: barcode_end + 7]

    # find barcodes - max 2
    # barcode_iterator = re.finditer('0*X+0*', library_architecture)
    # for instance in barcode_iterator:
    #     barcode_start = instance.start()
    #     barcode_end = instance.end()
    #     barcode_string = instance.group(0)

    # find library - max 1


    # list of indices that character occurs at    l = [m.start() for m in re.finditer('X', library_read)]
    # split list based on gaps np.split(l, np.add(np.where(np.diff(l) > 1), 1).tolist()[0])
    # return (lib_start_str, lib_end_str), (barcode_start_str, barcode_end_str)


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
    # barcode_desc_dict = parse_barcode_descriptor_file(barcode_descriptor_file)

    # decided to not make this general for now. going to hard code architectures because it is way less of a pain
    # const_library_tup, const_barcode_tup = parse_library_architecture(library_architecture)

    # barcode_seq_dict = parse_fastqs_to_dict(fastq_files, expected_error_max, const_library_tup, const_barcode_tup,
    #                                         barcode_desc_dict)
    # for barcode, desc in barcode_desc_dict.iteritems():
    #     output_lines = []
    #     for s, seq in enumerate(barcode_seq_dict[barcode].keys()):
    #         output_lines.append('>seq_{0}\n'.format(s))
    #         output_lines.append('{0}\n'.format(seq))
    #     with open('{0}.fasta'.format(desc), 'w') as o:
    #         o.writelines(output_lines)
    #
    #     output_df = pd.DataFrame.from_dict(barcode_seq_dict[barcode], orient='index')
    #     output_df.to_json('{0}.json'.format(desc))
        # output_df.to_csv('{0}.csv'.format(desc))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""script to parse library sequences from fastq file.
            input is fastq(s) of sequences and a tab delimited text file of barcodes and descriptors.
            for each barcode, a fasta file and pandas dataframe is output in json and csv format"""
    )
    required = parser.add_argument_group('required arguments')
    required.add_argument('-f', '--forward_fastq', required=True, help='forward read fastq.gz file')
    required.add_argument('-r', '--reverse_fastq', required=True, help='reverse read fastq.gz file')
    required.add_argument('-b', '--barcodes', required=True,
                          help='tab delimited text file with one column for barcodes and one with one word names.'
                               'can optionally give a sub-architecture if there is more than one level of barcoding')
    required.add_argument('-l', '--library_architecture', required=True,
                          help='text file with an example one-line library read used to specify sequence architecture.'
                               'library sequence should be denoted with "N"s, and barcodes should be denoted '
                               'with "X"s. currently this script is takes c residues before and after the library '
                               'and barcode sequences, assumed constant, to extract library and barcode sequences.'
                               'can have "*"s in architecture to specify random bases before or after lib or barcodes')
    parser.add_argument('-e', '--expected_error_max', default=1, type=int,
                        help='integer number of the number of expected wrong bases per read based on phred quality '
                             'score')
    parser.add_argument('-m', '--minimum', default=5, type=int,
                        help='minimum number of reads a sequence must have to be included in fluorescence calculations')
    parser.add_argument('-c', '--constant_to_match', default=6, type=int,
                        help='number of bases to match on either side of library and barcode sequences')
    args = parser.parse_args()
    main(args.fastq, args.barcodes, args.library_architecture, args.expected_error_max, args.minimum)