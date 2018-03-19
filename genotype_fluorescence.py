#!/usr/bin/env python
import argparse
import collections
import gzip
import itertools
import numpy as np
import pandas as pd
import re
import sys
import time
from functools import wraps
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


def fn_timer(function):
    @wraps(function)
    def function_timer(*args, **kwargs):
        t0 = time.time()
        result = function(*args, **kwargs)
        t1 = time.time()
        print ("Total time running %s: %s seconds" % (function.func_name, str(t1-t0)))
        return result
    return function_timer


class Bin(object):
    def __init__(self, barcode_sequence_str, median_fluor_float, fraction_float):
        self._sequence = None
        self._median_fluor = None
        self._fraction = None
        self.sequence = barcode_sequence_str
        self.median_fluor = float(median_fluor_float)
        self.fraction = float(fraction_float)

    @property
    def sequence(self):
        return self._sequence

    @sequence.setter
    def sequence(self, sequence_str):
        self._sequence = str(sequence_str).upper()

    @property
    def median_fluor(self):
        return self._median_fluor

    @median_fluor.setter
    def median_fluor(self, median_fluor_float):
        self._median_fluor = float(median_fluor_float)

    @property
    def fraction(self):
        return self._fraction

    @fraction.setter
    def fraction(self, fraction_float):
        if fraction_float > 1:
            self._fraction = float(fraction_float) / 100
        else:
            self._fraction = float(fraction_float)


@fn_timer
def parse_barcode_file(barcode_files):
    label_dict = collections.OrderedDict(
        [('barcode', None), ('protein', None), ('state', None), ('fluor', None), ('percent', None), ('licate', None)])
    protein_replicate_state_dict = collections.defaultdict(lambda: collections.defaultdict(dict))
    barcode_seq_list = []
    for barcode_file in barcode_files:
        with open(barcode_file, 'r') as f:
            first_line_list = f.next().rstrip().lower().split('\t')
            labels = label_dict.keys()
            for index, item in enumerate(first_line_list):
                for label in labels:
                    if label in item:
                        label_dict[label] = index
                        labels.remove(label)
                        break
            if None in label_dict.values():
                raise Exception('One of the essential labels in the barcode input file was not found')
            for line in f:
                if not line.isspace() and not line.startswith('#'):
                    line_list = line.rstrip().split('\t')
                    barcode_seq = line_list[label_dict['barcode']].upper()
                    protein = line_list[label_dict['protein']]
                    state = line_list[label_dict['state']].lower()
                    median_fluor = line_list[label_dict['fluor']]
                    percent = line_list[label_dict['percent']]
                    replicate = int(line_list[label_dict['licate']])

                    if 'rep' in state:
                        state = 'repressed'
                    elif 'ind' in state:
                        state = 'induced'
                    elif 'free' in state:
                        state = 'free'
                    else:
                        raise Exception('States can only be repressed, induced, or free')

                    bin_obj = Bin(barcode_seq, median_fluor, percent)
                    barcode_seq_list.append(barcode_seq)

                    if state not in protein_replicate_state_dict[protein][replicate].keys():
                        protein_replicate_state_dict[protein][replicate][state] = []
                    protein_replicate_state_dict[protein][replicate][state].append(bin_obj)

    return protein_replicate_state_dict, barcode_seq_list


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


@fn_timer
def parse_fastqs_to_dict(fastq_files, max_expected_incorrect, const_library_tup, const_barcode_tup, barcode_seq_list):
    expected_errors = []
    total_sequence_count = 0
    sequences_passing_qc = 0
    sequences_with_matching_const_region = 0
    sequence_dict = collections.defaultdict(collections.Counter)
    barcode_read_df = pd.DataFrame(0, index=barcode_seq_list, columns=[16, 17, 18, 19])
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
                    if library_seq and barcode_seq:
                        library_seq = library_seq.group(1)
                        barcode_seq = barcode_seq.group(1)
                        rev_comp_lib_seq = str(Seq(library_seq, generic_dna).reverse_complement())
                        if rev_comp_lib_seq < library_seq:
                            library_seq = rev_comp_lib_seq
                        sequence_dict[barcode_seq][library_seq] += 1
                        sequences_with_matching_const_region += 1

                        if barcode_seq in barcode_seq_list and len(library_seq) in [16, 17, 18, 19]:
                            barcode_read_df.loc[barcode_seq][len(library_seq)] += 1

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


@fn_timer
def calculate_fluorescence_from_seq_barcode_dict(barcode_sequence_dict, protein_replicate_state_dict, minimum_count):
    sys.stdout.write('dropping items from dict\n')
    t0 = time.time()
    drop = 0
    for barcode, seq_counter in barcode_sequence_dict.iteritems():
        for seq, count in itertools.dropwhile(lambda barcode_count: barcode_count[1] >= minimum_count,
                                              seq_counter.most_common()):
            del barcode_sequence_dict[barcode][seq]
            drop += 1
    t1 = time.time()
    sys.stdout.write('dropped {0} items\n'.format(drop))
    sys.stdout.write('time to drop: {0}\n'.format(t1 - t0))

    for protein in protein_replicate_state_dict.keys():
        replicate_matrices = []
        replicate_keys = []
        replicate_count_matrices = []
        for replicate in protein_replicate_state_dict[protein].keys():
            replicate_keys.append('replicate {0}'.format(replicate))
            state_matrices = []
            rename_dict = {}
            state_count_matrices = {}
            for state in protein_replicate_state_dict[protein][replicate].keys():
                state_fractions = []
                state_median_fluors = []
                state_seqs = []
                for bin_obj in protein_replicate_state_dict[protein][replicate][state]:
                    state_fractions.append(bin_obj.fraction)
                    state_median_fluors.append(bin_obj.median_fluor)
                    state_seqs.append(bin_obj.sequence)
                    rename_dict[bin_obj.sequence] = state

                state_median_fluor_vec = pd.Series(state_median_fluors, index=state_seqs)
                state_fraction_vec = pd.Series(state_fractions, index=state_seqs)
                state_fraction_vec = state_fraction_vec.div(state_fraction_vec.sum())
                state_matrix = pd.DataFrame([barcode_sequence_dict[s] for s in state_seqs], index=state_seqs).T

                state_matrix.dropna(axis='index', how='all', inplace=True)
                fraction_matrix = state_fraction_vec * state_matrix.div(state_matrix.sum(axis='index'), axis='columns')
                norm_fraction_matrix = fraction_matrix.div(fraction_matrix.sum(axis='columns'), axis='index')
                log_fluor_matrix = norm_fraction_matrix * np.log10(state_median_fluor_vec)
                seq_vec = log_fluor_matrix.sum(axis='columns')
                state_mean_fluor = np.power(10, seq_vec)
                state_mean_fluor.rename(state, inplace=True)
                state_matrices.append(state_mean_fluor)

                state_count_matrices[state] = state_matrix

            replicate_matrix = pd.concat(state_matrices, axis='columns')
            replicate_matrix['fold_induction'] = replicate_matrix['induced'] / replicate_matrix['repressed']
            replicate_matrix['fold_repression'] = replicate_matrix['free'] / replicate_matrix['repressed']
            replicate_matrices.append(replicate_matrix)

            replicate_count_matrix = pd.concat(state_count_matrices, axis='columns')
            replicate_count_matrix['sum'] = replicate_count_matrix.sum(axis='columns')
            replicate_count_matrices.append(replicate_count_matrix)

        protein_matrix = pd.concat(replicate_matrices, axis='columns', keys=replicate_keys)
        sort_tup = zip(replicate_keys, ['fold_induction'] * len(replicate_keys))
        protein_matrix.sort_values(by=sort_tup, ascending=[False] * len(replicate_keys), inplace=True)

        with open('{0}.out'.format(protein), 'w') as o:
            o.write('{0}\n'.format(protein_matrix.to_string()))
        sys.stdout.write('{0}\n'.format(protein))
        sys.stdout.write('{0}\n'.format(protein_matrix.info()))
        sys.stdout.write('\n')

        protein_count_matrix = pd.concat(replicate_count_matrices, axis='columns', keys=replicate_keys)
        sort_count = zip(replicate_keys, ['sum'] * len(replicate_keys))
        protein_count_matrix.sort_values(by=sort_count, ascending=[False] * len(replicate_keys), inplace=True)

        with open('{0}_counts.csv'.format(protein), 'w') as o:
            o.write('{0}\n'.format(protein_count_matrix.to_csv()))


def main(fastq_files, barcode_files, library_architecture, max_expected_incorrect, minimum_count):
    protein_replicate_state_dict, barcode_seqs = parse_barcode_file(barcode_files)
    const_library_tup, const_barcode_tup = parse_from_library_architecture(library_architecture)

    barcode_sequence_dict = \
        parse_fastqs_to_dict(fastq_files, max_expected_incorrect, const_library_tup, const_barcode_tup, barcode_seqs)

    calculate_fluorescence_from_seq_barcode_dict(barcode_sequence_dict, protein_replicate_state_dict, minimum_count)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""script to calculate the mean fluorescence value for a specific genotype present in sequence data.
        input is fastq(s) of sequences and a tab delimited text file of barcodes, median fluorescence values,
        population percentages, and replicate numbers for each barcode."""
    )
    required = parser.add_argument_group('required arguments')
    required.add_argument('-f', '--fastq', nargs='*', required=True, help='text file(s) of sequences')
    required.add_argument('-b', '--barcodes', nargs='*', required=True,
                          help='tab delimited text file with one column for barcodes, protein, states, '
                               'median fluorescence values, population percentages, and replicate number.')
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
