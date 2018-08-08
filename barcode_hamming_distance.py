#!/usr/bin/env python

# this file is part of the Raman Lab github repository: https://github.com/raman-lab/library_analysis
# author: nwhoppe
# created: 7/10/18

import argparse
import collections
import itertools
import numpy as np
import pandas as pd
import random
import sys


def compute_hamming_distance(string1, string2):
    """Return the Hamming distance between equal-length sequences"""
    if len(string1) != len(string2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(el1 != el2 for el1, el2 in zip(string1, string2))


def parse_barcode_output(barcode_pair_file):
    min_hamming_distance = False
    barcode_pairs = []
    with open(barcode_pair_file, 'r') as f:
        for line in f:
            if line.startswith('Largest minimum hamming distance:'):
                min_hamming_distance = int(line.rstrip().split()[-1])
            if min_hamming_distance and line[0].isdigit():
                line_list = line.rstrip().split()
                barcode_str = '-'.join([line_list[1], line_list[3]])
                barcode_pairs.append(barcode_str)
    return barcode_pairs, min_hamming_distance


def barcode_hamming_distance(input_barcode_csv, number_of_pairs, iterations, barcode_pair_file):
    barcode_df = pd.read_csv(input_barcode_csv, index_col=False)
    nucleotide_set = set('ACGTacgt')
    barcode_list = []
    for barcode in barcode_df.iloc[:, 0].values:
        if set(barcode) <= nucleotide_set:
            barcode_list.append(barcode)
        else:
            raise Exception('{0} is not a dna sequence'.format(barcode))

    barcode_pair_str_list = []
    for barcode1, barcode2 in itertools.product(barcode_list, repeat=2):
        barcode_pair_str_list.append('-'.join([barcode1, barcode2]))

    barcode_pair_distance_dict = {}
    sys.stdout.write('Computing hamming distance for all pairs\n')
    for barcode_pair1, barcode_pair2 in itertools.combinations_with_replacement(barcode_pair_str_list, 2):
        if barcode_pair1 == barcode_pair2:
            barcode_pair_distance_dict[(barcode_pair1, barcode_pair2)] = 0
        else:
            hamming_distance = compute_hamming_distance(barcode_pair1, barcode_pair2)
            barcode_pair_distance_dict[(barcode_pair1, barcode_pair2)] = hamming_distance
            barcode_pair_distance_dict[(barcode_pair2, barcode_pair1)] = hamming_distance

    if not barcode_pair_file:
        # need to consider all pairwise distances in group of N
        largest_average_barcode_pairs = []
        largest_average_hamming_distance = -1
        min_for_largest_average = 0

        largest_minimum_barcode_pairs = []
        largest_minimum_hamming_distance = -1
        avg_for_largest_minimum = 0

        # previously_tried_sets = []
        sys.stdout.write('Sampling random barcode combinations\n')
        i = 0
        while i < iterations:
            random.shuffle(barcode_pair_str_list)
            random_barcode_pairs = barcode_pair_str_list[0:number_of_pairs]

            hamming_distances = []
            for barcode_pair1, barcode_pair2 in itertools.combinations(random_barcode_pairs, 2):
                hamming_distances.append(barcode_pair_distance_dict[barcode_pair1, barcode_pair2])
            average_hamming_distance = np.average(hamming_distances)
            minimum_hamming_distance = min(hamming_distances)

            if average_hamming_distance > largest_average_hamming_distance:
                largest_average_hamming_distance = average_hamming_distance
                largest_average_barcode_pairs = random_barcode_pairs
                min_for_largest_average = minimum_hamming_distance

            if minimum_hamming_distance > largest_minimum_hamming_distance:
                largest_minimum_hamming_distance = minimum_hamming_distance
                largest_minimum_barcode_pairs = random_barcode_pairs
                avg_for_largest_minimum = average_hamming_distance

            if i % 100000 == 0:
                sys.stdout.write('Iteration: {0}\n'.format(i))

            # previously_tried_sets.append(random_barcode_pairs)
            i += 1

        # output largest average
        sys.stdout.write('\nLargest average hamming distance: {0}\n'.format(largest_average_hamming_distance))
        sys.stdout.write('Minimum hamming distance: {0}\n'.format(min_for_largest_average))
        sys.stdout.write('Corresponding barcodes (number and sequence):\n')
        for barcode_pair_str in largest_average_barcode_pairs:
            barcode1, barcode2 = barcode_pair_str.split('-')
            index1 = barcode_list.index(barcode1) + 1
            index2 = barcode_list.index(barcode2) + 1
            sys.stdout.write('{0}\t{1}\t{2}\t{3}\n'.format(index1, barcode1, index2, barcode2))

        # output largest minimum
        sys.stdout.write('\nLargest minimum hamming distance: {0}\n'.format(largest_minimum_hamming_distance))
        sys.stdout.write('Average hamming distance: {0}\n'.format(avg_for_largest_minimum))
        sys.stdout.write('Corresponding barcodes (number and sequence):\n')
        for barcode_pair_str in largest_minimum_barcode_pairs:
            barcode1, barcode2 = barcode_pair_str.split('-')
            index1 = barcode_list.index(barcode1) + 1
            index2 = barcode_list.index(barcode2) + 1
            sys.stdout.write('{0}\t{1}\t{2}\t{3}\n'.format(index1, barcode1, index2, barcode2))

    else:
        previous_pair_str_list, previous_min_hamming_distance = parse_barcode_output(barcode_pair_file)

        sys.stdout.write('Sampling new barcode pair additions\n')
        i = 0
        while i < iterations:
            random.shuffle(barcode_pair_str_list)
            random_barcode_pairs = barcode_pair_str_list[0:number_of_pairs]

            hamming_distances = []
            for barcode_pair1, barcode_pair2 in itertools.combinations(random_barcode_pairs, 2):
                hamming_distances.append(barcode_pair_distance_dict[barcode_pair1, barcode_pair2])
            minimum_hamming_distance = min(hamming_distances)

            if minimum_hamming_distance >= previous_min_hamming_distance:
                new_barcode_pairs = []
                for new_pair in random_barcode_pairs:
                    take_new_pair = True
                    for previous_pair in previous_pair_str_list:
                        if barcode_pair_distance_dict[new_pair, previous_pair] < previous_min_hamming_distance:
                            take_new_pair = False
                            break
                    if take_new_pair:
                        new_barcode_pairs.append(new_pair)

                if len(new_barcode_pairs) == number_of_pairs:
                    break

            if i % 100000 == 0:
                sys.stdout.write('Iteration: {0}\n'.format(i))
            i += 1

        if i < iterations:
            sys.stdout.write('\nAdditional barcodes with minimum hamming distance {0}\n'.format(
                previous_min_hamming_distance))
            for new_barcode_pair_str in new_barcode_pairs:
                barcode1, barcode2 = new_barcode_pair_str.split('-')
                index1 = barcode_list.index(barcode1) + 1
                index2 = barcode_list.index(barcode2) + 1
                sys.stdout.write('{0}\t{1}\t{2}\t{3}\n'.format(index1, barcode1, index2, barcode2))

        else:
            sys.stdout.write('\nMax iterations reached\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""script that writes pairs of barcodes with increasing hamming 
    distance to stdout given an input csv file of barcodes. first column of csv must contain the barcode sequences""")
    parser.add_argument('-t', '--iterations', type=int, default=1e6)
    parser.add_argument('-n', '--number_of_pairs', type=int, default=87)
    parser.add_argument('-p', '--pairs', required=True,
                        help='input txt file of containing barcode pairs output by barcode_hamming_distance.py ')
    required = parser.add_argument_group('required')
    required.add_argument('-i', '--input_barcode_csv', required=True,
                          help="input csv with first column containing barcode sequences in 5'-3' direction")
    args = parser.parse_args()
    barcode_hamming_distance(args.input_barcode_csv, args.number_of_pairs, args.iterations, args.pairs)
