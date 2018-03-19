#!/usr/bin/env python
import argparse
import collections
import math


def get_data_by_column(axes, txt_files, seq_list):
    possible_lengths = [16, 17, 18, 19]

    length_tuple_dict = collections.OrderedDict()
    for length in possible_lengths:
        length_tuple_dict[length] = []
    total_matched = 0
    total_seqs = 0
    seqs_with_vals = 0
    seqs_with_vals_matched = 0
    for data_file in txt_files:
        indices = {}
        indices_found = False
        with open(data_file, 'r') as f:
            while not indices_found:
                line = f.readline().rstrip()
                if axes[0] and axes[1] in line:
                    split_line = line.split()
                    x_index = split_line.index(axes[0])
                    y_index = split_line.index(axes[1])
                    indices[axes[0]] = x_index
                    indices[axes[1]] = y_index
                    indices_found = True
            for line in f:
                split_line = line.rstrip().split()
                try:
                    x = float(split_line[indices[axes[0]]])
                    y = float(split_line[indices[axes[1]]])
                    seq = split_line[0]
                    seq_len = len(seq)
                    total_seqs += 1
                    if seq in seq_list:
                        total_matched += 1
                except IndexError:
                    continue
                if x == 0.0 or y == 0.0:
                    continue
                if not math.isnan(x) and not math.isnan(y) and seq_len in possible_lengths:
                    seqs_with_vals += 1
                    length_tuple_dict[seq_len].append((x, y))
                    if seq in seq_list:
                        seqs_with_vals_matched += 1

    print "matched / total: {0} / {1} = {2}".format(float(total_matched), total_seqs, float(total_matched) / total_seqs)
    print "with values: {0} / {1} = {2}".format(
        float(seqs_with_vals_matched), seqs_with_vals, float(seqs_with_vals_matched) / seqs_with_vals)
    return length_tuple_dict


def get_seqs_from_data_file(data_file, allowable_lengths):
    seq_list = []
    with open(data_file, 'r') as f:
        f.next()
        f.next()
        for line in f:
            seq = line.split()[0]
            if len(seq) in allowable_lengths:
                seq_list.append(seq)

    return seq_list


def get_seqs_from_txt_file(txt_files):
    seq_list = []
    for txt in txt_files:
        with open(txt, 'r') as f:
            f.next()
            for line in f:
                seq_list.append(line.split()[0])
    return seq_list


def compare_seqs(txt_files, data_file):
    allowable_lengths = [16, 17, 18, 19]
    data_list = get_seqs_from_data_file(data_file, allowable_lengths)
    seq_list = get_seqs_from_txt_file(txt_files)

    n_to_d = 0.0
    d_to_n = 0.0

    for seq in data_list:
        if seq in seq_list:
            n_to_d += 1
    for seq in seq_list:
        if seq in data_list:
            d_to_n += 1

    print "{0}".format(data_file.split('/')[-1].split('.')[0])
    print "compare nick to devesh: {0} / {1} = {2}".format(n_to_d, len(data_list), n_to_d/len(data_list))
    print "compare devesh to nick: {0} / {1} = {2}".format(d_to_n, len(seq_list), d_to_n/len(seq_list))
    print ""

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="quick script to compare ngs analysis")
    parser.add_argument('-t', nargs='*')
    parser.add_argument('-d')
    args = parser.parse_args()
    compare_seqs(args.t, args.d)
