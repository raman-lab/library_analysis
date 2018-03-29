#!/usr/bin/env python
import argparse
import pandas as pd
import sys
from plotly_csv_scatter import dataframe_from_file


def split_dataframe(input_file, fold_induction_range, repressed_value):
    allowable_lengths = [16, 17, 18, 19]
    fi_min, fi_max = fold_induction_range.split('-')
    input_df = dataframe_from_file(input_file)

    seq_list = list(input_df.index)
    seq_length_list = []
    for x in seq_list:
        if type(x) is str:
            seq_length_list.append(len(x))
        else:
            seq_length_list.append(0)
    input_df['length'] = pd.Series(seq_length_list, index=input_df.index)
    input_df = input_df[input_df['length'].isin(allowable_lengths)]

    if repressed_value:
        if not fi_max:
            output_df = input_df[
                (input_df.repressed < repressed_value) &
                (input_df.fold_induction > float(fi_min))
                ]
            output_name = '{0}_repressed_fi_min_cutoff.out'.format(input_file.split('/')[-1].split('.')[0])
        else:
            output_df = input_df[
                (input_df.repressed < repressed_value) &
                (input_df.fold_induction > float(fi_min)) &
                (input_df.fold_induction < float(fi_max))
            ]
            output_name = '{0}_repressed_fi_range_cutoff.out'.format(input_file.split('/')[-1].split('.')[0])

    else:
        if not fi_max:
            output_df = input_df[(input_df.fold_induction > float(fi_min))]
            output_name = '{0}_fi_min_cutoff.out'.format(input_file.split('/')[-1].split('.')[0])
        else:
            output_df = input_df[
                (input_df.fold_induction > float(fi_min)) &
                (input_df.fold_induction < float(fi_max))
            ]
            output_name = '{0}_fi_range_cutoff.out'.format(input_file.split('/')[-1].split('.')[0])

    with open(output_name, 'w') as o:
        o.write('{0}\n'.format(output_df.to_string()))

    for s, seq in enumerate(output_df.index):
        sys.stdout.write('>seq_{0}\n'.format(s))
        sys.stdout.write('{0}\n'.format(seq))


if __name__ == '__main__':
    parser = argparse.ArgumentParser("""script to split up dataframe into two files from genotype_fluorescence.py
    based on fold_induction values and repression percentiles""")
    required = parser.add_argument_group('required')
    required.add_argument('-i', '--input_file', required=True, help='file from genotype_fluorescence.py')
    required.add_argument('-fi', '--fold_induction_range', required=True, type=str,
                          help='range of fold induction values to output. put hyphen between numbers. if no upper '
                               'limit, do not put number after hyphen')
    parser.add_argument('-r', '--repressed', type=float,
                        help='repressed value that sequence must be under to be output')
    args = parser.parse_args()
    split_dataframe(args.input_file, args.fold_induction_range, args.repressed)
