#!/usr/bin/env python

# this file is part of the Raman Lab github repository: https://github.com/raman-lab/library_analysis
# author: nwhoppe
# created: 4/24/18

import argparse
import pandas as pd


def trim_csv(csv, output_csv_name, required_column_labels, minimum_index_length, maximum_index_length,
             minimum_column_values, maximum_column_values):
    df = pd.read_csv(csv, index_col=0)
    df = df[df.index.notnull()]
    if minimum_index_length:
        df = df[df.index.map(len) >= minimum_index_length]
    if maximum_index_length:
        df = df[df.index.map(len) <= maximum_index_length]
    if required_column_labels:
        df.dropna(axis='index', how='any', subset=required_column_labels, inplace=True)
    if minimum_column_values:
        for column_value_pair in minimum_column_values:
            column, value = column_value_pair.split(':')
            df = df[df[column] >= int(value)]
    if maximum_column_values:
        for column_value_pair in maximum_column_values:
            column, value = column_value_pair.split(':')
            df = df[df[column] >= int(value)]

    df.to_csv(output_csv_name)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""script to edit input csv files according to given arguments and 
    output new csv file. more functionality will be added as needed""")
    parser.add_argument('-min_il', '--minimum_index_length', type=int, help='minimum length each index string must be')
    parser.add_argument('-max_il', '--maximum_index_length', type=int, help='maximum length each index string must be')
    parser.add_argument('-r', '--required_columns', nargs='*',
                        help='space separated list of column labels. for all column labels given,'
                             'rows will be required to have a non NA value present to be output')
    parser.add_argument('-min_cvals', '--minimum_column_values', nargs='*',
                        help='space separated pairs of <column_index_label>:<min_integer_value>')
    parser.add_argument('-max_cvals', '--maximum_column_values', nargs='*',
                        help='space separated pairs of <column_index_label>:<max_integer_value>')
    required = parser.add_argument_group('required')
    required.add_argument('-c', '--csv', required=True)
    required.add_argument('-n', '--name', required=True, help='name of output csv')
    args = parser.parse_args()
    trim_csv(args.csv, args.name, args.required_columns, args.minimum_index_length, args.maximum_index_length,
             args.minimum_column_values, args.maximum_column_values)

