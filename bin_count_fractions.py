#!/usr/bin/env python
import argparse
import pandas as pd
import sys


def bin_count_fractions(input_file):
    count_df = pd.read_csv(input_file, header=[0, 1, 2], index_col=0)
    state_bin_list = []
    bin_number_dict = {}
    for state in ['repressed', 'induced']:
        state_bin_count = count_df.loc[:, ('replicate 1', state)].count(axis='columns')
        state_bin_count.rename(state, inplace=True)
        state_bin_list.append(state_bin_count)
        bin_number_dict[state] = len(count_df.loc[:, ('replicate 1', state)].columns)

    bin_count_df = pd.concat(state_bin_list, axis='columns')
    bin_count_df = bin_count_df[(bin_count_df > 0).all(axis='columns')]

    shared_bin = pd.DataFrame(data=0, index=range(1, bin_number_dict['repressed'] + 1),
                              columns=range(1, bin_number_dict['induced'] + 1))

    for seq in bin_count_df.index:
        rep_count, ind_count = bin_count_df.loc[seq]
        shared_bin[rep_count][ind_count] += 1

    shared_bin = shared_bin.divide(shared_bin.sum().sum())
    expectation = 0
    variance = 0
    for i in shared_bin.index:
        for j in shared_bin.columns:
            expectation += (i + j) * shared_bin.loc[i][j]
            square_ij = pow((i + j), 2)
            variance += square_ij * shared_bin.loc[i][j]
    variance -= pow(expectation, 2)
    sys.stdout.write('{0}\n'.format(shared_bin.to_string()))
    sys.stdout.write('expected sum: {0}'.format(expectation))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""script to calculate the fraction of sequences that appear in
    multiple bins during sorting. input is csv file from genotype_fluorescence.py""")
    required = parser.add_argument_group('required')
    required.add_argument('-i', '--input_file', required=True)
    args = parser.parse_args()
    bin_count_fractions(args.input_file)
