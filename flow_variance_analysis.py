#!/usr/bin/env python

# this file is part of the Raman Lab github repository: https://github.com/raman-lab/library_analysis
# author: nwhoppe
# created: 7/5/18

import argparse
import collections
import itertools
import numpy as np
import pandas as pd
from scipy import stats
import sys


def flow_variance_analysis(input_csvs, column, significance):
    column_data_dict = collections.OrderedDict()
    variances = []
    max_p = 0
    for csv in input_csvs:
        df = pd.read_csv(csv, index_col=False)
        df = df[df[column] > 0]
        col_list = df[column].tolist()

        # natural log transformation
        log_transformed_list = np.log(col_list)
        variances.append(np.var(log_transformed_list))

        # test normality of data
        skew_kurtosis_stat, p_norm = stats.normaltest(log_transformed_list)
        if p_norm > max_p:
            max_p = p_norm
        sys.stdout.write('Normality test for {0} gives p = {1}\n'.format(csv, p_norm))
        column_data_dict[csv] = log_transformed_list

    variance_ratio = max(variances) / min(variances)
    sys.stdout.write('\nRatio of max and min variance: ')
    sys.stdout.write('{0} / {1} = {0}\n\n'.format(max(variances), min(variances), variance_ratio))
    kruskal_wallis = False

    if max_p > significance and variance_ratio > 3:
        sys.stdout.write('Normality and variance assumption failed. Kruskal-Wallis will be run\n')
        kruskal_wallis = True
        post_hoc = 'Wilcoxon rank sum'
    elif max_p > significance and variance_ratio <= 3:
        sys.stdout.write('Normality assumption failed. Kruskal-Wallis will be run\n')
        kruskal_wallis = True
        post_hoc = 'Wilcoxon rank sum'
    elif max_p < significance and variance_ratio > 3:
        sys.stdout.write('Variance assumption failed. Kruskal-Wallis will be run.\n')
        kruskal_wallis = True
        post_hoc = 't-test with unequal variance'
    else:
        sys.stdout.write('Normality and variance assumption satisfied. ANOVA will be run\n')
        post_hoc = 't-test with equal variance'
    if kruskal_wallis:
        test_stat, p_test = stats.kruskal(*column_data_dict.values())
    else:
        test_stat, p_test = stats.f_oneway(*column_data_dict.values())

    sys.stdout.write('\nStatistical test gives p = {0}\n'.format(p_test))
    if p_test <= significance:
        sys.stdout.write('Significant difference between samples\n\n')
        sys.stdout.write('Post hoc tests: {0}\n'.format(post_hoc))

        for file1, file2 in itertools.combinations(column_data_dict.keys(), 2):
            if post_hoc == 't-test with equal variance':
                test_stat, p_post_hoc = stats.ttest_ind(column_data_dict[file1], column_data_dict[file2],
                                                        equal_var=True)
            elif post_hoc == 't-test with unequal variance':
                test_stat, p_post_hoc = stats.ttest_ind(column_data_dict[file1], column_data_dict[file2],
                                                        equal_var=False)
            else:
                test_stat, p_post_hoc = stats.ranksums(column_data_dict[file1], column_data_dict[file2])

            if p_post_hoc <= significance:
                sys.stdout.write('Significant difference between {0} and {1} (p = {2})\n'.format(
                    file1, file2, p_post_hoc))
    else:
        sys.stdout.write('NO significant difference between samples\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""script to run an anova on input csv files. csv files should have 
    a column of gfp events titled: ''
    values will be log transformed and checked to see if anova assumptions are satisfied:
        Each sample is from a normally distributed population
        The population standard deviations of the groups are all equal (the rule of thumb that 
            max(variance) / min(variance) < 3 will be used)
    it is assumed input files are independent.
    if assumptions are not satisfied, a Kruskal-Wallis test will be run.
    if assumptions are satisfied, pairwise t-tests will be done.
    """)
    required = parser.add_argument_group('required')
    parser.add_argument('-c', '--column', default='EGFP-A-Compensated',
                        help='header of column containing data to be analyzed')
    parser.add_argument('-s', '--significance', default=0.05, type=float,
                        help='max p value to reject null and conclude there is a significant difference')
    required.add_argument('-i', '--input_csvs', nargs='*')
    args = parser.parse_args()
    flow_variance_analysis(args.input_csvs, args.column, args.significance)
