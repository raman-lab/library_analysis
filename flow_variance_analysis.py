#!/usr/bin/env python

# this file is part of the Raman Lab github repository: https://github.com/raman-lab/library_analysis
# author: nwhoppe
# created: 7/5/18

import argparse
import pandas as pd

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
    required.add_argument('-i', '--input_csvs', nargs='*')
    args = parser.parse_args()
    