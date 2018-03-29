#!/usr/bin/env python
import argparse
import pandas as pd
import sys


def line_list_from_txt(txt_file):
    with open(txt_file, 'r') as f:
        line_list = [line.rstrip() for line in f]
    return line_list


def library_coverage_from_input_files(input_csv_files, library_txt, column_index):
    csv_set = set()
    for input_csv in input_csv_files:
        input_df = pd.DataFrame.from_csv(input_csv)
        csv_descriptors = input_df.iloc[:, column_index]
        csv_set.update(csv_descriptors)

    library_descriptors = line_list_from_txt(library_txt)
    library_set = set(library_descriptors)
    intersection = set.intersection(csv_set, library_set)
    library_difference_csv = library_set - csv_set
    coverage = float(len(intersection)) / len(library_set)
    sys.stdout.write('coverage: {0} / {1} = {2}'.format(len(intersection), len(library_set), coverage))
    output_list = ['{0}\n'.format(missing_item) for missing_item in library_difference_csv]
    sys.stdout.writelines(output_list)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""script to compare the number of similar elements between string 
    descriptors in an input csv file to string descriptors in an input text file. the coverage fraction and missing 
    descriptors are printed to stdout""")
    parser.add_argument('-c', '--column_index', type=int, default=3,
                        help='zero-indexed column which contains string descriptors')
    required = parser.add_argument_group('required')
    required.add_argument('-i', '--input_csv_files', required=True, nargs='*',
                          help='csv file containing string descriptors from a library')
    required.add_argument('-l', '--library_txt', required=True,
                          help='text file containing string descriptors from a library. must be one per line')
    args = parser.parse_args()
    library_coverage_from_input_files(args.input_csv_files, args.library_txt, args.column_index)
