#!/usr/bin/env python
import argparse
import pandas as pd


def column_ratio_from_csv(input_csv_list, name, denominator, numerator):
    if len(input_csv_list) == 1:
        input_csv = input_csv_list[0]
        input_df = pd.DataFrame.from_csv(input_csv)
        input_df[name] = input_df[numerator] / input_df[denominator]
        input_df.to_csv(input_csv)
    elif len(input_csv_list) == 2:
        input_df1 = pd.DataFrame.from_csv(input_csv_list[0])
        input_df2 = pd.DataFrame.from_csv(input_csv_list[1])
        column1 = '{0}_{1}'.format(input_csv_list[0].split('/')[-1].split('.csv')[0], numerator)
        column2 = '{0}_{1}'.format(input_csv_list[0].split('/')[-1].split('.csv')[0], denominator)
        input_df1.rename(index=str, columns={numerator: column1}, inplace=True)
        input_df2.rename(index=str, columns={denominator: column2}, inplace=True)
        output_df = pd.concat([input_df1[column1], input_df2[column2]], axis='columns')
        output_df['ratio'] = input_df1[column1] / input_df2[column2]
        output_df.dropna(how='any', inplace=True)
        output_df.to_csv('{0}'.format(name))
    else:
        raise Exception('script only accepts one or two input files')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""script to calculate ratio of two columns from input csv file(s).
    if one csv is given, the column is appended. if two csvs are given, a new file is written""")
    parser.add_argument('-den', '--denominator', default='repressed', help='name of denominator column')
    parser.add_argument('-num', '--numerator', default='induced', help='name of numerator column')
    parser.add_argument('-n', '--name', default='fi',
                        help='name of ratio column if one file given, name of output file if two files given')
    required = parser.add_argument_group('required')
    required.add_argument('-i', '--input_csv', nargs='*', required=True,
                          help='if one csv is given, the column is appended. '
                               'if two csvs are given, a new file is written. first csv is used as numerator.'
                               'second csv is used as denominator')
    args = parser.parse_args()
    column_ratio_from_csv(args.input_csv, args.name, args.denominator, args.numerator)
