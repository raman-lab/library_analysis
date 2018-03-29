#!/usr/bin/env python
import argparse
import pandas as pd
import plotly


def html_scatter_plot_from_csv(input_csv, x_axis, y_axis, text_column, filter_string, x_scale_log, y_scale_log):
    input_df = pd.DataFrame.from_csv(input_csv)
    input_df.dropna(how='any', inplace=True)
    if filter_string:
        column, value = filter_string.split(':')
        input_df = input_df.loc[input_df.iloc[:, int(column)] == float(value)]

    if x_scale_log and y_scale_log:
        fig = {
            'data': [
                {
                    'x': input_df.iloc[:, x_axis],
                    'y': input_df.iloc[:, y_axis],
                    'text': input_df.loc[:, text_column]
                }
            ],
            'layout': {
                'xaxis': {'title': input_df.columns[x_axis], 'type': 'log'},
                'yaxis': {'title': input_df.columns[y_axis], 'type': 'log'}
            }
        }
    elif x_scale_log:
        fig = {
            'data': [
                {
                    'x': input_df.iloc[:, x_axis],
                    'y': input_df.iloc[:, y_axis],
                    'text': input_df.loc[:, text_column]
                }
            ],
            'layout': {
                'xaxis': {'title': input_df.columns[x_axis], 'type': 'log'},
                'yaxis': {'title': input_df.columns[y_axis]}
            }
        }
    elif y_scale_log:
        fig = {
            'data': [
                {
                    'x': input_df.iloc[:, x_axis],
                    'y': input_df.iloc[:, y_axis],
                    'text': input_df.loc[:, text_column]
                }
            ],
            'layout': {
                'xaxis': {'title': input_df.columns[x_axis]},
                'yaxis': {'title': input_df.columns[y_axis], 'type': 'log'}
            }
        }
    else:
        fig = {
            'data': [
                {
                    'x': input_df.iloc[:, x_axis],
                    'y': input_df.iloc[:, y_axis],
                    'text': input_df.loc[:, text_column]
                }
            ],
            'layout': {
                'xaxis': {'title': input_df.columns[x_axis]},
                'yaxis': {'title': input_df.columns[y_axis], 'type': 'log'}
            }
        }
    name = input_csv.split('/')[-1].split('.csv')[0]
    plotly.offline.plot(fig, filename='{0}.html'.format(name), auto_open=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="script generates scatter plot of induction data from "
                                                 "genotype_fluorescence.py")
    parser.add_argument("-x", "--x_axis", default=1, help="zero-indexed column to be plotted on x-axis")
    parser.add_argument("-y", "--y_axis", default=2, help="zero-indexed column to be plotted on y-axis")
    parser.add_argument("-t", "--text_column", default=3, help="zero-indexed column to with text label for data points")
    parser.add_argument("-f", "--filter", help="zero-indexed column:value to select")
    parser.add_argument("--x_scale_log", action='store_true')
    parser.add_argument("--y_scale_log", action='store_true')
    required = parser.add_argument_group('required arguments')
    required.add_argument("-i", "--input_csv", required=True, help="csv file containing data to be plotted")
    args = parser.parse_args()
    html_scatter_plot_from_csv(args.input_csv, args.x_axis, args.y_axis, args.text_column, args.filter,
                               args.x_scale_log, args.y_scale_log)
