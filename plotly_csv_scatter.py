#!/usr/bin/env python
import argparse
import pandas as pd
import plotly


def html_scatter_plot_from_csv(input_csv_files, x_axis, y_axis, text_column, filter_string, x_scale_log, y_scale_log,
                               name):
    df_dict = {}
    for input_csv in input_csv_files:
        input_df = pd.DataFrame.from_csv(input_csv)
        input_df.dropna(how='any', inplace=True)
        if filter_string:
            column, value = filter_string.split(':')
            input_df = input_df.loc[input_df.iloc[:, int(column)] == float(value)]
        df_dict[input_csv.split('/')[-1].split('.csv')[0]] = input_df
    keys = df_dict.keys()
    if x_scale_log and y_scale_log:
        fig = {
            'data': [
                {
                    'x': df_dict[key].iloc[:, x_axis],
                    'y': df_dict[key].iloc[:, y_axis],
                    'text': df_dict[key].iloc[:, text_column],
                    'mode': 'markers',
                    'name': key
                } for key in df_dict.keys()
            ],
            'layout': {
                'xaxis': {'title': df_dict[keys[0]].columns[x_axis], 'type': 'log'},
                'yaxis': {'title': df_dict[keys[0]].columns[y_axis], 'type': 'log'}
            }
        }
    elif x_scale_log:
        fig = {
            'data': [
                {
                    'x': df_dict[key].iloc[:, x_axis],
                    'y': df_dict[key].iloc[:, y_axis],
                    'text': df_dict[key].iloc[:, text_column],
                    'mode': 'markers',
                    'name': key
                } for key in df_dict.keys()
            ],
            'layout': {
                'xaxis': {'title': df_dict[keys[0]].columns[x_axis], 'type': 'log'},
                'yaxis': {'title': df_dict[keys[0]].columns[y_axis], 'type': 'log'}
            }
        }
    elif y_scale_log:
        fig = {
            'data': [
                {
                    'x': df_dict[key].iloc[:, x_axis],
                    'y': df_dict[key].iloc[:, y_axis],
                    'text': df_dict[key].iloc[:, text_column],
                    'mode': 'markers',
                    'name': key
                } for key in df_dict.keys()
            ],
            'layout': {
                'xaxis': {'title': df_dict[keys[0]].columns[x_axis], 'type': 'log'},
                'yaxis': {'title': df_dict[keys[0]].columns[y_axis], 'type': 'log'}
            }
        }
    else:
        fig = {
            'data': [
                {
                    'x': df_dict[key].iloc[:, x_axis],
                    'y': df_dict[key].iloc[:, y_axis],
                    'text': df_dict[key].iloc[:, text_column],
                    'mode': 'markers',
                    'name': key
                } for key in df_dict.keys()
            ],
            'layout': {
                'xaxis': {'title': df_dict[keys[0]].columns[x_axis], 'type': 'log'},
                'yaxis': {'title': df_dict[keys[0]].columns[y_axis], 'type': 'log'}
            }
        }
    if not name:
        name = '{0}.html'.format(input_csv_files[0].split('/')[-1].split('.csv')[0])
    plotly.offline.plot(fig, filename='{0}'.format(name), auto_open=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="script generates scatter plot of induction data from "
                                                 "genotype_fluorescence.py")
    parser.add_argument("-x", "--x_axis", default=0, help="zero-indexed column to be plotted on x-axis")
    parser.add_argument("-y", "--y_axis", default=1, help="zero-indexed column to be plotted on y-axis")
    parser.add_argument("-t", "--text_column", default=3, help="zero-indexed column to with text label for data points")
    parser.add_argument("-f", "--filter", help="zero-indexed column:value to select")
    parser.add_argument("-xl", "--x_scale_log", action='store_true')
    parser.add_argument("-yl", "--y_scale_log", action='store_true')
    parser.add_argument("-n", "--name")
    required = parser.add_argument_group('required arguments')
    required.add_argument("-i", "--input_csv_files", required=True, nargs='*',
                          help="csv file containing data to be plotted")
    args = parser.parse_args()
    html_scatter_plot_from_csv(args.input_csv_files, args.x_axis, args.y_axis, args.text_column, args.filter,
                               args.x_scale_log, args.y_scale_log, args.name)
