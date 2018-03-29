#!/usr/bin/env python
import argparse
import pandas as pd
import plotly
import plotly.graph_objs as go


def dms_heatmap(input_csv_files, xaxis_range, name):
    aa_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    df_list = []
    for input_csv in input_csv_files:
        input_df = pd.DataFrame.from_csv(input_csv)
        input_df.dropna(how='any', inplace=True)
        point_mut_df = input_df.loc[input_df['Number of Mutations'] == 1]
        df_list.append(point_mut_df)
    dms_df = pd.concat(df_list)
    start_pos, end_pos = map(int, xaxis_range.split('-'))
    fi_df = pd.DataFrame(float('NaN'), index=aa_list, columns=range(start_pos, end_pos + 1))

    for index in dms_df.index:
        mutation_str = dms_df.loc[index, 'Mutation']
        repressed = dms_df.loc[index, 'repressed']
        induced = dms_df.loc[index, 'induced']
        position = int(mutation_str[1:-1])
        aa = mutation_str[-1]
        fi_df.loc[aa, position] = float(induced) / repressed

    fi_df.dropna(axis='columns', how='all', inplace=True)
    trace = go.Heatmap(z=fi_df.values.tolist(), x=fi_df.columns.tolist(), y=fi_df.index.tolist(), colorscale='Viridis')
    data = [trace]
    if not name:
        name = '{0}_heatmap.html'.format(input_csv_files[0].split('/')[-1].split('.csv')[0])
    plotly.offline.plot(data, filename='{0}'.format(name), auto_open=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""script to make heatmap from dms data within input csv files.
    script is sensitive to formatting of csv files. csv file must have the following columns: repressed, induced, 
    Mutation, Number of Mutations. heatmaps will be colored by the ration of induced/repressed values""")

    parser.add_argument('-x', '--xaxis_range', default='1-300',
                        help='positions to be included on the x axis. must be two integers separated by hyphen. '
                             'positions without any values will be dropped')
    parser.add_argument('-n', '--name')
    required = parser.add_argument_group('required')
    required.add_argument('-i', '--input_csv_files', required=True, nargs='*',
                          help='csv file containing data to be plotted')
    args = parser.parse_args()
    dms_heatmap(args.input_csv_files, args.xaxis_range, args.name)