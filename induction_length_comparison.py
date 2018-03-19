#!/usr/bin/env python
import sys
import pandas as pd
import numpy as np
from induction_plotly import dataframe_from_file
import plotly
import plotly.graph_objs as go
import collections


def compare_lengths(data_file):
    allowable_lengths = [16, 17, 18, 19]
    axes = 'fold_induction'
    dataframe = dataframe_from_file(data_file)
    seq_list = list(dataframe.index)
    seq_length_list = []
    for x in seq_list:
        if type(x) is str:
            seq_length_list.append(len(x))
        else:
            seq_length_list.append(0)
    dataframe['length'] = pd.Series(seq_length_list, index=dataframe.index)
    dataframe = dataframe[dataframe['length'].isin(allowable_lengths)]
    dataframe.dropna(subset=[axes], how='any', inplace=True)
    dataframe_16 = dataframe[dataframe['length'] == 16]
    columns = dataframe.columns.values
    columns = np.append(columns, 'delta')
    protein_dataframe = pd.DataFrame(columns=columns)
    plotting_dict = {
        17: {
            'seq': [],
            'ratio': [],
            'name': []
        },
        18: {
            'seq': [],
            'ratio': [],
            'name': []
        },
        19: {
            'seq': [],
            'ratio': [],
            'name': []
        }
    }

    for seq_16 in dataframe_16.index:
        individual_df = pd.DataFrame(columns=columns)
        individual_df.loc[seq_16] = dataframe_16.loc[seq_16]
        for seq_length in [17, 18, 19]:
            index_seq_list = list(dataframe[dataframe.length == seq_length].index)
            matching_list = [seq for seq in index_seq_list if seq_16 in seq]
            for matching_seq in matching_list:

                individual_df.loc[matching_seq] = dataframe.loc[matching_seq]
                ind_long_seq = individual_df.loc[matching_seq]['fold_induction']
                ind_16 = individual_df.loc[seq_16]['fold_induction']
                individual_df.loc[matching_seq]['delta'] = ind_long_seq / ind_16

                plotting_dict[seq_length]['seq'].append(seq_16)
                plotting_dict[seq_length]['ratio'].append(ind_long_seq / ind_16)
                plotting_dict[seq_length]['name'].append(matching_seq)

        if len(individual_df.index) > 1:
            protein_dataframe = pd.concat([protein_dataframe, individual_df])
    name = '{0}_length_comparison'.format(data_file.split('.')[0])
    with open('{0}.txt'.format(name), 'w') as o:
        o.write('{0}\n'.format(protein_dataframe.to_string()))

    fig = {
        'data': [
            {
                'x': plotting_dict[length]['seq'],
                'y': plotting_dict[length]['ratio'],
                'text': plotting_dict[length]['name'],
                'name': length,
                'mode': 'markers',
            } for length in plotting_dict.keys()
            ],
        'layout': {
            'xaxis': {'title': 'seq'},
            'yaxis': {'title': 'ratio of fold induction'},
            'title': data_file.split('.')[0]}
    }

    plotly.offline.plot(fig, filename='{0}.html'.format(name), auto_open=False)

if __name__ == "__main__":
    data_files = sys.argv[1:]
    for data_file in data_files:
        compare_lengths(data_file)
