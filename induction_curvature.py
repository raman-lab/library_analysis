#!/usr/bin/env python
import argparse
from dnacurve import CurvedDNA
import pandas as pd
import plotly


def dataframe_from_file(data_file):
    columns = ['repressed', 'induced', 'free', 'fold_induction', 'fold_repression']
    with open(data_file, 'r') as f:
        dataframe = pd.read_table(f, sep='\s+', header=None, skiprows=2, index_col=0)
        dataframe.columns = columns
        dataframe.index.name = 'sequence'
    return dataframe


def induction_curvature(data_file):
    dataframe = dataframe_from_file(data_file)
    dataframe.dropna(subset=['fold_induction', 'repressed'], how='any', inplace=True)
    allowable_lengths = [16, 17, 18, 19]

    plotting_dict = {
        '0-1': {
            'fold induction': [],
            'curvature': [],
            'repression': [],
            'text': []
        },
        '1-1.5': {
            'fold induction': [],
            'curvature': [],
            'repression': [],
            'text': []
        },
        '1.5-2': {
            'fold induction': [],
            'curvature': [],
            'repression': [],
            'text': []
        },
        '2-': {
            'fold induction': [],
            'curvature': [],
            'repression': [],
            'text': []
        }

    }

    for sequence in dataframe.index:
        if sequence is not None and len(sequence) in allowable_lengths:
            seq_curve_obj = CurvedDNA("TTGACA{0}TATAAT".format(sequence), "NUCLEOSOME", name="Example")
            curvature = seq_curve_obj.curvature[0, :]
            curvature = curvature[curvature.nonzero()]
            curvature_sum = curvature.sum()
            fi = dataframe.loc[sequence]['fold_induction']
            rep = dataframe.loc[sequence]['repressed']
            if curvature_sum <= 1:
                curve_bin = '0-1'
            elif curvature_sum <= 1.5:
                curve_bin = '1-1.5'
            elif curvature_sum <= 2:
                curve_bin = '1.5-2'
            else:
                curve_bin = '2-'
            plotting_dict[curve_bin]['fold induction'].append(fi)
            plotting_dict[curve_bin]['repression'].append(rep)
            plotting_dict[curve_bin]['curvature'].append(curvature_sum)
            plotting_dict[curve_bin]['text'].append('{0}<br>{1}'.format(sequence, curvature_sum))

    base_name = '{0}'.format(data_file.split('/')[-1].split('.')[0])
    fig = {
        'data': [
            {
                'x': plotting_dict[curve_bin]['repression'],
                'y': plotting_dict[curve_bin]['fold induction'],
                'text': plotting_dict[curve_bin]['text'],
                'name': curve_bin,
                'mode': 'markers',
            } for curve_bin in plotting_dict.keys()
        ],
        'layout': {
            'xaxis': {'title': 'repression', 'type': 'log'},
            'yaxis': {'title': 'fold induction'},
            'title': base_name}
    }

    plotly.offline.plot(fig, filename='{0}_rep_fi_curve.html'.format(base_name), auto_open=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""script to plot log repression and fold induction from
        sequence-induction datafiles. points are colored by curvature""")
    required = parser.add_argument_group('required')
    required.add_argument('-d', '--data_file', required=True)
    args = parser.parse_args()
    induction_curvature(args.data_file)
