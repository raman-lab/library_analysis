#!/usr/bin/env python
import argparse
import itertools
import json
import pandas as pd
import plotly
from induction_plotly import dataframe_from_file


def seq_heatmap_maker(data_file, msa, msa_csi, csi_json, allowable_lengths, column):
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
    dataframe.dropna(subset=[column], how='any', inplace=True)

    with open(csi_json, 'r') as f:
        seq_csi_score_dict = json.load(f)

    bases = ['A', 'C', 'G', 'T']
    fi_sum = 0
    csi_sum = 0
    msa_count = 0

    with open(msa, 'r') as f:
        for tag, alignment in itertools.izip_longest(f, f, fillvalue=None):
            alignment = alignment.rstrip()
            if msa_count == 0:
                length = len(alignment)
                positions = range(1, length + 1)
                position_fi_df = pd.DataFrame(data=0.0, index=bases, columns=positions)
            msa_count += 1
            seq = alignment.replace('-', '')
            data_value = dataframe.loc[seq][column]
            fi_sum += data_value
            for index, base in enumerate(alignment):
                if base is not '-':
                    position_fi_df.loc[base][index + 1] += data_value

    msa_count = 0
    with open(msa_csi, 'r') as f:
        for tag, alignment in itertools.izip_longest(f, f, fillvalue=None):
            alignment = alignment.rstrip()
            if msa_count == 0:
                length = len(alignment)
                positions = range(1, length + 1)
                position_csi_df = pd.DataFrame(data=0.0, index=bases, columns=positions)
            msa_count += 1
            seq = alignment.replace('-', '')
            csi_value = seq_csi_score_dict[seq]
            for index, base in enumerate(alignment):
                if base is not '-':
                    position_csi_df.loc[base][index + 1] += csi_value
            csi_sum += csi_value

    position_fi_df = position_fi_df.divide(fi_sum * length)
    position_csi_df = position_csi_df.divide(csi_sum * length)
    subed_df = position_fi_df.subtract(position_csi_df)

    name = '{0}_{1}'.format(data_file.split('/')[-1].split('.')[0], column)
    with open('{0}.txt'.format(name), 'w') as o:
        o.write('ID Matrix\n')
        o.write('BF\n')
        o.write('PO')
        o.write('{0}'.format(position_fi_df.transpose().to_string()))
    heatmap_data = [plotly.graph_objs.Heatmap(y=list(position_fi_df.index), z=position_fi_df.values.tolist(),
                                              colorscale='Viridis')]
    plotly.offline.plot(heatmap_data, filename='{0}.html'.format(name), auto_open=False)

    name = '{0}_csi'.format(data_file.split('/')[-1].split('.')[0])
    with open('{0}.txt'.format(name), 'w') as o:
        o.write('ID Matrix\n')
        o.write('BF\n')
        o.write('PO')
        o.write('{0}'.format(position_csi_df.transpose().to_string()))
    heatmap_data = [plotly.graph_objs.Heatmap(y=list(position_csi_df.index), z=position_csi_df.values.tolist(),
                                              colorscale='Viridis')]
    plotly.offline.plot(heatmap_data, filename='{0}.html'.format(name), auto_open=False)

    name = '{0}_{1}_sub_csi'.format(data_file.split('/')[-1].split('.')[0], column)
    with open('{0}.txt'.format(name), 'w') as o:
        o.write('ID Matrix\n')
        o.write('BF\n')
        o.write('PO')
        o.write('{0}'.format(subed_df.transpose().to_string()))
    heatmap_data = [plotly.graph_objs.Heatmap(y=list(subed_df.index), z=subed_df.values.tolist(),
                                              colorscale='Viridis')]
    plotly.offline.plot(heatmap_data, filename='{0}.html'.format(name), auto_open=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-l', '--lengths', nargs='*', default=[16, 17, 18, 19],
                        help='one heatmap will be made per length given')
    parser.add_argument('-col', '--column', default='fold_induction', help="column to use in data file")
    required = parser.add_argument_group('required arguments')
    required.add_argument('-d', '--data_file', required=True,
                          help="columns of data separated by white space. columns used to generate plot must be "
                               "named the same as the x and y axes")
    required.add_argument('-m', '--msa', required=True,
                          help='multiple sequence alignment in fasta format of sequences in data file')
    required.add_argument('-mc', '--msa_csi', required=True)
    required.add_argument('-c', '--csi_json', required=True,
                          help='json file of csi sequences and associated values')

    args = parser.parse_args()
    seq_heatmap_maker(args.data_file, args.msa, args.msa_csi, args.csi_json, args.lengths, args.column)