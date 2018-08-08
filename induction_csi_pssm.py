#!/usr/bin/env python
import argparse
import numpy as np
import pandas as pd
import scipy.interpolate
import subprocess
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns


def dataframe_from_file(data_file):
    columns = ['repressed', 'induced', 'free', 'fold_induction', 'fold_repression']
    with open(data_file, 'r') as f:
        dataframe = pd.read_table(f, sep='\s+', header=None, skiprows=2, index_col=0)
        dataframe.columns = columns
        dataframe.index.name = 'sequence'
    return dataframe


def contour_plot(xy_dataframe):
    almost_gray = '#808080'
    sns_plot = sns.jointplot(x='pssm score', y='fold induction', data=xy_dataframe, kind='kde', n_levels=25)
    sns_plot.plot_joint(plt.scatter, c=almost_gray, s=5, marker='o')
    sns_plot.savefig("sns_test.png")


def score_seq_by_pssm(pssm, aligned_seq, non_gap_positions):
    score = 0.0
    for position in non_gap_positions:
        base = aligned_seq[position - 1]
        pssm_value = pssm.loc[base][position]
        score += pssm_value
    return score


def progress(count, total, suffix=''):
    bar_len = 60
    filled_len = int(round(bar_len * count / float(total)))
    percents = round(100 * count / float(total), 1)
    bar = '=' * filled_len + '-' * (bar_len - filled_len)
    sys.stdout.write('[%s] %s%s ...%s\r' % (bar, percents, '%', suffix))
    sys.stdout.flush()


def scatter_plot_with_percentile(xy_dataframe, percentile, name):
    almost_gray = '#808080'
    almost_black = '#262626'
    color_set = {
        'CmeR': '#ab62c0',
        'DesT': '#50ab6d',
        'NalC': '#c95573',
        'PmeR': '#929d3d',
        'SmeT': '#648ace',
        'TtgR': '#ca743e'
    }
    # For bold font:
    plt.rcParams["font.weight"] = "bold"
    plt.rcParams["axes.labelweight"] = "bold"
    # note original tick label size = 12, fontsize = 15

    xy_dataframe = xy_dataframe[xy_dataframe['pssm score'] >= 0]
    # xy_dataframe = xy_dataframe[xy_dataframe['pssm score'] >= -40]
    x = np.asarray(xy_dataframe['pssm score'])
    y = np.asarray(xy_dataframe['fold induction'])
    q90 = np.percentile(y, 90)
    x_below_q90 = x[np.where(y < q90)]
    y_below_q90 = y[np.where(y < q90)]
    x_above_q90 = x[np.where(y >= q90)]
    y_above_q90 = y[np.where(y >= q90)]

    # x_bins = np.linspace(min(x), max(x), 16)
    # x_bins = range(int(min(x)) / 5, int(max(x)) + 5, 5)
    # bin_indices = np.digitize(x, x_bins)
    # y_bin_avg = []
    # y_bin_q75 = []
    # for i in range(1, len(x_bins)):
    #     if i in bin_indices:
    #         y_bin = y[np.where(bin_indices == i)]
    #         q25, q75, q90 = np.percentile(y, [25, 75, 90])
    #         iqr = q75 - q25
    #         avg_fi = np.average(y_bin[np.where(y_bin >= q75 + 1.5 * iqr)])
    #         avg_fi_gt_q75 = np.average(y_bin[np.where(y_bin >= q90)])
    #         # y_above_percentile = [k for k in y_bin if k >= q]
    #         # y_above_percentile_median = np.median(y_above_percentile)
    #         # y_percentile_medians.append(y_above_percentile_median)
    #         y_bin_avg.append(avg_fi)
    #         y_bin_q75.append(avg_fi_gt_q75)
    #     else:
    #         x_bins = np.delete(x_bins, i - 1)
    #
    # x_avg_bins = []
    # for i in range(0, len(x_bins) - 1):
    #     avg = (x_bins[i] + x_bins[i + 1]) / 2.0
    #     x_avg_bins.append(avg)

    # x_smooth = np.linspace(min(x_avg_bins), max(x_avg_bins), 500)
    # cs = scipy.interpolate.CubicSpline(x_avg_bins, y_bin_avg)
    # cs2 = scipy.interpolate.CubicSpline(x_avg_bins, y_bin_q75)
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    # ax1.scatter(x, y, c=almost_gray, marker='o', alpha=0.5, edgecolor=almost_black, linewidth=0.15)
    # ax1.plot(x_smooth, cs(x_smooth), color=almost_black, alpha=0.75)
    # ax1.plot(x_smooth, cs2(x_smooth), color='blue', alpha=0.75)
    ax1.scatter(x_below_q90, y_below_q90, c=almost_gray, marker='o', alpha=0.5, edgecolor=almost_black, linewidth=0.15)
    ax1.scatter(x_above_q90, y_above_q90, c=color_set[name.split('_')[0]], marker='o', alpha=0.75,
                edgecolor=almost_black, linewidth=0.15)

    ax1.set_ybound(lower=0.0)
    # ax1.set_xbound(lower=min(x) - 1, upper=max(x) + 1)
    ax1.tick_params(labelsize=13)
    plt.xlabel('Pssm Score', fontsize=16)
    plt.ylabel('Fold Induction', fontsize=16)
    plt.text(max(x), max(y), s=name.split('_')[0], fontsize=16, horizontalalignment='right')

    spines_to_remove = ['top', 'right']
    for spine in spines_to_remove:
        ax1.spines[spine].set_visible(False)
    ax1.xaxis.set_ticks_position('none')
    ax1.yaxis.set_ticks_position('none')
    spines_to_keep = ['bottom', 'left']
    for spine in spines_to_keep:
        ax1.spines[spine].set_linewidth(0.5)
        ax1.spines[spine].set_color(almost_black)
    ax1.xaxis.label.set_color(almost_black)
    ax1.yaxis.label.set_color(almost_black)

    plt.savefig('{0}.pdf'.format(name))
    plt.close()


def violin_plot(xy_dataframe, basename):
    pd.options.mode.chained_assignment = None
    xy_dataframe = xy_dataframe[xy_dataframe['pssm score'] >= 0]
    # xy_dataframe = xy_dataframe[xy_dataframe['pssm score'] >= -40]
    min_x = int(min(xy_dataframe['pssm score'])) / 5
    max_x = int(max(xy_dataframe['pssm score'])) + 5
    bins = range(min_x, max_x, 5)
    labels = []
    for i in range(0, len(bins) - 1):
        label = '{0}-{1}'.format(bins[i], bins[i + 1])
        labels.append(label)
    # print pd.cut(xy_dataframe.loc[:, 'pssm score'], bins=bins)
    xy_dataframe['bins'] = pd.cut(xy_dataframe['pssm score'], bins=bins, labels=labels)
    sns.set_style("whitegrid")
    ax = sns.violinplot(x="bins", y="fold induction", data=xy_dataframe[['bins', 'fold induction']],
                        cut=0, scale='width', palette='muted')
    fig = ax.get_figure()
    fig.savefig('{0}_violin.pdf'.format(basename))


def induction_csi_psiblast(dataframe_file, pssm_json_file, profile, clustal_path, filter_bool, output_json, basename):
    percentile = 90
    if not basename:
        basename = profile.split('/')[-1].split('.')[0]

    if not output_json:
        allowable_lengths = [16, 17, 18, 19]
        input_df = dataframe_from_file(dataframe_file)
        input_df = input_df[input_df.fold_induction > 0]
        pssm = pd.read_json(pssm_json_file)

        if filter_bool:
            input_df = input_df[(input_df.repressed < input_df.repressed.quantile(0.8)) & (input_df.fold_induction > 2)]
        output_df = pd.DataFrame(columns=['pssm score', 'fold induction'])

        total = len(input_df.index)

        for s, sequence in enumerate(input_df.index):
            progress(s, total, '')
            if sequence is not None and len(sequence) in allowable_lengths:
                # top_psiblast_hit, e_val = run_clustal_from_shell(clustal_path, sequence, e_val=0.01)
                # if e_val <= 0.01:
                #     csi_value = seq_csi_score_dict[top_psiblast_hit]
                #     fold_ind = input_df.loc[sequence]['fold_induction']
                #     output_df.loc[sequence] = {'csi score': np.log10(csi_value), 'fold induction': fold_ind}
                with open('seq.fasta', 'w') as o:
                    o.write('>seq\n')
                    o.write('{0}\n'.format(sequence))

                clustal_command = [clustal_path, '--profile1', 'seq.fasta', '--profile2', '{0}'.format(profile)]
                p = subprocess.Popen(clustal_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                     stdin=subprocess.PIPE)
                output, error = p.communicate()
                aligned_seq = output.split('\n')[1]
                non_gap_positions = [i + 1 for i, pos in enumerate(aligned_seq)
                                     if (pos != '-' and i + 1 in pssm.columns)]

                aligned_seq_score = score_seq_by_pssm(pssm, aligned_seq, non_gap_positions)
                fi = input_df.loc[sequence]['fold_induction']
                output_df.loc[sequence] = {'pssm score': aligned_seq_score, 'fold induction': fi}
        sys.stdout.write('\n')
        output_df.to_json('{0}.json'.format(basename))
    else:
        output_df = pd.read_json(output_json)
    scatter_plot_with_percentile(output_df, percentile, basename)
    # violin_plot(output_df, basename)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""script to compare dna sequences from induction dataframe to csi
    enrichment sequences""")
    parser.add_argument('-o', '--output_json', help='json file with pssm scores snd fold induction values from '
                                                    'previous run. if provided, script jumps to plotting.')
    parser.add_argument('-n', '--name_base', help='prefix for file naming')

    parser.add_argument('-e', '--executable', default='/Users/nwhoppe/raman_lab/scripts/clustalo',
                        help='path to clustal command line executable.')
    parser.add_argument('-f', '--filter', action='store_true',
                        help='enforce filter on induction sequences used. the filter requires sequences have a fold '
                             'induction greater than 2 and a repressed fluorescence value less than the 80th '
                             'percentile. in the future, the values used by the filter may be flexible')
    required = parser.add_argument_group('required')
    required.add_argument('-d', '--data', required=True, help='induction dataframe from genotype_fluorescence.py')
    required.add_argument('-m', '--pssm_json', required=True,
                          help='json file containing of pssm from msa')
    required.add_argument('-p', '--profile', required=True, help='alignment of csi sequences for use in clustal')
    args = parser.parse_args()
    induction_csi_psiblast(args.data, args.pssm_json, args.profile, args.executable, args.filter, args.output_json,
                           args.name_base)
