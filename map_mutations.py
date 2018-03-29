#!/usr/bin/env python
import argparse
import pandas as pd
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


def compare_sequence_to_reference(sequence, reference, position_list=False):
    mutation_string_list = []
    if len(sequence) != len(reference):
        raise Exception('sequence and reference must be the same length')
    if not position_list:
        position_list = range(1, len(reference) + 1)
    for index, ref_char in enumerate(reference):
        if sequence[index] != ref_char:
            position = position_list[index]
            mutation_string = '{0}{1}{2}'.format(ref_char, position, sequence[index])
            mutation_string_list.append(mutation_string)
    return '_'.join(mutation_string_list)


def map_mutations(input_csv, wild_type_txt, position_range, leading_to_drop, trailing_to_drop):
    with open(wild_type_txt, 'r') as f:
        wild_type_dna = f.readline().rstrip().upper()
    if len(wild_type_dna) % 3 != 0:
        raise Exception('wild type dna sequence must be divisible 3. it cannot have a partial codon')
    wild_type_aa = str(Seq(wild_type_dna, generic_dna).translate())
    if position_range:
        start, stop = map(int, position_range.split('-'))
        position_list = range(start, stop + 1)
        if len(wild_type_aa) != len(position_list):
            raise Exception('length of wild type amino acid sequence and position range must be the same')
    else:
        position_list = range(1, len(wild_type_aa) + 1)

    input_df = pd.DataFrame.from_csv(input_csv)
    input_df.dropna(how='any', inplace=True)
    input_df = input_df[input_df.index.notnull()]
    mutation_series = pd.Series(0, index=input_df.index, name='Mutation')
    input_df = pd.concat([input_df, mutation_series], axis='columns')

    for seq in input_df.index:
        original_seq = seq
        if leading_to_drop:
            seq = seq[leading_to_drop:]
        if trailing_to_drop:
            seq = seq[:-trailing_to_drop]

        if len(seq) == len(wild_type_dna):
            if seq == wild_type_dna:
                input_df.loc[original_seq, 'Mutation'] = 'wild type'
            else:
                seq_aa = str(Seq(seq, generic_dna).translate())
                if '*' in seq_aa:
                    input_df.loc[original_seq, 'Mutation'] = 'stop'
                elif seq_aa == wild_type_aa:
                    input_df.loc[original_seq, 'Mutation'] = 'silent'
                else:
                    # mutation_string = compare_sequence_to_reference(seq_aa, wild_type_aa, position_list)
                    input_df.loc[original_seq, 'Mutation'] = 'missense'
        else:
            if (abs(len(seq) - len(wild_type_dna))) % 3 == 0:
                seq_aa = str(Seq(seq, generic_dna).translate())
                if '*' in seq_aa:
                    input_df.loc[original_seq, 'Mutation'] = 'stop'
                if len(seq) > len(wild_type_dna):
                    input_df.loc[original_seq, 'Mutation'] = 'insertion'
                else:
                    input_df.loc[original_seq, 'Mutation'] = 'deletion'
            else:
                input_df.loc[original_seq, 'Mutation'] = 'shift'

    input_df.to_csv('{0}_mutations.csv'.format(input_csv.split('/')[-1].split('.csv')[0]))

    missense_df = input_df.loc[input_df['Mutation'] == 'missense']
    missense_series = pd.Series(0, index=missense_df.index, name='Number of Mutations')
    missense_df = pd.concat([missense_df, missense_series], axis='columns')
    for seq in missense_df.index:
        original_seq = seq
        if leading_to_drop:
            seq = seq[leading_to_drop:]
        if trailing_to_drop:
            seq = seq[:-trailing_to_drop]
        seq_aa = str(Seq(seq, generic_dna).translate())
        mutation_string = compare_sequence_to_reference(seq_aa, wild_type_aa, position_list)
        mutation_count = mutation_string.count('_') + 1
        missense_df.loc[original_seq, 'Mutation'] = mutation_string
        missense_df.loc[original_seq, 'Number of Mutations'] = mutation_count

    missense_df.to_csv('{0}_missense.csv'.format(input_csv.split('/')[-1].split('.csv')[0]))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""script to map amino acid mutations from an input wild type dna 
    sequence to sequences in an input csv file. output is a csv file with the sequences replaced by their 
    corresponding mutations from the wild type sequence in the format of 
    <wild type amino acid>_<position number>_<mutant amino acid>. output is first ordered by increasing number of 
    mutations, then in/dels w/o frame shift, then in/dels w/ frame shift""")

    parser.add_argument('-r', '--position_range',
                        help='start and end amino acid position which wild type dna sequence corresponds to. start '
                             'and end position should be separated by a hyphen')
    parser.add_argument('-l', '--leading_to_drop', type=int,
                        help='integer number of leading characters to drop from each input sequence')
    parser.add_argument('-t', '--trailing_to_drop', type=int,
                        help='integer number of trailing characters to drop from each input sequence')
    required = parser.add_argument_group('required')
    required.add_argument('-i', '--input_csv', required=True,
                          help='csv file with first column containing dna sequences to be compared to wild type '
                               'sequence')
    required.add_argument('-w', '--wild_type', required=True,
                          help='txt file with wild type dna sequence')
    args = parser.parse_args()
    map_mutations(args.input_csv, args.wild_type, args.position_range, args.leading_to_drop, args.trailing_to_drop)
