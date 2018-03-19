#!/usr/bin/env python
import sys
import math


def write_fasta(input_file):
    allowed_lengths = [16, 17, 18, 19]
    starting_chars = ['A', 'C', 'G', 'T']
    input_name = input_file.split('.')[0]
    with open(input_file, 'r') as f:
        for l, line in enumerate(f):
            if line[0] in starting_chars:
                sequence = line.split()[0]
                fold_ind = float(line.split()[4])
                if not math.isnan(fold_ind) and len(sequence) in allowed_lengths:
                    # sys.stdout.write('>{0}_{1}_{2}\n'.format(input_name, l, fold_ind))
                    sys.stdout.write('>{0}_{1}\n'.format(input_name, l))
                    sys.stdout.write('{0}\n'.format(sequence))


if __name__ == '__main__':
    input_file = sys.argv[1]
    write_fasta(input_file)
