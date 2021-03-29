#! /usr/bin/env python3
import argparse
import os
import re
import multiprocessing as mp
from itertools import zip_longest
#from Bio.Seq import Seq
#from Bio.Alphabet import generic_dna

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--align')
    return parser.parse_args()

def seq_analyzer(seq):
    nucs = ['A', 'C', 'G', 'T', 'N']
    gaps = ['-']
    nuc_count = 0
    nuc_positions = []
    start_stops = []
    split_seq = seq.split('\n', 1)
    if len(split_seq) == 2:
        name = split_seq[0]
        print(name)
        # print(split_seq[1])
        previous_type = 'start'
        dna_seq = split_seq[1].replace('\n', '')
        # print(dna_seq)
        for num, nuc in enumerate(dna_seq):
            nuc_type = ''
            if nuc.upper() in nucs:
                nuc_type = 'nuc'
                nuc_count+=1
                nuc_positions.append(num)
                # print(nuc)
            elif nuc.upper() in gaps:
                nuc_type = 'gap'
                # print(nuc)

            if nuc_type != previous_type:
                start_stops.append(num)
                # print(nuc_type)
                # print(previous_type)

            previous_type = nuc_type

            
        start_stops.append(len(dna_seq))
                

    
    
    print(nuc_count)
    print(start_stops)
            

def main():
    args = parse_args()

    input_file = open(args.align, 'r')
    read_input = input_file.read()

    split_input = read_input.split('>')

    for num, split_data in enumerate(split_input):
        if len(split_data) >= 1:
            print('\n\n\n')
            print("CURRENT SEQ NUMBER: ", num)
            seq_analyzer(split_data)

    
    




if __name__ == '__main__':
    main()