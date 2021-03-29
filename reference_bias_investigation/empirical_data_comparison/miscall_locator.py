#! /usr/bin/env python3
import argparse
import os
import re


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_file')
    #parser.add_argument('--prefix', default='')
    return parser.parse_args()


def main():
    args = parse_args()

    with open(args.input_file) as myfile:
            head = [next(myfile) for x in range(12)]
            #print(head)
            identical_nucs = head[2]
            miscalled_bases = head[4]
            gaps = head[6]
            total_nucs = head[8]
            identical_positions = head[10]

            positions = 0

            identical_positions = ' '.join(identical_positions)

            positions_list = []
            for item in identical_positions.split(','):
                positions_list.append(item.strip().strip('[').strip(']').replace(' ', ''))
                

            for num, pos in enumerate(positions_list):
                # print(pos)
                if positions != 0:
                    if int(pos) != (positions + 1):
                        print(pos)
                positions = int(pos)



if __name__ == '__main__':
    main()