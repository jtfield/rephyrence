#! /usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import sys
import os
import statistics
import numpy
import dendropy
import nested_clade

import pandas as pd
import numpy as np
import scipy
# from sklearn.preprocessing import StandardScaler
import scipy.stats
import matplotlib.pyplot as plt


VERSION = "0.0.5"

class Rephyrence(object):
    """A Rephyrence object that holds the original input tree, distribution of branch lengths
    and best references for various clades"""
    
    def __init__(self, tree):
        """initialized object."""
        self.tree = ''
        self.tns = dendropy.TaxonNamespace()
        self.dendropy_tree = None
        self.main = 1
        self.branch_lengths = []
        if tree == 'test':
            self.tree = '(taxon_10:0.000490504536564,\
            (taxon_15:0.00147793977114,\
            ((taxon_22:0.00775296610854,\
            (taxon_11:0.00233219469227,\
            taxon_19:0.00292833772938):0.00661072260997)\
            :0.00437388771198,\
            ((taxon_17:0.00238143585578,\
            ((taxon_20:0.00040300210567,\
            ((taxon_27:0.000706168118649,taxon_23:0.000100735490224):1.00000050003e-06,\
            (taxon_25:0.000503957096912,\
            (taxon_26:1.00000050003e-06,\
            taxon_21:0.000100729540367):0.000403230214661)\
            :1.00000050003e-06):1.00000050003e-06):0.000200876170437,\
            (taxon_28:0.0011141033398,taxon_24:9.95483145232e-05)\
            :0.000507160142818):0.00184911046943):0.00339928596404,\
            (taxon_16.ref:0.000400637478349,\
            (taxon_13:0.000100794665489,(taxon_14:0.000706692317158,\
            (taxon_18:1.00000050003e-06,taxon_1:1.00000050003e-06)\
            :0.000302521801106):1.00000050003e-06):0.00061029754872)\
            :0.00511761213428):0.000906053315062):0.0040106698136)\
            0.00276938494406,taxon_12:0.00214526732704):0.0;'
            sys.stdout.write("Running Rephyrence using the incuded test phylogeny\n")
        else:
            if os.path.isfile(tree):
                self.tree = tree
                sys.stdout.write("Running Rephyrence using tree file {}\n".format(self.tree))
            else:
                sys.stderr.write("Tree file '{}' not found. Exiting.\n".format(tree))
                self._exit_handler()

        # self.print_tree()
        self.read_tree()
        self.print_tree_length()
        self.get_branch_lengths(self.dendropy_tree.seed_node)
        print(self.branch_lengths)

    def _exit_handler(self):
        '''makes debugging interactively easier, by not exiting on errors'''
        if self.main:
            sys.exit()
        else:
            pass

    def print_tree(self):
        '''prints the tree string'''
        print(self.tree)

    def read_tree(self):
        '''reads the tree string into dendropy and assigns that dendropy tree object to self.dendropy_tree'''
        self.dendropy_tree = dendropy.Tree.get(data=self.tree, schema='newick', taxon_namespace=self.tns)

    def print_tree_length(self):
        print("tree length is {}\n".format(self.dendropy_tree.length()))


    def get_branch_lengths(self, node):
        if node.edge.length != None:
            self.branch_lengths.append(node.edge.length)
        
        for child in node.child_nodes():
            self.get_branch_lengths(child)




    

        





# def parse_args():
#     parser = argparse.ArgumentParser()
#     parser.add_argument('--tree', help='phylogenetic tree file in newick format')
#     parser.add_argument('--test', default=False, help='included test tree for debugging and development')
#     parser.add_argument('--dist', default='short', help='branch length distance leading to references (DEFAULT:short)(Options: short, long')
#     return parser.parse_args()

# def identify_nests(tree, taxa):
#     mle_len = tree.length()
#     print("tree length is", mle_len)

    # for split in tree.encode_bipartitions():
        # print(split.leafset_bitmask)
        # print(dir(split))
        # for split2 in tree.encode_bipartitions():
        #     if split.is_nested_within(split2):
        #         print(split.leafset_as_newick_string(taxa))
        #         print(split2.leafset_as_newick_string(taxa))
        #         print("is nested")

# def process_node(tree_start_node):
#     print('XXXXXXXXXXXXXXXXXX')
#     print('prime_node: ', tree_start_node)
#     print(tree_start_node.edge.length)
#     for child in tree_start_node.child_nodes():
#         process_node(child)

# def main():
#     args = parse_args()

#     #Establish namespace.
#     taxa = dendropy.TaxonNamespace()

#     mle = ''


#     mle_1 = '(taxon_10:0.000490504536564,\
#         (taxon_15:0.00147793977114,\
#         ((taxon_22:0.00775296610854,\
#         (taxon_11:0.00233219469227,\
#         taxon_19:0.00292833772938):0.00661072260997)\
#         :0.00437388771198,\
#         ((taxon_17:0.00238143585578,\
#         ((taxon_20:0.00040300210567,\
#         ((taxon_27:0.000706168118649,taxon_23:0.000100735490224):1.00000050003e-06,\
#         (taxon_25:0.000503957096912,\
#         (taxon_26:1.00000050003e-06,\
#         taxon_21:0.000100729540367):0.000403230214661)\
#         :1.00000050003e-06):1.00000050003e-06):0.000200876170437,\
#         (taxon_28:0.0011141033398,taxon_24:9.95483145232e-05)\
#         :0.000507160142818):0.00184911046943):0.00339928596404,\
#         (taxon_16.ref:0.000400637478349,\
#         (taxon_13:0.000100794665489,(taxon_14:0.000706692317158,\
#         (taxon_18:1.00000050003e-06,taxon_1:1.00000050003e-06)\
#         :0.000302521801106):1.00000050003e-06):0.00061029754872)\
#         :0.00511761213428):0.000906053315062):0.0040106698136)\
#         0.00276938494406,taxon_12:0.00214526732704):0.0;'


#     #Handle whether working with test tree or actual input tree.
#     if args.test == "True":
#         mle = dendropy.Tree.get(data=mle_1, schema='newick', taxon_namespace=taxa)
#     else:
#         mle = dendropy.Tree.get(path=args.tree_file, schema='newick', taxon_namespace=taxa)

#     # nest_list = identify_nests(mle, taxa)

#     process_node(mle.seed_node)

    #TODO: implement branch length distribution fitting to assess best standard deviation cuttoff


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""Rephyrence - a Python program to identify the best reference for various clades in a phylogeny""",
        epilog="""Still in development, contact - jtoscanifield@ucmerced.edu - with bugs or issues"""
    )
    parser.add_argument('--tree', default='test', help='phylogenetic tree file in newick format')
    # parser.add_argument('--test', default=False, help='included test tree for debugging and development')
    # parser.add_argument('--dist', default='short', help='branch length distance leading to references (DEFAULT:short)(Options: short, long')
    args = parser.parse_args()
    rephy = Rephyrence(tree=args.tree)