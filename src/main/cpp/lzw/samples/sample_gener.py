import subprocess
import argparse
import random
import numpy as np
import os
import sys

# -----------------------------------------------------------------------------------------

def split_file(args):
    filename = args.source
    with open(filename) as fin:
        fout = open(args.destination_dir + "/0.samp", "w")
        for i, line in enumerate(fin):
            fout.write(line)
            if (i + 1) % args.n_of_lines_on_file == 0:
                fout.close()
                fout = open(args.destination_dir + "/%d.samp"% (i / args.n_of_lines_on_file + 1),"w")

        fout.close()

# -----------------------------------------------------------------------------------------

def add_parse_arguments(parser):
    parser.add_argument("-nl",  "--n_of_lines_on_file", type = int, default = 50)
    parser.add_argument('-src',  "--source", type = str, default =
                    os.path.dirname(os.path.abspath(__file__)) + '/big_samples/big_src.samp')
    parser.add_argument('-dest',  '--destination_dir',          type = str, default =
                    os.path.dirname(os.path.abspath(__file__)) + '/src_samples/')

# -----------------------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    add_parse_arguments(parser)
    args = parser.parse_args()

    split_file(args)
