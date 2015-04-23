#!/usr/bin/env python
# Author:  csiu
# Created:
import argparse
from utils import get_value_from_keycolonvalue_list, ensure_dir

usage = """
"""
def read_data(somefile):
    pass

def main(infile, outfile):
    print outfile
    pass

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=usage,
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-i', '--infile', dest='infile',
                        required=True,
                        help='''path to input file''')

    parser.add_argument('-o', '--outfile', dest='outfile',
                        help='''specify path to output file''')

    ##get at the arguments
    args = parser.parse_args()

    if args.outfile == None:
        outfile = '%s.plt' % args.infile
    else:
        outfile = args.outfile

    ## do something..
    main(args.infile, outfile)
