#!/usr/bin/env python
# Author:  csiu
# Created:
import argparse
from utils import get_value_from_keycolonvalue_list, ensure_dir

usage = """
"""
def _filter_predictions(infile, outdir, keep='prom'):
    pass

def _read_data(somefile):
    pass

def main(infile, outdir):
    print outdir
    pass

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=usage,
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-i', '--infile', dest='infile',
                        required=True,
                        help='''path to input file''')

    parser.add_argument('-o', '--outdir', dest='outdir',
                        default='.',
                        help='''specify path to output directory''')

    ##get at the arguments
    args = parser.parse_args()

    ## do something..
    main(args.infile, args.outdir)
