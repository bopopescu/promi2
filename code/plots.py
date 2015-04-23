#!/usr/bin/env python
# Author:  csiu
# Created:
import argparse
import os
from utils import get_value_from_keycolonvalue_list, ensure_dir

usage = """
"""
def _filterPredictionsByClass_reformat2gff(infile, outdir, keep='prom'):
    outfile = os.path.join(outdir, os.path.basename(infile)+'.filtered')
    with open(outfile, 'w') as out:
        with open(infile) as f:
            for l in f:
                l = l.strip().split('\t')
                classlabel = l[13]

                if classlabel == keep:
                    newinfo = ';'.join([l[8],
                                        'prior_prom:' + l[9],
                                        'prior_back:' + l[10],
                                        'prob_prom:' + l[11],
                                        'prob_back:' + l[12],
                                        'class:' + classlabel])
                    newline = '\t'.join(l[0:8] + [newinfo])
                    out.write(newline + '\n')
    return outfile

def _read_data(gff_infile):
    with open(gff_infile) as f:
        for l in f:
            pass
    pass

def main(infile, outdir):
    outdir = os.path.abspath(outdir)
    ensure_dir(outdir)

    infile = _filterPredictionsByClass_reformat2gff(infile, outdir)

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
