#!/usr/bin/env python
# Author:  csiu
# Created:
import argparse
import re
import os
import pandas as pd
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

def _read_dat(gff_infile):
    dat = {}
    n = 0
    with open(gff_infile) as f:
        for l in f:
            n += 1

            l = l.strip().split('\t')
            chrom = l[0]
            tstart = l[3]
            tstop  = l[4]
            strand = l[6]
            tss = ','.join([chrom, tstart, tstop, strand])

            info  = l[8].split(';')
            mirna = get_value_from_keycolonvalue_list('mirbase_id', info)
            label = get_value_from_keycolonvalue_list('mirna_label', info)
            if label == '': label = 'NA'

            features = l[7].split(';')
            corr = get_value_from_keycolonvalue_list('corr', features)
            if get_value_from_keycolonvalue_list('mirna_prox', features) != '0':
                distance = get_value_from_keycolonvalue_list('distance', info)

            dat[n] = [tss, mirna, label, distance, corr]

    dat = pd.DataFrame.from_dict(dat, orient='index')
    dat.columns = ['tss', 'mirna', 'label', 'distance', 'correlation']
    return dat


def main(infile, outdir):
    outdir = os.path.abspath(outdir)
    ensure_dir(outdir)

    infile = _filterPredictionsByClass_reformat2gff(infile, outdir)
    dat = _read_dat(infile)
    print dat


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
