#!/usr/bin/env python
# Author:  csiu
# Created:
import argparse
import re
import os
import sys
import pandas as pd
import numpy  as np

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects.lib import ggplot2

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
            if get_value_from_keycolonvalue_list('mirna_prox', features) != 0:
                distance = get_value_from_keycolonvalue_list('distance', info)

            dat[n] = [tss, mirna, label, distance, corr]

    dat = pd.DataFrame.from_dict(dat, orient='index')
    dat.columns = ['tss', 'mirna', 'label', 'distance', 'correlation']
    return dat

def _plt_pie(dat, pdf, title='', rm_na=False, col="label"):
    x = dat[col]
    if rm_na: x = x[x != 'NA']
    x = x.value_counts()

    plt.figure()
    plt.pie(x,
            labels = x.index,
            autopct='%1.1f%%',
            startangle=90)
    plt.text(-1,-1, '[Total = %s]' % sum(x))
    plt.axis('equal')
    plt.title(title)
    pdf.savefig()
    plt.close()
    return

def _plt_countdistr(dat, fname, col, title='', pfill='label'):
    df = dat[dat[pfill] != 'NA'] ## remove invalid pairs
    df = {col: robjects.FloatVector(list(df[col])),
          pfill: robjects.StrVector(list(df[pfill]))}
    df = robjects.DataFrame(df)

    grdevices = importr('grDevices')
    grdevices.pdf(file=fname)

    pp = ggplot2.ggplot(df) + \
        ggplot2.aes_string(x=col, fill=pfill) + \
        ggplot2.ggtitle(title)

    if col == 'distance':
        pp = pp + \
            ggplot2.geom_histogram(binwidth=1000, alpha=.5, position='identity', origin=-500) + \
            ggplot2.xlim(-1000, 51000)
    else:
        pp = pp + \
            ggplot2.geom_histogram(alpha=.5, position='identity')

    pp.plot()
    grdevices.dev_off()
    return

def _plt_dendistr(dat, fname, col, title='', pfill='label'):
    df = dat[dat[pfill] != 'NA'] ## remove invalid pairs
    df = {col: robjects.FloatVector(list(df[col])),
          pfill: robjects.StrVector(list(df[pfill]))}
    df = robjects.DataFrame(df)

    grdevices = importr('grDevices')
    grdevices.pdf(file=fname)

    pp = ggplot2.ggplot(df) + \
        ggplot2.aes_string(x=col, fill=pfill, y='..density..') + \
        ggplot2.geom_density(alpha=.5, origin=-500) + \
        ggplot2.ggtitle(title)

    if col == 'distance':
        pp = pp + \
            ggplot2.geom_histogram(binwidth=1000, alpha=.33, position='identity', origin=-500) + \
            ggplot2.xlim(-1000, 51000)
    else:
        pp = pp + \
            ggplot2.geom_histogram(alpha=.33, position='identity')

    pp.plot()
    grdevices.dev_off()
    return

def main(infile, outdir):
    outdir = os.path.abspath(outdir)
    ensure_dir(outdir)

    infile = _filterPredictionsByClass_reformat2gff(infile, outdir)
    dat = _read_dat(infile)

    ## FIXME
    pdf_outfile = 'test.df'
    pdf_outfile_distr_dist = 'xdistance.pdf'
    pdf_outfile_distr_corr = 'xcorrelation.pdf'
    pdf_outfile_distr_ddist = 'xdistanceden.pdf'
    pdf_outfile_distr_dcorr = 'xcorrelationden.pdf'

    with PdfPages(pdf_outfile) as pdf:
        _plt_pie(dat, pdf, 'All TSS-[miRNA,NA] pairs')
        _plt_pie(dat, pdf, 'All valid TSS-miRNA pairs', True)
    _plt_countdistr(dat, pdf_outfile_distr_dist, 'distance',    'All valid tss-miRNA pairs')
    _plt_dendistr(dat, pdf_outfile_distr_ddist,  'distance',    'All valid tss-miRNA pairs')
    _plt_countdistr(dat, pdf_outfile_distr_corr, 'correlation', 'All valid tss-miRNA pairs')
    _plt_dendistr(dat, pdf_outfile_distr_dcorr,  'correlation', 'All valid tss-miRNA pairs')

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
