#!/usr/bin/env python
# Author:  csiu
# Created: 2015-04-24
import argparse
from ConfigParser import SafeConfigParser
import re
import os
import sys
import pandas as pd

import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects.lib import ggplot2
from rpy2.robjects import r

from utils import get_value_from_keycolonvalue_list, ensure_dir
import label as lb

usage = """Generate plots (pie + histogram) with R via rpy2

Depends on:
- pandas     (python module)
- rpy2       (python module)
- ggplot2    (R package)
"""
def _check_labelling(infile, labelfile):
    ## simple check
    with open(infile) as f:
        for l in f:
            info = l.strip().split('\t')[8].split(';')
            label = get_value_from_keycolonvalue_list('mirna_label', info)
            if label == '':
                isLabelled = False
            else:
                isLabelled = True
            break

    if isLabelled:
        return infile
    else:
        print '## No labelling is found, proceed with labelling...'
        outfile = '%s.label' % infile

        lb.main(infile, labelfile, outfile)
        return outfile

def _verify_valid_distance(infile):
    out_good = infile + '.validdistance'
    out_bad = infile + '.badpair'

    with open(out_bad, 'w') as outB:
        with open(out_good, 'w') as outG:
            with open(infile) as f:
                for l in f:
                    l = l.strip().split('\t')

                    info = l[8].split(';')
                    d    = get_value_from_keycolonvalue_list('distance', info)

                    if d == 'NA':
                        chrom  = l[0]
                        start  = l[3]
                        stop   = l[4]
                        strand = l[6]
                        mirna  = get_value_from_keycolonvalue_list('mirna_query', info)

                        badpair = 'chr%s:%s..%s,%s\t%s' % (chrom, start, stop, strand,
                                                           mirna)
                        outB.write(badpair + '\n')
                    else:
                        outG.write('\t'.join(l) + '\n')

    if os.stat(out_bad).st_size != 0:
        print "## There are some bad positions in your input file:"
        print "## chromosome or strand differences between TSS and miRNA pair"
        print out_bad
    else:
        os.remove(out_bad)

    return out_good

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
            mirbase_id = get_value_from_keycolonvalue_list('mirbase_id', info)
            mstart     = get_value_from_keycolonvalue_list('mirna_start', info)
            mstop      = get_value_from_keycolonvalue_list('mirna_start', info)
            label      = get_value_from_keycolonvalue_list('mirna_label', info)

            if label == '': label = 'NA'
            mirna = ','.join([chrom, mstart, mstop, strand])

            features = l[7].split(';')
            corr = get_value_from_keycolonvalue_list('corr', features)
            if get_value_from_keycolonvalue_list('mirna_prox', features) != 0:
                distance = get_value_from_keycolonvalue_list('distance', info)
                if distance == '': distance = 0

            dat[n] = [tss, mirna, mirbase_id, label, distance, abs(float(distance)), corr]

    dat = pd.DataFrame.from_dict(dat, orient='index')
    dat.columns = ['tss', 'mirna', 'mirbase_id', 'label', 'Distance', 'distance', 'correlation']
    return dat

def _item_findClosestPartner(dat, item, ignoreCorr):
    df = dat[dat['label'] != 'NA'] ## remove invalid pairs

    xindex = []
    for i in df[item].unique():
        subm = df[df[item] == i]

        ## min distance
        submd = subm[subm['distance'] == subm['distance'].min()]

        if ignoreCorr:
            l = submd
        else:
            ## absmax correlation
            l = submd['correlation'].astype(float).abs()
            l = l[l == l.max()]

        xindex.append(l.index[0])
    return dat[dat.index.isin(xindex)]

def _plt_percountr(dat, independentpdf=False, fname='xpercount.pdf'):
    def _filt_dat(dat, item, getlabel=True):
        df = pd.DataFrame(dat[item].value_counts())
        df.columns = ['count']
        if getlabel: df['label'] = [list(dat[dat[item] == i]['label'])[0] for i in df.index]
        n  = len(df)
        mx = max(df['count'])
        return df, n, mx

    dat = dat[dat['label'] != 'NA']

    ## NUMBER OF MIRNA PER TSS
    df, n, mx = _filt_dat(dat, 'tss', False)
    df = {'count': robjects.IntVector(df['count'])}
    df = robjects.DataFrame(df)

    pt = ggplot2.ggplot(df) + \
        ggplot2.geom_histogram(binwidth=1, origin=-.5, alpha=.5, position="identity") + \
        ggplot2.xlim(-.5, mx+1) + \
        ggplot2.aes_string(x='count') + \
        ggplot2.ggtitle('TSS [Total = %s]' % n) + \
        ggplot2.labs(x='Number of miRNA per TSS (max = %s)' % mx)

    pt_den = ggplot2.ggplot(df) + \
        ggplot2.aes_string(x='count', y='..density..') + \
        ggplot2.geom_density(binwidth=1, alpha=.5, origin=-.5) + \
        ggplot2.geom_histogram(binwidth=1, alpha=.33, position='identity', origin=-.5) + \
        ggplot2.ggtitle('TSS [Total = %s]' % n) + \
        ggplot2.labs(x='Number of miRNA per TSS (max = %s)' % mx)

    ## NUMBER OF TSS PER MIRNA
    df, n, mx = _filt_dat(dat, 'mirna')
    df = {'count': robjects.IntVector(df['count']),
          'label': robjects.StrVector(df['label']) }
    df = robjects.DataFrame(df)

    _pm = ggplot2.ggplot(df) + \
        ggplot2.geom_histogram(binwidth=1, origin=-.5, alpha=.5, position="identity") + \
        ggplot2.xlim(-.5, mx+1) + \
        ggplot2.ggtitle('miRNA [Total = %s]' % n)

    _pm_den = ggplot2.ggplot(df) + \
        ggplot2.geom_density(binwidth=1, alpha=.5, origin=-.5) + \
        ggplot2.geom_histogram(binwidth=1, alpha=.33, position='identity', origin=-.5) + \
        ggplot2.ggtitle('miRNA [Total = %s]' % n)

    ## not split by label
    pm     = _pm     + ggplot2.aes_string(x='count')
    pm_den = _pm_den + ggplot2.aes_string(x='count', y='..density..')

    ## split by label
    pms     = _pm     + ggplot2.aes_string(x='count', fill='label')
    pm_dens = _pm_den + ggplot2.aes_string(x='count', fill='label', y='..density..')

    ## add xlabelling (need to be added after aes_string)
    _xlab = ggplot2.labs(x='Number of TSS per miRNA (max = %s)' % mx)
    pm      += _xlab
    pm_den  += _xlab
    pms     += _xlab
    pm_dens += _xlab

    if independentpdf:
        grdevices = importr('grDevices')
        grdevices.pdf(fname)
        pt.plot()
        pt_den.plot()
        pm.plot()
        pm_den.plot()
        pms.plot()
        pm_dens.plot()
        grdevices.dev_off()
    else:
        pt.plot()
        pt_den.plot()
        pm.plot()
        pm_den.plot()
        pms.plot()
        pm_dens.plot()
    return

def _plt_distr(dat, col, title='', splitBy_pfill=True, pfill='label', independentpdf=False, fname='xdistr.pdf'):
    df = dat[dat[pfill] != 'NA'] ## remove invalid pairs
    n  = len(df)
    df = {col: robjects.FloatVector(list(df[col])),
          pfill: robjects.StrVector(list(df[pfill]))}
    df = robjects.DataFrame(df)

    pp = ggplot2.ggplot(df) + \
        ggplot2.ggtitle('%s [Total = %s]' % (title, n))

    ## Plot1: counts
    if splitBy_pfill:
        p1 = pp + ggplot2.aes_string(x=col, fill=pfill)
    else:
        p1 = pp + ggplot2.aes_string(x=col)

    ## Plot2: density
    if splitBy_pfill:
        p2 = pp + ggplot2.aes_string(x=col, fill=pfill, y='..density..')
    else:
        p2 = pp + ggplot2.aes_string(x=col, y='..density..')
    p2 = p2 + ggplot2.geom_density(alpha=.5, origin=-500)

    if col == 'distance':
        p1 = p1 + \
            ggplot2.geom_histogram(binwidth=1000, alpha=.5, position='identity', origin=-500) + \
            ggplot2.xlim(-1000, 51000)

        p2 = p2 + \
            ggplot2.geom_histogram(binwidth=1000, alpha=.33, position='identity', origin=-500) + \
            ggplot2.xlim(-1000, 51000)
    else:
        p1 = p1 + \
            ggplot2.geom_histogram(alpha=.5, position='identity')

        p2 = p2 + \
            ggplot2.geom_histogram(alpha=.33, position='identity')

        if col == 'correlation':
            p1 = p1 + ggplot2.xlim(-1.1, 1.1)
            p2 = p2 + ggplot2.xlim(-1.1, 1.1)

    if independentpdf:
        grdevices = importr('grDevices')
        grdevices.pdf(file=fname)
        p1.plot()
        p2.plot()
        grdevices.dev_off()
    else:
        p1.plot()
        p2.plot()
    return

def _plt_pier(dat, title='', rm_na=False, col='label',
              independentpdf=False, fname='xpire.pdf'):
    x = dat[col]
    if rm_na: x = x[x != 'NA']
    x = x.value_counts()
    _sum = sum(x)

    r.assign('n', len(x))

    slices = robjects.FloatVector(x)
    pct  = ['%0.1f%%' % float(float(i)*100/_sum) for i in x]
    lbls = ['%s %s' % i  for i in zip(x.index, pct)]
    lbls   = robjects.StrVector(lbls)

    ggcolors = r('''
      gg_color_hue <- function(n){
            hues = seq(15, 375, length=n+1)
            hcl(h=hues, l=65, c=100)[1:n]
            }

      rev(gg_color_hue(n))
    ''')

    if independentpdf:
        grdevices = importr('grDevices')
        grdevices.pdf(pdf)

    r.pie(slices,
          labels=lbls,
          init_angle=90,
          col=ggcolors,
          main=title)
    r.text(-1, -1, '[Total = %s]' % _sum)

    if independentpdf:
        grdevices.dev_off()
    return

def main(infile, outdir, config):
    outdir = os.path.abspath(outdir)
    ensure_dir(outdir, False)

    cparser = SafeConfigParser()
    cparser.read(config)
    labelfile = cparser.get('configs', 'labelfile')

    infile = _check_labelling(infile, labelfile)
    infile = _verify_valid_distance(infile)
    infile = _filterPredictionsByClass_reformat2gff(infile, outdir)

    bname = os.path.basename(infile)
    pdf_rplots = os.path.join(outdir, bname + '.plots.pdf')

    dat = _read_dat(infile)

    if all(dat['correlation'] == ''):
        ignoreCorr = True
    else:
        ignoreCorr = False

    dat_mirna = _item_findClosestPartner(dat, 'mirna', ignoreCorr)
    dat_tss   = _item_findClosestPartner(dat, 'tss',   ignoreCorr)

    print '## Generating plot file...'
    grdevices = importr('grDevices')
    grdevices.pdf(file=pdf_rplots)

    _plt_pier(dat, 'All predicted TSS-miRNA pairs', True)
    _plt_distr(dat, 'distance', 'All predicted tss-miRNA pairs', False)
    _plt_distr(dat, 'distance', 'All predicted tss-miRNA pairs')
    if not ignoreCorr:
        _plt_distr(dat, 'correlation', 'All predicted tss-miRNA pairs', False)
        _plt_distr(dat, 'correlation', 'All predicted tss-miRNA pairs')

    _plt_percountr(dat)

    _plt_pier(dat_tss, 'TSS (label from closest miRNA)')
    _plt_distr(dat_tss, 'distance', 'TSS to closest miRNA', False)
    _plt_distr(dat_tss, 'distance', 'TSS to closest miRNA')
    if not ignoreCorr:
        _plt_distr(dat_tss, 'correlation', 'TSS to closest miRNA', False)
        _plt_distr(dat_tss, 'correlation', 'TSS to closest miRNA')

    _plt_pier(dat_mirna, 'miRNA')
    _plt_distr(dat_mirna, 'distance', 'miRNA to closest TSS', False)
    _plt_distr(dat_mirna, 'distance', 'miRNA to closest TSS')
    if not ignoreCorr:
        _plt_distr(dat_mirna, 'correlation', 'miRNA to closest TSS', False)
        _plt_distr(dat_mirna, 'correlation', 'miRNA to closest TSS')

    grdevices.dev_off()
    print '## Plot file:'
    print pdf_rplots
    return pdf_rplots

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=usage,
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-i', '--infile', dest='infile',
                        required=True,
                        help='''path to input file;
e.g. output of "label.py"
"mirna_label:" should be contained in the info column (9)''')

    parser.add_argument('-o', '--outdir', dest='outdir',
                        default='../Testout-plot',
                        help='''specify path to output directory''')

    parser.add_argument('-c', '--config', dest='config',
                        default='config.ini',
                        help='path to config file; default="config.ini"')

    ##get at the arguments
    args = parser.parse_args()

    ## do something..
    main(args.infile, args.outdir, args.config)
