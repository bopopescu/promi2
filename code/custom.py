#!/usr/bin/env python
# Author:  csiu
# Created:
import argparse
from ConfigParser import SafeConfigParser
import sys
import os
import re
import linecache

from utils import get_value_from_keycolonvalue_list, ensure_dir, random_string
import tc_normalization
import correlation
import mirna_proximity
import features
import promi2

usage = """
- you will need to set parameters in tc-normalization
"""
def _reformat_infile_gff2tcnorm(infile, outfile):
    with open(outfile, 'w') as out:
        with open(infile) as f:
            for l in f:
                l = l.strip().split('\t')
                chrom = l[0]
                start = l[3]
                stop  = l[4]
                strand = l[6]

                newline = 'chr%s:%s..%s,%s' % (chrom, start, stop, strand)
                out.write(newline + '\n')
    return

def _index_tcnorm(f_ids):
    tc_index = {}
    with open(f_ids) as f:
        c = 5
        for tc_id in f:
            tc_id = tc_id.strip()
            c += 1

            tc_index[tc_id] = ['%s:%s' % (tc_id, c)]
    return tc_index

def _interpret_tss_mirna_pairings(gff_infile, gff_mirna, pair_pos):
    mirna_lookup = {}
    with open(gff_mirna) as f:
        for l in f:
            chrom, _, mirbaseid, start, stop, _, strand, r, info = l.strip().split('\t')
            srnaseqid = info.split(':')[1]

            mirna_annot = ','.join(['title=' + srnaseqid,
                                    'mirbase_id='  + mirbaseid,
                                    'mirna_start=' + start,
                                    'mirna_stop='  + stop
                                    ])
            mirna_lookup[srnaseqid] = [r, mirna_annot]
    mirna_lookup_keys =  mirna_lookup.keys()

    with open(pair_pos, 'w') as out:
        with open(gff_infile) as f:
            n = 0
            for l in f:
                n += 1

                chrom, _, _, start, stop, _, strand, _, mirna = l.strip().split('\t')
                tss_annot = 'title=chr%s:%s..%s,%s' % (chrom, start, stop, strand)

                mirna = mirna.lower()
                for k in [m for m in mirna_lookup_keys if re.match('%s(-.*$|$)' % mirna, m)]:
                    r, mirna_annot = mirna_lookup[k]
                    pairid = 'pair_id=' + '.'.join([chrom, start, stop, strand,
                                                    re.sub('-[35]p$', '', k)])

                    newline = '\t'.join([tss_annot, str(n),
                                         ','.join([mirna_annot, pairid]), r])

                    out.write(newline + '\n')
    return

def _index_corr_pairid(gff_corr):
    pairid_index = {}
    with open(gff_corr) as f:
        c = 0
        for l in f:
            c += 1

            chrom, _, _, start, stop, _, strand, _, info = l.strip().rsplit('\t')
            info    = info.split(';')

            pid = get_value_from_keycolonvalue_list('pair_id', info)
            if pid == '': pid = '.'.join([chrom, start, stop, strand])

            try:
                pairid_index[pid].append(c)
            except KeyError:
                pairid_index[pid] = [c]

    return pairid_index

def _reformat_tss_to_1kb(gff_infile, gff1kb_infile):
    with open(gff1kb_infile, 'w') as out:
        with open(gff_infile) as f:
            for l in f:
                chrom, _, _, start, stop, _, strand, _, _ = l.strip().split()

                newinfo = ';'.join(['start:' + start,
                                    'stop:'  + stop])

                start = int(start)
                stop  = int(stop)
                mid = (start + stop) / 2

                newline = '\t'.join([chrom, '.', '.',
                                     str(mid - 500),
                                     str(mid + 500),
                                     '.', strand, '.', newinfo])
                out.write(newline + '\n')
    return

def _index_1kbfeatures(gff_1kbfeatures):
    feature_index = {}
    with open(gff_1kbfeatures) as f:
        c = 0
        for l in f:
            c += 1
            chrom, _, _, start, stop, _, strand, features, _ = l.split('\t')

            pos = '.'.join([chrom, start, stop, strand])

            try:
                feature_index[pos].append(c)
            except KeyError:
                feature_index[pos] = [c]

    return feature_index

def extractFeatures_given_posPairs(config, gff_infile, outdir, has_mirna):
    cparser = SafeConfigParser()
    cparser.read(config)

    tc_config = cparser.get('configs', 'tcconfig')
    m_mirna   = cparser.get('correlation', 'srnaseqmatrix')

    ## PART1: tc normalization
    ## 1a. setup infile
    outdir_tc = os.path.join(outdir, 'tc-norm')
    f1_pos    = os.path.join(outdir, 'f1_pos.txt')
    ensure_dir(outdir_tc)

    ## 1b. reformat infile so that can be read by tc-quantify
    _reformat_infile_gff2tcnorm(gff_infile, f1_pos)

    ## 1c. run
    fo_bed = tc_normalization.main(tc_config, f1_pos, outdir_tc)

    ncount_dict = {}
    with open(fo_bed) as f:
        for l in f:
            l = l.strip().split('\t')
            _, chrom, start, _, stop, strand = re.split('[r:.,]', l[3])
            pos = '.'.join([chrom, start, stop, strand])
            ncount_dict[pos] = l[6]

    ## 1d. setup outfile
    f_rle = re.sub('max_tpm.bed$', 'tpm_rle.matrix', fo_bed)

    tcparser = SafeConfigParser()
    tcparser.read(tc_config)
    f_ids = tcparser.get('tc_normalization', 'ids')


    ## PART2: correlation setup
    outdir_corr = os.path.join(outdir, 'corr/')
    ensure_dir(outdir_corr)

    mirbase_gff2 = cparser.get('mirbase', 'gff2')
    corrmethod   = cparser.get('correlation', 'corrmethod')
    gff_mirna = os.path.join(outdir_corr, '4corr_mirna.gff')
    gff_tss   = os.path.join(outdir_corr, '4corr_tss.gff')
    pair_pos    = os.path.join(outdir_corr, '4corrPair_row_pos_tss-mirna.txt')
    pair_sample = os.path.join(outdir_corr, '4corrPair_col_sample_CAGE-sRNAseq.txt')
    fo_corr = os.path.join(outdir_corr, 'features_correlation-%s.gff' % corrmethod)

    ## position pair:
    correlation._find_miRNA_pos(m_mirna, mirbase_gff2, gff_mirna)
    correlation._get_tss_pos(f1_pos, gff_tss)
    if has_mirna:
        _interpret_tss_mirna_pairings(gff_infile, gff_mirna, pair_pos)
    else:
        correlation._get_tss_mirna_pairings(gff_tss, gff_mirna, pair_pos)

    ## sample pair:
    srnaseq_index = correlation._index_srnaseq(m_mirna)
    cage_index    = _index_tcnorm(f_ids)
    correlation._get_sample_pairings(cage_index, srnaseq_index, pair_sample)

    ## compute correlation:
    correlation._compute_correlation(pair_pos, pair_sample,
                                     f_rle, m_mirna,
                                     fo_corr, corrmethod, '.')

    pairid_index = _index_corr_pairid(fo_corr)

    ## PART3: compute features
    ## compute cpg, cons, tata ...
    outdir_seqfeatures = os.path.join(outdir, 'seqfeatures/')
    ensure_dir(outdir_seqfeatures)

    f_fasta      = cparser.get('genome','fasta')
    f_chromsizes = cparser.get('genome','chromsizes')
    d_phastcons  = cparser.get('cons','phastcons')
    TRAP         = cparser.get('tata','trap')
    f_psemmatrix = cparser.get('tata','psem')

    gff1kb_infile   = os.path.join(outdir_seqfeatures, 'infile_1kbseq.gff')
    gff_1kbfeatures = os.path.join(outdir_seqfeatures, 'features_1kbseq.gff')

    _reformat_tss_to_1kb(gff_infile, gff1kb_infile)

    features.main(gff1kb_infile, outdir_seqfeatures,
                  f_fasta, f_chromsizes, d_phastcons, TRAP, f_psemmatrix,
                  gff_1kbfeatures)

    features_index = _index_1kbfeatures(gff_1kbfeatures)

    ## start consolidating features ...
    gff_allfeatures = os.path.join(outdir, 'features.gff')
    with open(gff_allfeatures, 'w') as out:
        with open(gff_infile) as f:
            for l in f:
                chrom, _, _, start, stop, _, strand, _, mirna = l.strip().split('\t')
                mirna = mirna.lower()

                ## setting ids...
                tssid  = '.'.join([chrom, start, stop, strand])

                if has_mirna:
                    pairid = '.'.join([tssid, mirna])
                else:
                    pairid = tssid

                ## getting info...
                ncount = ncount_dict[tssid]

                mirna_partner = []
                if pairid_index.has_key(pairid):
                    for n in pairid_index[pairid]:

                        ## feature: correlation (corr)
                        cline =  linecache.getline(fo_corr, n).strip().split('\t')
                        corr  = cline[5]
                        cinfo = cline[8].split(';')

                        ## feature: mirna (mprox)
                        mstart = int(get_value_from_keycolonvalue_list('mirna_start', cinfo))
                        mstop = int(get_value_from_keycolonvalue_list('mirna_stop', cinfo))

                        d = mirna_proximity.calculate_distance(start, stop, mstart, mstop, strand)
                        mprox = str(mirna_proximity.distance_score(d))

                        mirna_info = ';'.join([';'.join(cinfo[:-1]), 'distance:'+str(d)])
                        mirna_partner.append([corr, mprox, mirna_info])
                else:
                    mirna_partner.append(['0', '0', ''])

                ## features: cpg, cons, tata (cct)
                for m in mirna_partner:
                    corr, mprox, mirna_info = m

                    for t in features_index[tssid]:
                        fline = linecache.getline(gff_1kbfeatures, t).strip().split('\t')
                        cct = fline[7]
                        info_region = fline[8]

                        ## write out!
                        newfeatures = ';'.join([cct,
                                                'mirna_prox:' + mprox,
                                                'corr:'       + corr])

                        newinfo = ';'.join([info_region, mirna_info])

                        newline = '\t'.join([chrom, '.', '.',
                                             start, stop, ncount,
                                             strand, newfeatures, newinfo])
                        out.write(newline + '\n')

    return gff_allfeatures

def main(f_config, gff_infile, outdir, has_mirna):
    outdir = '../Testout-custom2'
    ensure_dir(outdir)

    cparser = SafeConfigParser()
    cparser.read(f_config)
    f_params       = cparser.get('promi2', 'params')
    listoffeatures = cparser.get('promi2', 'features').split(',')

    ## Extract features
    gff_allfeatures = extractFeatures_given_posPairs(f_config, gff_infile, outdir, has_mirna)

    ## Run Promirna
    fo_predictions = os.path.join(outdir,
                                  'Predictions.%s.txt' % os.path.basename(gff_infile))
    promi2.promi2(f_params, listoffeatures, gff_allfeatures, fo_predictions)

    ## Plot
    print 'good day'

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=usage,
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-i', '--infile', dest='infile',
                        required=True,
                        help='''Path to input gff input file.
Tab-separated columns should be like:
  1. chrom
  2. "." (source; not used)
  3. "." (feature; not used)
  4. putative_tss start
  5. putative_tss stop
  6. "." (score; not used)
  7. strand
  8. "." (frame; not used)
  9. "." (or miRNA if you used '-m' flag)
''')

    parser.add_argument('-m', dest='has_mirna',
                        action='store_true',
                        help='''Flag to designate that INFILE has miRNA in column 9
e.g. tss-mirna pairing is already provided &
program does not need to look for pairs

''')

    parser.add_argument('-c', '--config', dest='f_config',
                        default='config.ini',
                        help='Path to config file; default="config.ini"')

    parser.add_argument('-o', '--outdir', dest='outdir',
                        help='Specify output directory')

    ##get at the arguments
    args = parser.parse_args()

    ## do something..
    main(args.f_config, args.infile, args.outdir, args.has_mirna)
