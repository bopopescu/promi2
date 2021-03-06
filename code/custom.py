#!/usr/bin/env python
# Author:  csiu
# Created: 2015-04-22
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
import gff_unify_features
import promi2
import label

usage = """
Action:
- 0. Given gff infile containing a list of positions
     (and possibly miRNA partners)
- 1. Extracts cpg,cons,tata,mprox,corr features
- 2. Calculates probabilty of promoter
- 3. Labels position as intra/intergenic based on miRNA
- 4. Generate plots (when "-p" is given)

Notes:
- 'tcnorm.ini' needs to be filled in
- for enabling plot generation, add "-p"
- when miRNA partner is given, add "-m"

Examples:
python2.7 custom.py -i ../test/test-custom.gff -o ../Testout-custom -p [-m]
"""
def _verify_infile(infile):
    ## should not contain chromosome M
    with open(infile) as f:
        for l in f:
            chrom = l.split('\t')[0]
            chrom.upper()
            if 'M' in chrom:
                sys.exit('''## ERROR: please remove chrM from the input file\n%s''' %
                         l)
    return

def _reformat_infile_gff2tcnorm(infile, outfile):
    with open(outfile, 'w') as out:
        pos = []
        with open(infile) as f:
            oldtss = ''
            for l in f:
                l = l.strip().split('\t')
                chrom = l[0]
                start = l[3]
                stop  = l[4]
                strand = l[6]

                tss = 'chr%s:%s..%s,%s' % (chrom, start, stop, strand)

                ## crude way to remove duplicate entries
                if tss != oldtss:
                    pos.append(tss)
                    oldtss = tss

        ## finer way to remove duplicate entries
        for newline in list(set(pos)):
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

def _interpret_tss_mirna_pairings(gff_infile, gff_tss, gff_mirna, pair_pos):
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

    tss_lookup = {}
    with open(gff_tss) as f:
        for l in f:
            chrom, _, _, start, stop, _, strand, r, _ = l.split('\t')
            tss = 'chr%s:%s..%s,%s' % (chrom, start, stop, strand)
            tss_lookup[tss] = r

    with open(pair_pos, 'w') as out:
        with open(gff_infile) as f:
            for l in f:
                chrom, _, _, start, stop, _, strand, _, mirna = l.strip().split('\t')
                tss_annot = 'title=chr%s:%s..%s,%s' % (chrom, start, stop, strand)

                n = tss_lookup['chr%s:%s..%s,%s' % (chrom, start, stop, strand)]

                mirna = mirna.lower()
                for k in [m for m in mirna_lookup_keys if re.match('%s(-.*$|$)' % mirna, m)]:
                    r, mirna_annot = mirna_lookup[k]
                    pairid = 'pair_id=' + '.'.join([chrom, start, stop, strand,
                                                    re.sub('-[35]p$', '', k)])

                    newline = '\t'.join([tss_annot, str(n),
                                         ','.join([mirna_annot, 'mirna_query='+mirna, pairid]),
                                         r])

                    out.write(newline + '\n')
    return

def _index_feat(gff_ufeat, has_mirna):
    pairid_index = {}
    with open(gff_ufeat) as f:
        c = 0
        for l in f:
            c += 1

            chrom, _, _, start, stop, _, strand, _, info = l.strip().rsplit('\t')
            info    = re.split('[;@]', info)

            pid = '.'.join([chrom, start, stop, strand])
            if has_mirna:
                mirna = get_value_from_keycolonvalue_list('mirna_query', info)
                val = '%s:%s' % (mirna, c)
            else:
                val = c

            try:
                pairid_index[pid].append(val)
            except KeyError:
                pairid_index[pid] = [val]

    return pairid_index

def _reformat_tss_to_1kb(posfile, gff1kb_infile):
    with open(gff1kb_infile, 'w') as out:
        with open(posfile) as f:
            for l in f:
                l = l.strip().split('\t')[0]
                _, chrom, start, _, stop, strand = re.split('[r:.,]', l)

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

def _interpret_mprox(gff_infile, mirbase_gff2, outfile):
    mbase_lookup = {}
    with open(mirbase_gff2) as f:
        for l in f:
            if l.startswith('#'): continue

            l = l.strip().split('\t')
            chrom = re.sub('^[Cc]hr', '', l[0])
            start = l[3]
            stop  = l[4]
            strand = l[6]

            info = re.split('["]', l[8])
            macc = info[1]
            mid  = info[3].lower()

            try:
                mbase_lookup[mid].append([chrom, start, stop, strand, macc])
            except KeyError:
                mbase_lookup[mid] = [[chrom, start, stop, strand, macc]]

    mbase_mirna = mbase_lookup.keys()

    with open(outfile, 'w') as out:
        with open(gff_infile) as f:
            for l in f:
                chrom, _, _, tstart, tstop, _, strand, _, mirna = l.strip().split('\t')
                mirna = mirna.lower()

                if mirna in mbase_mirna:
                    for m in mbase_lookup[mirna]:
                        mchrom, mstart, mstop, mstrand, macc = m

                        if chrom != mchrom or strand != mstrand:
                            d = 'NA'
                            mprox = 0
                        else:
                            d = mirna_proximity.calculate_distance(tstart, tstop, mstart, mstop, strand)
                            mprox = mirna_proximity.distance_score(d)

                        newinfo = ';'.join(['distance:%s' % d,
                                            'mirna_acc:%s' % macc,
                                            'mirbase_id:%s' % mirna,
                                            'mirna_start:%s' % mstart,
                                            'mirna_stop:%s' % mstop,
                                            'mirna_query:%s' % mirna])

                        newline = '\t'.join([chrom, '.', mirna,
                                             tstart, tstop, str(mprox), strand,
                                             '.', newinfo])
                        out.write(newline + '\n')
                else:
                    ## Here is for when listed mirna is not found in mirbase
                    ## FIXME: right now assumes listed mirna is found in mirbase
                    pass
    return outfile

def extractFeatures_given_gff(config, gff_infile, outdir, has_mirna, is_consider_corr):
    cparser = SafeConfigParser()
    cparser.read(config)

    tc_config = cparser.get('configs', 'tcconfig')
    m_mirna   = cparser.get('correlation', 'srnaseqmatrix')

    f_fasta      = cparser.get('genome','fasta')
    f_chromsizes = cparser.get('genome','chromsizes')
    d_phastcons  = cparser.get('cons','phastcons')
    TRAP         = cparser.get('tata','trap')
    f_psemmatrix = cparser.get('tata','psem')

    mirbase_gff2 = cparser.get('mirbase', 'gff2')
    corrmethod   = cparser.get('correlation', 'corrmethod')

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
        for line in f:
            l = line.strip().split('\t')
            try:
                _, chrom, start, _, stop, strand = re.split('[r:.,]', l[3])
                pos = '.'.join([chrom, start, stop, strand])
                ncount_dict[pos] = l[6]
            except ValueError:
                print '#[tcBedSpltErr]: %s' % line,

    ## 1d. setup outfile
    f_rle = re.sub('max_tpm.bed$', 'tpm_rle.matrix', fo_bed)

    tcparser = SafeConfigParser()
    tcparser.read(tc_config)
    f_ids = tcparser.get('tc_normalization', 'ids')

    ## PART2: compute cpg, cons, tata ...
    outdir_seqfeatures = os.path.join(outdir, 'seqfeatures/')
    ensure_dir(outdir_seqfeatures)

    gff1kb_infile   = os.path.join(outdir_seqfeatures, 'infile_1kbseq.gff')
    gff_1kbfeatures = os.path.join(outdir_seqfeatures, 'features_1kbseq.gff')

    _reformat_tss_to_1kb(f1_pos, gff1kb_infile)

    features.main(gff1kb_infile, outdir_seqfeatures,
                  f_fasta, f_chromsizes, d_phastcons, TRAP, f_psemmatrix,
                  gff_1kbfeatures)

    ## PART3: compute mprox ...
    outdir_tmp = os.path.join(outdir, 'intermediates')
    ensure_dir(outdir_tmp, False)

    gff_mproxfeatures = os.path.join(outdir_tmp, 'features_mprox.gff')
    gff_ufeat1 = os.path.join(outdir_tmp, 'features.1kb.mprox.gff')

    if has_mirna:
        _interpret_mprox(gff_infile, mirbase_gff2, gff_mproxfeatures)
    else:
        mirna_proximity.main(gff1kb_infile, mirbase_gff2, gff_mproxfeatures)

    gff_unify_features.main(gff_1kbfeatures, gff_mproxfeatures, 'mirna_prox', '0', gff_ufeat1, True)

    ## PART4: compute corr
    if is_consider_corr:
        ## correlation setup:
        outdir_corr = os.path.join(outdir, 'corr')
        ensure_dir(outdir_corr, False)

        gff_mirna = os.path.join(outdir_corr, '4corr_mirna.gff')
        gff_tss   = os.path.join(outdir_corr, '4corr_tss.gff')
        pair_pos    = os.path.join(outdir_corr, '4corrPair_row_pos_tss-mirna.txt')
        pair_sample = os.path.join(outdir_corr, '4corrPair_col_sample_CAGE-sRNAseq.txt')
        fo_corr = os.path.join(outdir_corr, 'features_correlation-%s.gff' % corrmethod)
        gff_ufeat2 = os.path.join(outdir_tmp, 'features.1kb.mprox.corr.gff')

        ## position pair:
        correlation._find_miRNA_pos(m_mirna, mirbase_gff2, gff_mirna)
        correlation._get_tss_pos(f1_pos, gff_tss)
        if has_mirna:
            _interpret_tss_mirna_pairings(gff_infile, gff_tss, gff_mirna, pair_pos)
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

        gff_unify_features.main(gff_ufeat1, fo_corr, 'corr', '0', gff_ufeat2, True)

        gff_ufeat = gff_ufeat2
    else:
        gff_ufeat = gff_ufeat1

    findex = _index_feat(gff_ufeat, has_mirna)

    ## PART4: start consolidating features ...
    gff_allfeatures = os.path.join(outdir, 'features.gff')
    with open(gff_allfeatures, 'w') as out:
        with open(gff_infile) as f:
            for l in f:
                chrom, _, _, start, stop, _, strand, _, mirna = l.strip().split('\t')
                mirna = mirna.lower()

                ## setting ids...
                tssid  = '.'.join([chrom, start, stop, strand])

                ## getting info...
                try:
                    ncount = ncount_dict[tssid]
                except KeyError:
                    ncount = '0'

                if findex.has_key(tssid):
                    for n in findex[tssid]:
                        if has_mirna:
                            m, n = n.split(':')
                            if m != mirna: continue

                        newline = linecache.getline(gff_ufeat, int(n))
                        newline = newline.split('\t')

                        newline[2] = mirna
                        newline[5] = ncount

                        out.write('\t'.join(newline))
    return gff_allfeatures

def _filter_keepValidPairs(gff_features):
    f_out = gff_features
    f_all = gff_features + '.intermediate'
    os.rename(gff_features, f_all)

    with open(f_out, 'w') as out:
        with open(f_all) as f:
            for l in f:
                info = l.strip().split('\t')[8]
                if 'mirbase_id' in info:
                    out.write(l)
    return f_out

def main(f_config, gff_infile, outdir, has_mirna, make_plots):
    ensure_dir(outdir)

    cparser = SafeConfigParser()
    cparser.read(f_config)
    f_params       = cparser.get('promi2', 'params')
    listoffeatures = cparser.get('promi2', 'features').split(',')
    labelfile = cparser.get('configs', 'labelfile')

    if 'corr' in listoffeatures:
        is_consider_corr = True
    else:
        is_consider_corr = False

    ## Make sure no chrM in infile
    _verify_infile(gff_infile)

    ## Extract features
    gff_allfeatures = extractFeatures_given_gff(f_config, gff_infile, outdir, has_mirna, is_consider_corr)

    ## Don't consider TSS which does not have a partner miRNA
    gff_allfeatures = _filter_keepValidPairs(gff_allfeatures)

    ## Run Promirna
    fo_predictions = os.path.join(outdir,
                                  'Predictions.%s.txt' % os.path.basename(gff_infile))
    promi2.promi2(f_params, listoffeatures, gff_allfeatures, fo_predictions)

    ## Label predictions
    fo_labelledpredictions = fo_predictions + '.label'
    label.main(fo_predictions, labelfile, fo_labelledpredictions)

    ## Generate plots
    if make_plots:
        import plots
        outdir_plt = os.path.join(outdir, 'plots')
        plots.main(fo_labelledpredictions, outdir_plt, f_config)

    return fo_labelledpredictions

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=usage,
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-i', '--infile', dest='infile',
                        required=True,
                        help='''Path to input gff input file.
Should NOT include chromosome M.
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
Note: Listed miRNA must be in miRBase
''')

    parser.add_argument('-p', dest='make_plots',
                        action='store_true',
                        help='''Flag to enable plotting
This requires extra packages to be pre-installed:
- Python: pandas, matplotlib, rpy2
- R: ggplot2

''')

    parser.add_argument('-c', '--config', dest='f_config',
                        default='config.ini',
                        help='Path to config file; default="config.ini"')

    parser.add_argument('-o', '--outdir', dest='outdir',
                        help='Specify output directory')

    ##get at the arguments
    args = parser.parse_args()

    ## do something..
    main(args.f_config, args.infile, args.outdir, args.has_mirna, args.make_plots)
