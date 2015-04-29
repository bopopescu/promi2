#!/usr/bin/env python
# Author:  csiu
# Created: 2015-03-05
import argparse
from ConfigParser import SafeConfigParser
import sys
import os
import re
import glob
import random
import shutil

from utils import random_string, get_value_from_keycolonvalue_list, ensure_dir, line_count
import features
import linecache
import tc_normalization
from correlation import _find_miRNA_pos, _compute_correlation
import mirna_proximity
import gff_unify_features

usage = """
Positive set
------------
- putative tss in at least min_lib
- putative tss is within 50kb upstream of pre-miRNA

Negative set
------------
- no repeats
- no sequences near exon1
- tc normalization (see tc_normalization.py)/overlap with peak
- cpg <= 0.5 and cons <= 0.2 and tata <= 0.1
- mirna_prox = 0

Memo
----
- specifying positive set @arg
- specifying negative set @tcnorm.ini
"""
def _position_counter(somefile, counter):
    with open(somefile) as f:
        for l in f:
            l = l.split('\t')
            chrom  = l[0]
            start  = l[3]
            stop   = l[4]
            strand = l[6]

            position = '%s,%s,%s,%s' % (chrom, start, stop, strand)

            try:
                counter[position] += 1
            except KeyError:
                counter[position] = 1
    return counter

def _readcount_finder(somefile, somelist, get_id=False):
    with open(somefile) as f:
        for l in f:
            l = l.split('\t')
            chrom  = l[0]
            start  = l[3]
            stop   = l[4]
            strand = l[6]

            position = '%s,%s,%s,%s' % (chrom, start, stop, strand)

            if position in somelist:
                x = l[5]

                if get_id:
                    info = l[8].strip().split(';')
                    sid = get_value_from_keycolonvalue_list('id', info)
                    somelist[position].append(sid+':'+x)
                else:
                    somelist[position].append(x)
    return somelist

def _50kb_upstream(pos_candidates, f_mirna, outfile):
    f_tmp = outfile + random_string(6)
    f_overlap = f_tmp + '.overlap'

    ## Pre-process
    with open(f_tmp, 'w') as out:
        for c in pos_candidates:
            chrom, start, stop, strand = c.split(',')
            chrom = 'chr' + chrom

            line = '\t'.join([chrom, '.', '.', start, stop, '.', strand, '.', '.'])
            out.write(line + '\n')

    ## Find those within 50kb upstream pre-miRNA
    cmd = 'bedtools window -a '+f_tmp+' -b '+f_mirna+' -l 0 -r 50000 -sw -sm -u > '+f_overlap
    os.system(cmd)

    ## Post-process
    new_candidates = []
    p_chr = re.compile('^chr')
    with open(f_overlap) as f:
        for l in f:
            l = l.split('\t')
            chrom  = p_chr.sub('', l[0])
            start  = l[3]
            stop   = l[4]
            strand = l[6]

            pos = ','.join([chrom, start, stop, strand])
            new_candidates.append(pos)

    ## Cleanup
    os.system('rm %s %s' % (f_tmp, f_overlap))

    return new_candidates

def _determine_region(start, stop):
    start = int(start)
    stop  = int(stop)

    mid = (start+stop)/2
    region_start = str(mid-500)
    region_stop  = str(mid+500)

    return region_start, region_stop

def create_positiveset(percent_lib, listoffiles, f_mirbasegff,
                       N, outfile, is_get_id):
    min_lib = len(listoffiles) * percent_lib

    pos_candidates = {}
    for i in range(0, len(listoffiles)):
        pos_candidates = _position_counter(listoffiles[i], pos_candidates)

    ## candidate positions should be in more than min_lib files
    pos_candidates = [k for k,v in pos_candidates.iteritems() if v>=min_lib]

    ## ensure that candidates are within 50kb of pre-miRNA
    pos_candidates = _50kb_upstream(pos_candidates, f_mirbasegff, outfile)

    ## from candidate list, select sample of 1000 positions
    pos_candidates = random.sample(pos_candidates, N)

    ## get list of readcounts from the 1000 positions
    pos_candidates = {k:[] for k in pos_candidates}
    for i in range(0, len(listoffiles)):
        pos_candidates = _readcount_finder(listoffiles[i], pos_candidates, is_get_id)

    with open(outfile, 'w') as out:
        for k,v in pos_candidates.iteritems():
            chrom, start, stop, strand = k.split(',')

            ## count value is selected from random lib containing position
            if is_get_id:
                sid, count = random.choice(v).split(':')
                newinfo = ';'.join(['start:'+start, 'stop:'+stop,
                                    'libs:'+str(len(v)),
                                    'id:'+sid])
            else:
                count = random.choice(v)
                newinfo = ';'.join(['start:'+start, 'stop:'+stop,
                                    'libs:'+str(len(v))])

            region_start, region_stop = _determine_region(start, stop)

            newline = '\t'.join([chrom, 'CAGE', 'putative_tss',
                                 region_start, region_stop,
                                 count, strand, '0',
                                 newinfo])
            out.write(newline + '\n')
    return outfile


def extract_tss_from_ensembl(f_ensembl, f_fiveprimegff):
    ## extract annotated gene start sites (used for filtering putative TSS)
    cmd = 'grep -w exon '+f_ensembl+' | grep \'exon_number \"1\"\''+' > '+f_fiveprimegff
    os.system(cmd)
    return


def _read_chromsizes(f_chromsizes):
    chromosomes = {}
    with open(f_chromsizes) as f:
        for l in f:
            c, s = l.strip().split('\t')
            if (('random' not in c) and ('chrM' not in c) and ('chrUn' not in c)):
                c = re.sub('chr', '', c)
                chromosomes[c] = int(s)
    return chromosomes

def _generate_random_positions(N, chromsizes, outfile):
    with open(outfile, 'w') as out:
        for i in range(0,N):
            chrom  = random.choice(chromsizes.keys())
            strand = random.choice(['+', '-'])

            L     = random.randrange(1000,5000,1)
            start = random.randrange(1, chromsizes[chrom]-L, 1)
            stop  = start + L

            region_start, region_stop = _determine_region(start, stop)
            newinfo = ';'.join(['start:'+str(start), 'stop:'+str(stop), 'id:random'])

            newline = '\t'.join(['chr'+chrom, 'RANDOM', 'BACK',
                                 region_start, region_stop,
                                 '.', strand, '0',
                                 newinfo])
            out.write(newline + '\n')
    return


def create_negativeset(f_chromsizes,   ##to generate random positions
                       f_repeats,      ##to filter out repeats
                       f_fiveprimegff, ##to filter out genes
                       f_traincfg,     ##to normalize tag counts
                       N,
                       outfile):
    ## this is your background set
    ## it is created in 4 steps
    f1_candidates    = outfile+'.bg_tmp.gff'
    f2_norepeats     = f1_candidates + 'norepeats.gff'
    f2_norepeats_tmp = f2_norepeats + random_string(6)
    f3_norepeats_nogene = f2_norepeats + 'nogenes.gff'
    f4_tcfile        = f3_norepeats_nogene + '.tcfile.txt'
    f4_tcfile_outdir = os.path.join(os.path.dirname(outfile),
                                    'tc-norm_negSet')

    chromsizes = _read_chromsizes(f_chromsizes)

    ## generate 10000 random candidate sequences
    _generate_random_positions(N*10, chromsizes, f1_candidates)

    ## filter repeats
    cmd = 'intersectBed -a '+f1_candidates+' -b '+f_repeats+' -u > '+f2_norepeats_tmp
    os.system(cmd)

    p_chr = re.compile('^chr')
    with open(f2_norepeats, 'w') as out:
        with open(f2_norepeats_tmp) as f:
            for line in f:
                line = p_chr.sub('', line)
                out.write(line)

    ## filter genes
    cmd = 'bedtools window -a '+f2_norepeats+' -b '+f_fiveprimegff+' -l 1200 -r 50 -sm -v > '+f3_norepeats_nogene
    os.system(cmd)

    ## tag count normalization
    with open(f4_tcfile, 'w') as out:
        linenum = 1
        with open(f3_norepeats_nogene) as f:
            for l in f:
                l = l.split('\t')
                chrom  = l[0]
                start  = l[3]
                stop   = l[4]
                strand = l[6]

                pos = 'chr%s:%s..%s,%s' % (chrom, start, stop, strand)
                out.write('%s\t%s\n' % (pos, linenum))
                linenum += 1

    f4_tpm = tc_normalization.main(f_traincfg, f4_tcfile, f4_tcfile_outdir)

    pos_dict = {}
    with open(f4_tcfile) as f:
        for l in f:
            pos, lnum = l.strip().split('\t')
            pos_dict[pos] = int(lnum)
    with open(outfile, 'w') as out:
        with open(f4_tpm) as f:
            for l in f:
                l = l.strip().split('\t')
                pos = l[3]
                tpm = float(l[6])

                if tpm != 0:
                    lnum = pos_dict[pos]
                    line = linecache.getline(f3_norepeats_nogene, lnum)
                    line = line.split('\t')
                    line[5] = str(tpm)
                    out.write('\t'.join(line))

    ## cleanup
    os.system('rm %s' % f1_candidates)
    os.system('rm %s %s' % (f2_norepeats, f2_norepeats_tmp))
    os.system('rm %s' % f3_norepeats_nogene)
    os.system('rm %s' % f4_tcfile)
    return

def feature_closest_corr(f_querygff,
                         f_mirbasegff, ##sRNAseq to gff format
                         m_mirna, m_tss, m_back, ## matrices with expression value
                         f_tcfilesinput,         ## will determine columnID for m_back
                         method,       ## correlation method
                         outfile,
                         verbose = False):
    ## files
    d = outfile + '_intermediates'
    ensure_dir(d, False)

    fo_mirnagff = os.path.join(d, 'srnaseq_pos.gff')
    fo_closest = os.path.join(d, 'closest_tss-mirna.txt')
    f_pos_pairing    = os.path.join(d, 'pairing_position.txt')
    f_sample_pairing = os.path.join(d, 'pairing_sample.txt')
    fo_corr = os.path.join(d, 'closest_corr.gff')

    ## 1. sRNAseq to gff format
    _find_miRNA_pos(m_mirna, f_mirbasegff, fo_mirnagff)

    ## 2a. find closest pair
    cmd = 'bedtools closest -a '+f_querygff+' -b '+fo_mirnagff+' -s -iu -D a -t first > '+fo_closest
    if verbose: print "STATUS: finding closest pair..."
    if verbose: print cmd
    os.system(cmd)

    ## 2b. get pairing info: position
    ## -> seq_id, seq_line, mirna_info, mirna_line, label
    if verbose: print 'STATUS: identifying pairing info: position...'
    with open(f_pos_pairing + '.posSet', 'w') as out_pos:
        with open(f_pos_pairing + '.negSet', 'w') as out_neg:

            cageseq_dict = {}
            with open(m_tss) as f:
                linenum = 0
                for l in f:
                    linenum += 1
                    if l.startswith('#') \
                           or l.startswith('00Annotation') \
                           or l.startswith('01STAT'):
                        continue

                    pos, _ = l.split('\t', 1)

                    chrom, start, _, stop, strand = re.split('[:.,]', pos)
                    start, stop = _determine_region(start, stop)

                    pos = '%s:%s..%s,%s' % (chrom, start, stop, strand)

                    cageseq_dict[pos] = linenum

            background_dict = {}
            with open(m_back) as f:
                linenum = 0
                for l in f:
                    linenum += 1
                    pos = l.split('\t')[3]
                    background_dict[pos] = linenum

            with open(fo_closest) as f:
                for l in f:
                    l    = l.strip()
                    _, d = l.rsplit('\t', 1)
                    d    = int(d)

                    if (d >= 0) and (d <= 50000):
                        l     = l.split('\t')
                        label = l[2]

                        pos    = 'chr%s:%s..%s,%s' % (l[0], l[3], l[4], l[6])
                        seq_id = 'title=%s' % pos

                        mirna_line = l[16]
                        mirna_info = ','.join(['title='      + l[17].split(':')[1],
                                               'mirbase_id=' + l[11],
                                               'mirna_start='+ l[12],
                                               'mirna_stop=' + l[13]])

                        if label == 'BACK':
                            info = l[8].split(';')
                            pos  = 'chr%s:%s..%s,%s' % (
                                l[0],
                                get_value_from_keycolonvalue_list('region_start', info),
                                get_value_from_keycolonvalue_list('region_stop', info),
                                l[6])
                            try:
                                seq_line = str(background_dict[pos])
                                newline = '\t'.join([seq_id, seq_line, mirna_info, mirna_line])
                                out_neg.write(newline + '\n')
                            except KeyError:
                                continue

                        else:
                            try:
                                seq_line = str(cageseq_dict[pos])
                                newline = '\t'.join([seq_id, seq_line, mirna_info, mirna_line])
                                out_pos.write(newline + '\n')
                            except KeyError:
                                continue

    ## 3. get pairing info: sample
    ## -> sampleID, cage_column_index, srnaseq_matrix_column_index, (cid,mid)
    if verbose: print 'STATUS: identifying pairing info: samples...'
    cage_id_pattern    = re.compile('^tpm.*(CNhs.*\..*)$')
    back_id_pattern    = re.compile('^.*(CNhs.*?\..*?)\..*$')
    srnaseq_id_pattern = re.compile('^.*(SRh.*?\..*?)\.')

    cage_index = {}
    with open(m_tss) as f:
        for l in f:
            if l.startswith('00Annotation'):
                l = l.strip().split('\t')
                c = 0
                for header in l:
                    if header.startswith('tpm'):
                        cage_sample_id = cage_id_pattern.match(header).group(1)
                        cage_id = cage_sample_id.split('.')[1]
                        try:
                            cage_index[cage_id].append('%s:%s' % (cage_sample_id, c))
                        except KeyError:
                            cage_index[cage_id] = ['%s:%s' % (cage_sample_id, c)]
                    c += 1
                break

    back_index = {}
    with open(f_tcfilesinput) as f:
        line = 6
        for l in f:
            cage_sample_id = back_id_pattern.match(l).group(1)
            cage_id = cage_sample_id.split('.')[1]
            try:
                back_index[cage_id].append('%s:%s' % (cage_sample_id, line))
            except KeyError:
                back_index[cage_id] = ['%s:%s' % (cage_sample_id, line)]
            line += 1

    srnaseq_index = {}
    with open(m_mirna) as f:
        for l in f:
            if l.startswith('ID'):
                l = l.strip().split('\t')
                c = 0
                for header in l:
                    if header.endswith('.bam'):
                        srnaseq_sample_id = srnaseq_id_pattern.match(header).group(1)
                        srnaseq_id = srnaseq_sample_id.split('.')[1]
                        try:
                            srnaseq_index[srnaseq_id].append('%s:%s' % (srnaseq_sample_id, c))
                        except KeyError:
                            srnaseq_index[srnaseq_id] = ['%s:%s' % (srnaseq_sample_id, c)]
                    c += 1
                break

    ## combine
    with open(f_sample_pairing + '.posSet', 'w') as out:
        sample_ids = set(cage_index.keys()).intersection(srnaseq_index.keys())
        for k in sample_ids:
            for c in cage_index[k]:
                for m in srnaseq_index[k]:
                    cid, cindex = c.split(':')
                    mid, mindex = m.split(':')
                    out.write('\t'.join([k, cindex, mindex,
                                         '%s,%s' % (cid, mid)]) +'\n')

    with open(f_sample_pairing + '.negSet', 'w') as out:
        sample_ids = set(back_index.keys()).intersection(srnaseq_index.keys())
        for k in sample_ids:
            for c in back_index[k]:
                for m in srnaseq_index[k]:
                    cid, cindex = c.split(':')
                    mid, mindex = m.split(':')
                    out.write('\t'.join([k, cindex, mindex,
                                         '%s,%s' % (cid, mid)]) +'\n')

    ## 4. compute correlation
    if verbose: print 'STATUS: computing correlation (method)...' % method
    _compute_correlation(f_pos_pairing +'.posSet', f_sample_pairing +'.posSet',
                         m_tss, m_mirna,
                         fo_corr +'.posSet', method, 'putative_tss')
    _compute_correlation(f_pos_pairing +'.negSet', f_sample_pairing +'.negSet',
                         m_back, m_mirna,
                         fo_corr +'.negSet', method, 'background')

    with open(fo_corr, 'w') as out:
        with open(fo_corr +'.negSet') as f:
            for l in f:
                out.write(l)
        with open(fo_corr +'.posSet') as f:
            for l in f:
                out.write(l)

    os.remove(fo_corr + '.negSet')
    os.remove(fo_corr + '.posSet')

    ## 5. unify
    if verbose: print 'STATUS: creating "%s"' % outfile
    gff_unify_features.main(f_querygff, fo_corr, 'corr', '0', outfile, True)

    return fo_corr

def main(files, outdir, N, percent_lib, is_get_id, f_config, verbose = False):
    if os.path.isdir(outdir): sys.exit('## ERROR: "%s" already exists' % outdir)

    cparser = SafeConfigParser()
    cparser.read(f_config)
    verbose = True

    f_mirbasegff = cparser.get('mirbase','gff2')
    f_chromsizes = cparser.get('genome','chromsizes')
    f_repeats    = cparser.get('genome','repeats')
    f_ensembl    = cparser.get('genome','ensemblgtf')
    f_fasta      = cparser.get('genome','fasta')
    d_phastcons  = cparser.get('cons','phastcons')
    TRAP         = cparser.get('tata','trap')
    f_psemmatrix = cparser.get('tata','psem')
    f_traincfg   = cparser.get('configs','tcconfig')
    m_mirna      = cparser.get('correlation','srnaseqmatrix')
    m_tss        = cparser.get('correlation','cageseqmatrix')
    corrmethod   = cparser.get('correlation','corrmethod')

    f_trainingset = os.path.join(outdir, 'TrainingSet.gff')
    outdir1 = f_trainingset + '_intermediates'

    ensure_dir(outdir, False)
    ensure_dir(outdir1, False)

    _files = glob.glob(files)

    ## creating auxillary file for negative set
    f_fiveprimegff = '../data/hsa.five_prime.gff'
    if not os.path.exists(f_fiveprimegff):
        if verbose: print 'STATUS: creating "%s" auxillary file...' % f_fiveprimegff
        extract_tss_from_ensembl(f_ensembl, f_fiveprimegff)


    ## create training set
    gff_ts_pos = os.path.join(outdir1, 'trainingset_pos.gff')
    gff_ts_neg = os.path.join(outdir1, 'trainingset_neg.gff')
    if verbose: print 'STATUS: creating positive candidate set...'
    create_positiveset(percent_lib, _files, f_mirbasegff, N, gff_ts_pos, is_get_id)
    if verbose: print 'STATUS: creating negative candidate set...'
    create_negativeset(f_chromsizes, f_repeats, f_fiveprimegff, f_traincfg,
                       N, gff_ts_neg)

    shutil.move(os.path.join(outdir1, 'tc-norm_negSet'),
                os.path.join(outdir, 'tc-norm_negSet'))

    ## feature extraction: cpg, cons, tata (features.py)
    if verbose: print 'STATUS: extracting features cpg/cons/tata...'
    gff_1kbfeatures_pos = os.path.join(outdir1, 'features1kb_ts_pos.gff')
    gff_1kbfeatures_neg = os.path.join(outdir1, 'features1kb_ts_neg.gff')

    features.main(gff_ts_pos, outdir1,
                  f_fasta, f_chromsizes, d_phastcons, TRAP, f_psemmatrix,
                  gff_1kbfeatures_pos)

    features.main(gff_ts_neg, outdir1,
                  f_fasta, f_chromsizes, d_phastcons, TRAP, f_psemmatrix,
                  gff_1kbfeatures_neg)

    ## feature extraction: mirna_proximity
    if verbose: print 'STATUS: extracting features mirna_proximity...'
    gff_mirnaprox_pos = os.path.join(outdir1, 'featureMprox_ts_pos.gff')
    gff_mirnaprox_neg = os.path.join(outdir1, 'featureMprox_ts_neg.gff')
    mirna_proximity.main(gff_ts_pos, f_mirbasegff, gff_mirnaprox_pos)
    mirna_proximity.main(gff_ts_neg, f_mirbasegff, gff_mirnaprox_neg)

    gff_features_pos = os.path.join(outdir1, 'Features_ts_pos.gff')
    gff_features_neg = os.path.join(outdir1, 'Features_ts_neg.gff')
    gff_unify_features.main(gff_1kbfeatures_pos, gff_mirnaprox_pos,
                            'mirna_prox', '0', gff_features_pos, True)
    gff_unify_features.main(gff_1kbfeatures_neg, gff_mirnaprox_neg,
                            'mirna_prox', '0', gff_features_neg, True)

    ## create final training set ...
    ## where background must pass criteria: cpg <= 0.5 and cons <= 0.2 and tata <= 0.1 and mirna_prox == 0
    if verbose: print 'STATUS: creating final training set...'
    good_background = gff_features_neg + '_cpglt0.5-conslt0.2-tatalt0.1-mproxeq0.gff'
    with open(good_background, 'w') as out:
        with open(gff_features_neg) as f:
            for line in f:
                info = line.strip().split('\t')[7].split(';')
                cpg  = float(get_value_from_keycolonvalue_list('cpg', info))
                cons = float(get_value_from_keycolonvalue_list('cons', info))
                tata = float(get_value_from_keycolonvalue_list('tata', info))
                mprx = float(get_value_from_keycolonvalue_list('mirna_prox', info))

                if cpg<=0.5 and cons<=0.2 and tata<=0.1 and mprx==0:
                    out.write(line)

    wc = line_count(good_background)
    selectedlines = random.sample(range(1,wc+1), N)

    with open(f_trainingset, 'w') as out:
        ## writing negative set
        for l in selectedlines:
            out.write(linecache.getline(good_background, l))

        ## writing positive set
        with open(gff_features_pos) as f:
            ## when mirna_prox extraction feature was used,
            ## extracted all pairs within 50kb upstream mirna
            ## -> single tss could have many mirna
            ## take pair with min distance
            ## -> essential first entry
            pos_list = []
            for line in f:
                l = line.split('\t')
                pos = ','.join([l[0], l[3], l[4], l[6]])
                if not (pos in pos_list):
                    pos_list.append(pos)
                    out.write(line)

    if not (os.path.isfile(m_mirna) and os.path.isfile(m_tss)):
        return f_trainingset

    ## create final training set with feature:correlation of closest tss->miRNA ...
    if verbose: print 'STATUS: creating final training set with correlation of closest tss->miRNA...'
    f_trainingset2 = os.path.join(outdir, 'TrainingSet-corr.gff')
    m_back         = glob.glob('%s/tc-norm_negSet/*tpm_rle.matrix' % outdir)[0]
    f_tcfilesinput = os.path.join(outdir, 'tc-norm_negSet', 'files.txt' )

    feature_closest_corr(f_trainingset, f_mirbasegff,
                         m_mirna, m_tss, m_back,
                         f_tcfilesinput,
                         corrmethod, f_trainingset2)

    return f_trainingset2

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=usage,
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-f', '--files', dest='files',
                        required=True,
                        help='''path to gff input files.
This is for the creating the positive set.
Am assuming these files contain all putative tss.
You will need quotes (\') when wildcard (*) is used.''')

    parser.add_argument('-o', '--outdir', dest='outdir',
                        default='trainingset_'+random_string(6),
                        help='specify path to outdirectory')

    parser.add_argument('-c', '--config', dest='f_config',
                        default='config.ini',
                        help='path to config file; default="config.ini"')

    parser.add_argument('-n', dest='n',
                        type=int, default=1000,
                        help='''number of positive (or negative) positions in training set.
Total number of training samples = 2n
Default = 1000''')

    ## options for positive set
    parser.add_argument('-l', '--min_lib', dest='l',
                        type=float, default=0.3,
                        help='''minimum percent of libs
for the putative tss to be in
to be considered as a candidate
for the positive set.
Default = 0.30''')

    ## other options
    parser.add_argument('-g', dest='get_id',
                        action='store_true',
                        help='''flag to get "id" from info column of input''')

    parser.add_argument('-v', dest='verbose',
                        action='store_true',
                        help='''flag for verbosity''')

    ##get at the arguments
    args = parser.parse_args()

    ## do something..
    main(args.files, args.outdir, args.n , args.l, args.get_id, args.f_config, args.verbose)

