#!/usr/bin/env python
# Author:  csiu
# Created: 2015-02-12
import argparse
from ConfigParser import SafeConfigParser
import sys
import os
import re
from utils import random_string, get_value_from_keycolonvalue_list, ensure_dir, random_string
import linecache
from scipy.stats.stats import pearsonr, spearmanr

usage = """Correlation feature extraction

!!! This is a custom script !!!

1. find miRNA positions from sRNAseq in miRBase   -> mirna.gff
2. transform CAGE tss ids/positions to gff format -> tss.gff
3. find overlap pairs between CAGE(tss.gff) & 50kb upstream of sRNA-seq(mirna.gff)
    -> pair-pos_tss-overlap-50kbUpmirna.txt
4. find sample pairs between CAGE and sRNA-seq
   (e.g. match columns in CAGE and sRNA-seq matrices)
    -> pair-sample_cage-srnaseq.txt
5. compute correlation for each TSS-miRNA pair
   given the column pairings
    -> Correlation.{pearson,spearman}.gff
"""
def correlation_pearson(X,Y):
    ## Calculates a Pearson correlation coefficient
    ## and the p-value for testing non-correlation.
    ## http://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.pearsonr.html#scipy.stats.pearsonr
    ## ... requires dataset be normally distributed
    return pearsonr(X,Y)

def correlation_spearman(a, b=None, axis=0):
    ## Calculates a Spearman rank-order correlation coefficient
    ## and the p-value to test for non-correlation.
    ## http://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.spearmanr.html#scipy.stats.spearmanr
    ## ... not assume X, Y are normally distributed
    return spearmanr(a, b, axis)


def _find_miRNA_pos(m_mirna, mirbase_gff2, fo_mirna_gff):
    mirna_dict = {}
    with open(mirbase_gff2) as f:
        for l in f:
            if l.startswith('#'):
                continue

            l = l.strip().split('\t')
            mirna = l[8].split('; ')[1].split('"')[1].lower()

            chrom  = re.sub('^chr', '', l[0])
            start  = l[3]
            stop   = l[4]
            strand = l[6]

            position = '%s,%s,%s,%s' % (chrom, start, stop, strand)
            mirna_dict[mirna] = position

    with open(fo_mirna_gff, 'w') as out:
        with open(m_mirna) as f:
            line = 0
            for l in f:
                line += 1
                if not l.startswith('hsa'):
                    continue

                mirna = l.split('\t')[0].lower()
                _mirna = re.sub('-[35]p$', '', mirna)

                try:
                    chrom, start, stop, strand = mirna_dict[_mirna].split(',')
                    newline = '\t'.join([chrom, 'mirbase', _mirna,
                                         start, stop, '.', strand, str(line),
                                         'srnaseq_id:'+mirna])
                    out.write(newline +'\n')
                except KeyError:
                    if not 'novelmir' in _mirna:
                        found = [s for s in mirna_dict.keys() if _mirna+'-' in s]
                        if found == []:
                            found = [s for s in mirna_dict.keys() \
                                     if s.startswith(_mirna) \
                                     if not s[len(_mirna)].isdigit()]

                        for m_id in found:
                            chrom, start, stop, strand = mirna_dict[m_id].split(',')
                            newline = '\t'.join([chrom, 'mirbase', m_id,
                                                 start, stop, '.', strand, str(line),
                                                 'srnaseq_id:'+mirna])
                            out.write(newline +'\n')
    return

def _get_tss_pos(m_tss, fo_tss_gff):
    with open(fo_tss_gff, 'w') as out:
        with open(m_tss) as f:
            line = 0
            for l in f:
                line += 1
                if l.startswith('#') or \
                   l.startswith('00Annotation') or \
                   l.startswith('01STAT'):
                    continue

                tid = re.sub('^chr', '', l.split('\t')[0]).strip()
                chrom, start, _, stop, strand = re.split('[:,.]', tid)

                newline = '\t'.join([chrom, 'CAGE', 'putative_tss',
                                     start, stop, '.', strand,
                                     str(line), '.'])
                out.write(newline +'\n')
    return

def _get_tss_mirna_pairings(f_tss_gff, f_mirna_gff, fo_pos_pairing):
    f_tmp = fo_pos_pairing + '.' + random_string(6)
    cmd = 'bedtools window -a '+f_tss_gff+' -b '+f_mirna_gff+' -l 0 -r 50000 -sw -sm > '+f_tmp
    os.system(cmd)

    with open(fo_pos_pairing, 'w') as out:
        with open(f_tmp) as f:
            for l in f:
                l = l.strip().split('\t')
                chrom     = l[0]
                tss_start = l[3]
                tss_stop  = l[4]
                strand    = l[6]
                cage_line = l[7]
                mirna_mirbase = l[11]
                mirna_start   = l[12]
                mirna_stop    = l[13]
                srnaseq_line  = l[16]

                tss_annot   = ','.join(['title='+'chr%s:%s..%s,%s' %
                                        (chrom, tss_start, tss_stop, strand)])

                mirna_annot = ','.join(['title='+l[17].split(':')[1],
                                        'mirbase_id='+mirna_mirbase,
                                        'mirna_start='+mirna_start,
                                        'mirna_stop='+mirna_stop])

                newline = '\t'.join([tss_annot, cage_line, mirna_annot, srnaseq_line])
                out.write(newline +'\n')
    os.system('rm '+f_tmp)
    return

def _index_cage(m_tss):
    cage_index    = {}
    cage_id_pattern    = re.compile('^tpm.*(CNhs.*\..*)$')
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
    return cage_index

def _index_srnaseq(m_mirna):
    srnaseq_index = {}
    srnaseq_id_pattern = re.compile('^.*(SRh.*?\..*?)\.')
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
    return srnaseq_index

def _get_sample_pairings(cage_index, srnaseq_index, fo_sample_pairing):
    with open(fo_sample_pairing, 'w') as out:
        sample_ids = set(cage_index.keys()).intersection(srnaseq_index.keys())
        for k in sample_ids:
            for c in cage_index[k]:
                for m in srnaseq_index[k]:
                    cid, cindex = c.split(':')
                    mid, mindex = m.split(':')
                    out.write('\t'.join([k, cindex, mindex,
                                         '%s,%s' % (cid, mid)]) +'\n')
    return

def _compute_correlation(f_pos_pairing, f_sample_pairing,
                         m_tss, m_mirna,
                         fo_corr, method="spearman", source="putative_tss"):
    ## sample matching
    cage_sample_order    = []
    srnaseq_sample_order = []
    with open(f_sample_pairing) as f:
        for l in f:
            _, cindex, mindex, _ = l.split('\t')
            cage_sample_order.append(int(cindex))
            srnaseq_sample_order.append(int(mindex))

    ## position matching
    with open(fo_corr, 'w') as out:
        with open(f_pos_pairing) as f:
            for l in f:
                tss_info, n_tss, mirna_info, n_mirna = l.strip().split('\t')

                line_tss   = linecache.getline(m_tss, int(n_tss)).strip().split('\t')
                line_mirna = linecache.getline(m_mirna, int(n_mirna)).strip().split('\t')

                X = [float(line_tss[i])   for i in cage_sample_order]
                Y = [float(line_mirna[i]) for i in srnaseq_sample_order]

                if method=="pearson":
                    corr, pval = correlation_pearson(X,Y)
                else:
                    corr, pval =  correlation_spearman(X,Y)

                if str(corr) == "nan":
                    ## Note: when X or Y is 0 -> corr == nan
                    if sum(X) == 0 or sum(Y) == 0:
                        continue

                ## writing output
                _, chrom, start, _, stop, strand = re.split('[r:.,]',tss_info)

                mirna_info = re.sub('=', ':', mirna_info)
                mirna_info = re.sub(',', ';', mirna_info)
                mirna_info = re.sub('^title', 'mirna_id', mirna_info)

                new_info = ';'.join(['pval:'+str(pval), mirna_info])
                newline  = '\t'.join([chrom, source, method,
                                      start, stop, str(corr), strand,
                                      '.', new_info])

                out.write(newline + '\n')

def main(mirbase_gff2, m_mirna, m_tss, method, outdir):
    ensure_dir(outdir, False)

    fo_mirna_gff = os.path.join(outdir, '4corr_mirna.gff')
    fo_tss_gff   = os.path.join(outdir, '4corr_tss.gff')
    fo_pos_pairing    = os.path.join(outdir, '4corr_pair-pos_tss-overlap-50kbUpmirna.txt')
    fo_sample_pairing = os.path.join(outdir, '4corr_pair-sample_cage-srnaseq.txt')
    fo_corr = os.path.join(outdir, 'featuresgff_correlation-' + method + '.gff')

    if os.path.exists(fo_pos_pairing) and os.path.exists(fo_sample_pairing):
        _compute_correlation(fo_pos_pairing, fo_sample_pairing, m_tss, m_mirna, fo_corr, method)
    else:
        ## to gff format
        _find_miRNA_pos(m_mirna, mirbase_gff2, fo_mirna_gff)
        _get_tss_pos(m_tss, fo_tss_gff)

        ## finding pairs
        _get_tss_mirna_pairings(fo_tss_gff, fo_mirna_gff, fo_pos_pairing)

        cage_index    = _index_cage(m_tss)
        srnaseq_index = _index_srnaseq(m_mirna)
        _get_sample_pairings(cage_index, srnaseq_index, fo_sample_pairing)


        ## computing correlation
        _compute_correlation(fo_pos_pairing, fo_sample_pairing, m_tss, m_mirna, fo_corr, method)

    return fo_corr

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=usage,
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-m', '--mirna', dest='f_mirna',
                        required=True,
                        help='path to small-RNA-seq matrix')

    parser.add_argument('-t', '--tss', dest='f_tss',
                        required=True,
                        help='path to CAGE matrix')

    parser.add_argument('-c', '--config', dest='config',
                        default='config.ini',
                        help='path config file')

    parser.add_argument('--method', dest='method',
                        default='spearman',
                        help='''Correlation method.
Choose 1 of {spearman, pearson}.
Default="spearman"''')

    parser.add_argument('-o', '--outdir', dest='outdir',
                        default='correlation_'+random_string(6),
                        help='specify output directory')

    ##get at the arguments
    args = parser.parse_args()

    cparser = SafeConfigParser()
    cparser.read(args.config)

    mirbase_gff2 = cparser.get('mirbase','gff2')

    ## do something..
    main(mirbase_gff2, args.f_mirna, args.f_tss, args.method, args.outdir)
