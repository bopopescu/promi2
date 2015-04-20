#!/usr/bin/env python
# Author:  csiu
# Created: 2015-02-01
import argparse
from ConfigParser import SafeConfigParser
import sys
import os
import math
import re
import linecache
from utils import lmean, random_string, line_count, all_same, ensure_dir, get_value_from_keycolonvalue_list

usage = """Feature extraction:
- CpG content
- Conservation
- TATA affinity

EXAMPLE:
python2.7 features.py -i ../test/test.gff  -o ../test_outdir
"""
def _CpG_content(seq):
    L = len(seq)
    CpG = seq.count('CG')+seq.count('cg')+seq.count('Cg')+seq.count('cG')
    C = seq.count('C')+seq.count('c')
    G = seq.count('G')+seq.count('g')
    content = (float(CpG)/float(L)) / pow(float(C+G)/float(2*L),2)
    return content

def _conservation(chromosomes, d_split, d_cons, phastcons_dir):
    [os.makedirs(d) for d in [d_cons] if not os.path.exists(d)]

    for c in chromosomes:
        gff = os.path.join(d_split, "tss_filtered_all_"+c+".gff")
        wig = os.path.join(phastcons_dir, "track_000.chr"+c+".wib")
        out = os.path.join(d_cons, "conservation_all_"+c+".txt")

        d = os.path.dirname(os.path.realpath(__file__))
        extract_phastcons2_pl = os.path.join(d, "external/extract_phastcons2.pl")
        os.system(extract_phastcons2_pl+" "+wig+" "+gff+" > "+out)

def _average_conservation(f_cons,f_aver_cons):
    with open(f_aver_cons, 'w') as out:
        with open(f_cons) as f:
            for l in f:
                l = l.split(',')

                info = l[0].split(';')
                start_pos = get_value_from_keycolonvalue_list('start', info)
                stop_pos  = get_value_from_keycolonvalue_list('stop', info)

                try:
                    av_conservation = lmean([float(i) for i in l[1:]])
                except:
                    av_conservation = 0.0

                out.write('\t'.join([start_pos, stop_pos, str(av_conservation)]) +'\n')

### =====================================================================================
def prep_work(in_fa, in_bed, out_fa):
    ## bed to fasta format
    cmd = 'fastaFromBed -fi '+in_fa+' -bed '+in_bed+' -fo '+out_fa
    print cmd
    os.system(cmd)

    ## sort input
    sorted_gff = os.path.join(os.path.dirname(out_fa), os.path.basename(in_bed)+'.sorted.tmp')
    cmd = "sort -k4,5 -n "+in_bed+" > "+sorted_gff
    os.system(cmd)
    return sorted_gff

def cpg_content(f_fasta, out_cpg):
    ## f_fasta = filtered fasta file containing sequence of interest

    f_o = open(out_cpg, "w")
    f   = open(f_fasta).read()
    f_list = f.split('>')
    for item in f_list[1:]:
        item_list = item.split('\n')
        seq = "".join(([i for i in item_list[1:]]))
        chromosome = item_list[0].split(':')[0]; #print chromosome
        start_pos  = item_list[0].split(':')[1].split('-')[0]; #print start_pos
        stop_pos   = item_list[0].split(':')[1].split('-')[1]; #print stop_pos
        try:
            CpG_value = _CpG_content(seq)
        except:
            CpG_value=0.0
        print >> f_o, "\t".join((chromosome, str(int(start_pos)+1), stop_pos, str(CpG_value)))
    f_o.close()

    ## sort CpG
    sorted_cpg = out_cpg + '.sorted.tmp'
    cmd = "sort -k2,3 -n "+out_cpg+" > "+sorted_cpg
    os.system(cmd)
    return sorted_cpg

def conservation_score(f_chromsizes, d_phastcons, in_gff, out_avcons):
    tmp = random_string(12)

    d_split = out_avcons + 'gff_by_chromosome_'+tmp
    d_cons =  out_avcons + 'conservation_'+tmp

    [os.makedirs(d) for d in [d_split, d_cons] if not os.path.exists(d)]

    f_cons      = out_avcons + 'conservation.txt'
    f_aver_cons = out_avcons

    ## get chromosomes
    chromosomes = []
    with open(f_chromsizes) as f:
        for l in f:
            c = l.split('\t')[0]
            if (('random' not in c) and ('chrM' not in c) and ('chrUn' not in c)):
                chromosomes.append(c[3:])

    ## separate infile by chromosome
    for c in chromosomes:
        f_out = os.path.join(d_split, 'tss_filtered_all_'+c+'.gff')
        with open(f_out, 'w') as out:
            with open(in_gff) as f:
                for line in f:
                    chrom = line.split('\t')[0]
                    if (chrom == c):
                        out.write(line)

    ## calculate conservation per chromosome
    _conservation(chromosomes, d_split, d_cons, d_phastcons)

    ## merge chromosomes
    os.system("cat "+d_cons+"/conservation_all_*txt > "+f_cons)

    ## get average conservation
    _average_conservation(f_cons, f_aver_cons)

    ## cleanup
    is_same = []
    for c in chromosomes:
        n_gff = line_count(os.path.join(d_split, 'tss_filtered_all_'+c+'.gff'))
        n_con = line_count(os.path.join(d_cons,  'conservation_all_'+c+'.txt'))
        is_same.append(n_gff == n_con)
    if all_same(is_same):
        os.system('rm -r %s %s %s' % (d_split, d_cons, f_cons))
    else:
        not_equal = [chromosomes[i] for i,v in enumerate(is_same) if not v]
        sys.exit('Error: Total number of positions does not match for chr: ' + ' '.join(not_equal))

    ## sort average conservation
    sorted_avcons = out_avcons + '.sorted.tmp'
    cmd = "sort -k1,2 -n "+out_avcons+" > "+sorted_avcons
    os.system(cmd)
    return sorted_avcons

def tata_affinity(TRAP, f_psemmatrix, f_fasta, out_tata):
    cmd_trap = TRAP+" -s "+f_fasta+" --psem "+f_psemmatrix+" -g 0.5"+" -o "+out_tata
    print cmd_trap
    os.system(cmd_trap)

    ## sort TATA affinity
    sorted_tata = out_tata + '.sorted.tmp'
    tmp_tata = sorted_tata + random_string(12)
    with open(tmp_tata, 'w') as out:
        with open(out_tata) as f:
            for l in f:
                if not l.startswith("#"):
                    l = l.split('\t')
                    chr   = l[0].split(':')[0]
                    start = l[0].split(':')[1].split('-')[0]
                    stop  = l[0].split(':')[1].split('-')[1]

                    l[0] = '\t'.join([chr, start, stop])
                    out.write('\t'.join(l))

    cmd = "sort -k2,3 -n "+tmp_tata+" > "+sorted_tata
    os.system(cmd)
    os.system("rm "+tmp_tata)
    return sorted_tata

def build_features_matrix(sorted_gff, sorted_cpg, sorted_avcons, sorted_tata, f_out):
    ## check that all in files contain same number of data lines
    n_g = line_count(sorted_gff)
    n_c = line_count(sorted_cpg)
    n_a = line_count(sorted_avcons)
    n_t = line_count(sorted_tata)
    if not all_same([n_g, n_c, n_a, n_t]):
        sys.exit('Error: line count of feature files are not all equal:%s,%s,%s,%s' %
            n_g, n_c, n_a, n_t)

    ## create matrix
    lcount = 0
    with open(f_out, 'w') as out:
        with open(sorted_gff) as f:
            for l in f:
                lcount += 1

                l = l.strip().split('\t')
                c      = l[0]
                region_up   = l[3] #500bp   upstream of start; not used
                region_down = l[4] #500bp downstream of start; not used
                count  = l[5]
                strand = l[6]

                info = l[8].split(';')
                #dist_score = '?'

                peak_start = get_value_from_keycolonvalue_list('start', info)
                peak_stop  = get_value_from_keycolonvalue_list('stop', info)

                CpG_value    = linecache.getline(sorted_cpg,lcount).strip().split('\t')[3]
                try:
                    conservation = linecache.getline(sorted_avcons,lcount).strip().split('\t')[2]
                except:
                    conservation = '0'

                affinity     = linecache.getline(sorted_tata,lcount).strip().split('\t')[7]

                features = ';'.join(['cpg:'+CpG_value, 'cons:'+conservation, 'tata:'+affinity])
                new_info = ';'.join(['region_start:'+region_up, 'region_stop:'+region_down])
                line = '\t'.join([c, l[1], l[2],
                                  peak_start, peak_stop, count, strand,
                                  features, new_info])
                out.write(line + '\n')

### =====================================================================================
def main(infile, outdir,
    f_fasta, f_chromsizes, d_phastcons, TRAP, f_psemmatrix,
    fo_outfile):
    [os.makedirs(d) for d in [outdir] if not os.path.exists(d)]

    id_infile = re.sub('.gff$', '' ,os.path.basename(infile))

    fo_filtered_fasta = os.path.join(outdir, id_infile+'.fa')
    fo_cpg            = os.path.join(outdir, id_infile+'.cpg')
    fo_avcons         = os.path.join(outdir, id_infile+'.avgcons')
    fo_tata           = os.path.join(outdir, id_infile+'.tata')

    sorted_gff = prep_work(f_fasta, infile, fo_filtered_fasta)
    sorted_cpg = cpg_content(fo_filtered_fasta, fo_cpg)
    sorted_avcons = conservation_score(f_chromsizes, d_phastcons, infile, fo_avcons)
    sorted_tata = tata_affinity(TRAP, f_psemmatrix, fo_filtered_fasta, fo_tata)

    build_features_matrix(sorted_gff, sorted_cpg, sorted_avcons, sorted_tata, fo_outfile)

if __name__ == '__main__':
    cparser = SafeConfigParser()
    cparser.read('config.ini')

    parser = argparse.ArgumentParser(description=usage,
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-i', '--infile', dest='infile',
                        required=True,
                        help='path to input file')

    parser.add_argument('-o', '--outdir', dest='outdir',
                        default=".",
                        help='''specify path to output directory.
Default is the current working directory.''')

    parser.add_argument('-f', '--fo_outfile', dest='fo_outfile',
                        help='''specify path to final output file.''')

    ## For prep work
    parser.add_argument('--f_fasta', dest='f_fasta',
                        default=cparser.get('genome','fasta'),
                        help='''Path to (hg19) fasta file''')

    ## For conservation score
    parser.add_argument('--f_chromsizes', dest='f_chromsizes',
                        default=cparser.get('genome','chromsizes'),
                        help='''Path to (hg19) chrom sizes files.
Format: <chr>\t<chr.size>''')

    parser.add_argument('--d_phastcons', dest='d_phastcons',
                        default=cparser.get('cons','phastcons'),
                        help='''Path to directory containing "track_000.chr<#>.wib" files for Phastcons''')

    ## For TATA affinity
    parser.add_argument('--TRAP', dest='TRAP',
                        default=cparser.get('tata','trap'),
                        help='''Path to TRAP ("ANNOTATE_v3") binary''')

    parser.add_argument('--f_psemmatrix', dest='f_psemmatrix',
                        default=cparser.get('tata','psem'),
                        help='''Path to TATA_box_jaspar.psem''')


    ## get at the arguments
    args = parser.parse_args()

    if args.fo_outfile == None:
        fo_outfile = os.path.join(args.outdir, 'output.gff')
    else:
        ensure_dir(args.fo_outfile)
        fo_outfile = args.fo_outfile

    ## do something..
    main(args.infile, args.outdir,
        args.f_fasta, args.f_chromsizes, args.d_phastcons, args.TRAP, args.f_psemmatrix,
        fo_outfile)
