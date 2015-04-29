#!/usr/bin/env python
# Author:  csiu
# Created: 2015-02-06
import argparse
from ConfigParser import SafeConfigParser
import sys
import os
import re
import math
from utils import get_value_from_keycolonvalue_list, ensure_dir

usage = """Feature extraction:
- miRNA proximity = (1 + d/1000)^-1
  where d = distance from putative_tss to pre-miRNA

EXAMPLE:
python2.7 mirna_proximity.py -i ../test/test.gff -o ../test_outdir/mirna_proximity.gff
"""
def distance_score(d):
    score = math.pow((1.0+float(d)/1000),-1)
    return score

def calculate_distance(peak_start, peak_stop,
                  mirna_start, mirna_stop,
                  strand):
    if strand == '+':
        peak_start  = int(peak_start)
        mirna_start = int(mirna_start)
        d = mirna_start - peak_start
    else:
        peak_stop  = int(peak_stop)
        mirna_stop = int(mirna_stop)
        d = peak_stop - mirna_stop
    return d

def _is_chr_in_dataline(somefile):
    with open(somefile) as f:
        for l in f:
            if l.startswith('#'):
                continue

            if l.startswith('chr'):
                return True
            else:
                return False

def _validate_mirbase(f_mirna):
    # check that chr is stripped off f_mirna
    if not _is_chr_in_dataline(f_mirna):
        return f_mirna

    f_out = f_mirna + '.strip_chr.gff2'
    if os.path.exists(f_out):
        if not _is_chr_in_dataline(f_out):
            return f_out

    with open(f_out, 'w') as out:
        with open(f_mirna) as f:
            for l in f:
                if l.startswith('#'):
                    out.write(l)
                else:
                    out.write(re.sub('^chr', '', l))
    return f_out

def _swap_columns(f_cage, f_out):
    with open(f_out, 'w') as out:
        with open(f_cage) as f:
            for l in f:
                l = l.strip().split('\t')
                tss_up   = l[3]
                tss_down = l[4]
                info = l[8].split(';')

                start = get_value_from_keycolonvalue_list('start', info)
                stop  = get_value_from_keycolonvalue_list('stop', info)

                ## new info column
                info = filter(lambda x: not (x.startswith('start:') or x.startswith('stop:')), info)
                info.append('region_start:%s;region_stop:%s' % (tss_up, tss_down))
                l[8] = ';'.join(info)

                ## new start & stop
                l[3] = start
                l[4] = stop

                out.write('\t'.join(l) + '\n')
    return

def _make_newline(l, d):
    ## output in gff format
    mirna_proximity = str(distance_score(d))

    chrom      = l[0]
    peak_start = l[3]
    peak_stop  = l[4]
    strand     = l[6]

    info = l[8].split(';')
    region_up   = get_value_from_keycolonvalue_list('region_start', info)
    region_down = get_value_from_keycolonvalue_list('region_stop', info)

    mirna_start = l[12]
    mirna_stop  = l[13]

    mirna_info  = re.sub(' |"', '', l[17]).strip().split(';')
    mirna_acc   = get_value_from_keycolonvalue_list('ACC', mirna_info, '=')
    mirbase_id  = get_value_from_keycolonvalue_list('ID', mirna_info, '=')


    new_info = ';'.join(['distance:'+str(d),
                         #'region_start:'+region_up, 'region_stop:'+region_down,
                         'mirna_acc:'+mirna_acc, 'mirbase_id:'+mirbase_id,
                         'mirna_start:'+mirna_start, 'mirna_stop:'+mirna_stop])

    newline = '\t'.join([chrom, l[1], l[2], #'overlap', 'putative_tss',
                         peak_start, peak_stop, mirna_proximity, strand, '.',
                         new_info])
    return newline

def mirna_proximity(fi_cage, f_mirna, f_out, swapcol=True):
    ## check that chr is stripped off f_mirna
    f_mirna = _validate_mirbase(f_mirna)

    if swapcol:
      ## start & stop columns represent the +/-500bp region around TSS
      ## made it so that these columns represent the actual TSS range
      f_cage = f_out  + '.actualTSSregion.gff'
      _swap_columns(fi_cage, f_cage)
    else:
      f_cage = fi_cage

    ## find TSS associated with the 50kbp upstream of each miRNA
    f_overlap = f_out + 'cage_50kbpUpmirna_overlap.tmp'
    cmd = 'bedtools window -a '+f_cage+' -b '+f_mirna+' -l 0 -r 50000 -sw -sm > '+f_overlap
    os.system(cmd)

    ## calculate proximity score & write output
    with open(f_out, 'w') as out:
        with open(f_overlap) as f:
            for l in f:
                l = l.split('\t')
                strand = l[6]

                d = calculate_distance(l[3], l[4], l[12], l[13], strand)

                out.write( _make_newline(l, d) + '\n')

    ## remove intermediates
    os.system('rm '+f_overlap)
    os.system('rm %s' % f_cage)
    return


def main(f_cage, f_mirna, f_out):
    f_out = os.path.abspath(f_out)
    ensure_dir(f_out)

    mirna_proximity(f_cage, f_mirna, f_out)
    return f_out

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=usage,
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-i', '--infile', dest='infile',
                        required=True,
                        help='path to input file (same as the one used in features.py)')

    parser.add_argument('-o', '--outfile', dest='outfile',
                        default='mirna_proximity.gff',
                        help='''specify path to outfile
Default = "./mirna_proximity.gff"''')

    ## miRBase files
    parser.add_argument('--gff2', dest='gff2',
                        default='../data/miRBase/v20/hsa.gff2',
                        help='miRBase hsa.gff2 file')

    ##get at the arguments
    args = parser.parse_args()

    ## do something..
    main(args.infile, args.gff2, args.outfile)
