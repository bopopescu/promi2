#!/usr/bin/env python
# Author:  csiu
# Created: 2015-02-17
import argparse
import sys
import os
from utils import random_string, get_value_from_keycolonvalue_list
usage = """Essentially running:
bedtools intersect -a <features.gff> -b <mirna_prox.gff> -s -f 1 -r -wao
"""
def bedtools_intersect(gff_a, gff_b, gff_out):
    ## unify
    cmd = 'bedtools intersect -a '+gff_a+' -b '+gff_b+' -s -f 1 -r -wao >'+gff_out
    print cmd
    os.system(cmd)
    return

def gff_unify_features(gff_a, gff_b, fname, dfvalue, f_out,
                       retainSourceFeature=False):
    ## unify
    f_out_tmp = f_out+'.tmp'
    bedtools_intersect(gff_a, gff_b, f_out_tmp)

    ## parse
    with open(f_out, 'w') as out:
        with open(f_out_tmp) as f:
            for l in f:
                l = l.strip().split('\t')

                chrom    = l[0]
                start    = l[3]
                stop     = l[4]
                count    = l[5]
                strand   = l[6]
                features = l[7]
                info_a   = l[8]
                _chrom   = l[9]

                if chrom == _chrom:
                    ## yes overlap of features w/ mirna_proximity
                    x_b      = l[14]
                    info_b   = l[17]
                    mirna_id = get_value_from_keycolonvalue_list('mirna_id',
                                                                 info_b.split(';'))
                else:
                    x_b        = dfvalue
                    info_b     = ''
                    mirna_id   = '.'

                features = '%s;%s:%s' % (features, fname, x_b)
                new_info = info_a + '@' + info_b

                if retainSourceFeature:
                    newline = '\t'.join([chrom, l[1], l[2], start, stop,
                                         count, strand, features, new_info])
                else:
                    newline = '\t'.join([chrom, 'putative_tss', mirna_id, start, stop,
                                         count, strand, features, new_info])
                out.write(newline + '\n')

    os.system('rm '+f_out_tmp)
    return

def _verify_mirbaseID(gff_infile, gff_outfile):
    with open(gff_outfile, 'w') as out:
        with open(gff_infile) as f:
            for l in f:
                info = l.strip().split('\t')[8].split('@')
                _x = info[-2].split(';')
                _y = info[-1].split(';')

                _x = get_value_from_keycolonvalue_list('mirbase_id', _x)
                _y = get_value_from_keycolonvalue_list('mirbase_id', _y)

                if _x == _y or _x == '' or _y == '':
                    out.write(l)
    return

def main(gff_a, gff_b, fname, dfvalue, f_out, retainSourceFeature=False):
    tmpfile = f_out + '.tmp'
    gff_unify_features(gff_a, gff_b, fname, dfvalue, tmpfile, retainSourceFeature)

    _verify_mirbaseID(tmpfile, f_out)
    os.remove(tmpfile)
    return f_out

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=usage,
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-a', dest='gff_a',
                        required=True,
                        help='''path to <features>.gff (this should contain
"cpg, cons, and tata" in the info column e.g. column 8)''')

    parser.add_argument('-b', dest='gff_b',
                        required=True,
                        help='''path to <mirna_proximity>.gff (this should contain
the feature you want in add to <GFF_A> in column 6 )''')

    parser.add_argument('-f', dest='fname',
                        required=True,
                        help='''name of feature''')

    parser.add_argument('-d', dest='dfvalue',
                        default='na',
                        help='''default feature value''')

    parser.add_argument('-o', dest='outfile',
                        default='all_features.gff',
                        help='specify outfile; default = "all_features.gff"')

    ##get at the arguments
    args = parser.parse_args()

    ## do something..
    main(args.gff_a, args.gff_b, args.fname, args.dfvalue, args.outfile)
