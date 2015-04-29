#!/usr/bin/env python
# Author:  csiu
# Created: 2015-04-19
import argparse
import sys
import os
import glob
import operator
from utils import random_string, get_value_from_keycolonvalue_list

usage = """Summarize samples from promi2 results
- considers tss-mirna pairing
- promoter in at least 1 sample
- take top promoter probability as representative probability for position

EXAMPLE:
python summarize-results.py -f '../test/test-summary/Predictions.*' -o ../Testout-summary.gff -s

The '-s' options is to print only results where mirna_prox score is > 0
"""

def _pull_putative_prom(somefile, somelist, is_strict):
    with open(somefile) as f:
        for l in f:
            l = l.strip().split('\t')
            label = l[13]

            if label.startswith('prom'):
                if is_strict:
                    mprox = float(get_value_from_keycolonvalue_list('mirna_prox',
                                                                    l[7].split(';')))
                    if mprox == 0:
                        continue

                chrom  = l[0]
                start  = l[3]
                stop   = l[4]
                count  = l[5]
                strand = l[6]
                prob_prom = l[11]

                info = l[8].split('@')[1].split(';')
                mirna = '%s-%s' % (get_value_from_keycolonvalue_list('mirna_start',
                                                                     info),
                                   get_value_from_keycolonvalue_list('mirna_stop',
                                                                     info))
                if mirna == '-':
                    mirna = 'NA-NA'

                try:
                    info = l[8].split('@')[2].split(';')
                    mirnaid = get_value_from_keycolonvalue_list('mirbase_id', info)
                    if mirnaid != '':
                        mirna = '%s:%s' % (mirna, mirnaid)
                    else:
                        mirna = '%s:%s' % (mirna, 'NA')
                except:
                    pass

                position = '%s,%s,%s,%s,%s' % (chrom, start, stop, strand, mirna)

                ## consider position with max prob
                if position in somelist:
                    old_count, old_prob, nlib, _ = somelist[position]
                    if (old_prob > prob_prom) or \
                           (old_prob == prob_prom and old_count >= count):
                        somelist[position][2] += 1
                        continue
                    else:
                        somelist[position] = [count, prob_prom, nlib+1, l[7]]
                else:
                    somelist[position] = [count, prob_prom, 1, l[7]]
    return somelist

def main(files, outfile, is_strict):
    listoffiles = glob.glob(files)

    putative_tss = {}
    for i in range(0, len(listoffiles)):
        putative_tss = _pull_putative_prom(listoffiles[i], putative_tss, is_strict)

    with open(outfile, 'w') as out:
        for k,v in putative_tss.iteritems():
            chrom, start, stop, strand, mirna = k.split(',')

            count, prob_prom, nlib, features = v

            ## writing out
            newline = '\t'.join([chrom, str(nlib), count, start, stop,
                                 prob_prom, strand, mirna,
                                 features])
            out.write(newline + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=usage,
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-f', '--files', dest='files',
                        required=True,
                        help='''path to input files.
These files should be the output of promi2.
You will need quotes when wildcard (*) is used.''')

    parser.add_argument('-s', dest='is_strict',
                        action="store_true",
                        help='''flag to consider only positions where
mirna_prox > 0''')

    parser.add_argument('-o', '--outfile', dest='outfile',
                        default='../promi2Summary_'+random_string(6)+'.gff',
                        help='''Specify path to output file. Columns are:
1. chrom
2. number of samples with position predicted as prom
3. count value associated with sample with highest promoter prob
4. start
5. stop
6. highest promoter probability
7. strand
8. mirna
9. {list of counts, list of probs}
''')

    ##get at the arguments
    args = parser.parse_args()

    ## do something..
    main(args.files, args.outfile, args.is_strict)
