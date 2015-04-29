#!/usr/bin/env python
# Author:  csiu
# Created: 2015-04-22
import argparse
from utils import get_value_from_keycolonvalue_list, ensure_dir
import re

usage = """
- Adds 'intragenic'/'intergenic' label info for mirna

- label 'NA' is when mirna_start and mirna_stop is
  not found in the info column
- label 'unknown' is when mirna_start and mirna_stop is
  found in the info column, but the mirna is
  not found in the annotation file
"""
def _labelfile2dict_matchPos(labelfile):
    label_dict = {}
    with open(labelfile) as f:
        for l in f:
            l = l.split('\t')
            chrom  = l[0]
            start  = l[3]
            stop   = l[4]
            strand = l[6]
            label  = l[5]
            label_dict[','.join([chrom,start,stop,strand])] = label
    return label_dict

def _labelfile2dict_matchMirna(labelfile):
    label_dict = {}
    with open(labelfile) as f:
        for l in f:
            l = l.split('\t')
            label = l[5]

            mirna = l[1].lower()
            mirna = re.match('^(\w*-\w*-\d*)', mirna).group(1)
            label_dict[mirna] = label
    return label_dict

def main(infile, labelfile, outfile):
    label_dict = _labelfile2dict_matchMirna(labelfile)

    with open(outfile, 'w') as out:
        with open(infile) as f:
            for line in f:
                line = line.strip().split('\t')
                chrom  = line[0]
                strand = line[6]

                info   = line[8].strip(';')
                info   = re.split('[;@]', info)

                mstart = get_value_from_keycolonvalue_list('mirna_start', info)
                mstop  = get_value_from_keycolonvalue_list('mirna_stop', info)

                #mirna = ','.join([chrom,mstart,mstop,strand])

                if mstart == '' and mstop == '':
                    label = 'NA'
                else:
                    mirna = get_value_from_keycolonvalue_list('mirbase_id', info)

                    mirna = re.match('^(\w*-\w*-\d*)', mirna).group(1)
                    if label_dict.has_key(mirna):
                        label = label_dict[mirna]
                    else:
                        label = 'unknown'

                info.append('mirna_label:'+label)
                line[8] = ';'.join(info)
                newline = '\t'.join(line)
                out.write(newline + '\n')
    print outfile

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=usage,
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-i', '--infile', dest='infile',
                        required=True,
                        help='''path to <features.gff> file;
- tab separated
- Expects 'mirna_start:<#>;mirna_stop:<#>' to be in col 9
- Expects chr to be in column 1
- Expects strand to be in column 7
''')

    parser.add_argument('-l', dest='labelfile',
                        default='../data/miriad_human_labels_v2014.gff',
                        help='''path to <labels>.gff
- tab separated
- expects chr, mirna start, mirna stop, strand, and label
  to be columns 1, 4, 5, 7, 6 respectively
''')

    parser.add_argument('-o', '--outfile', dest='outfile',
                        help='''specify path to output file''')

    ##get at the arguments
    args = parser.parse_args()

    if args.outfile == None:
        outfile = '%s.label' % args.infile
    else:
        outfile = args.outfile

    ## do something..
    main(args.infile, args.labelfile, outfile)
