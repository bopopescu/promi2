#!/usr/bin/env python
# Author:  csiu
# Created: 2015-02-26
import argparse
from ConfigParser import SafeConfigParser
import sys
import os
import re
from utils import ensure_dir
from shutil import copyfile

usage = """RUNS tc_quantify.sh FROM ROBIN ANDERSSON
-----------------------------------------

The tc_quantify.sh takes six arguments:

1) TC file
   (.txt file downloaded from FANTOM, or
   a file with 'chr:start..stop,strand' info from the first column)
2) A file with library ids
   (not really important that it is the correct ids,
   you could use a unique integer per row,
   a remnant from Robin's usage)
3) A file with the files (ctss.bed files) to be used
   including their paths
4) A file with the tag counts per library
5) A file with normalization factors for RLE normalization
6) out path

-----------------------------------------
Example:
python2.7 tc_normalization.py -i ../test/test-tc_quantify/f1_pos.txt -o ../Testout-tc
"""
def run_tc_quantify(f1, f2, f3, f4, f5, f6):
    d = os.path.dirname(os.path.realpath(__file__))
    tc_quantify = os.path.join(d, 'external/tc_quantify/bin/tc_quantify.sh')

    cmd = ' '.join([tc_quantify, f1, f2, f3, f4, f5, f6])
    os.system(cmd)

    copyfile(f3, os.path.join(f6, 'files.txt'))

    id = re.sub('.txt', '', os.path.basename(f1))
    #fo_count = os.path.join(f6, id + '_expression_max_count.bed')
    fo_tpm   = os.path.join(f6, id + '_expression_max_tpm.bed')
    return fo_tpm

def main(f_config, f1_pos, f6_outdir):
    f6_outdir = f6_outdir + '/'
    ensure_dir(f6_outdir, False)

    cparser = SafeConfigParser()
    cparser.read(f_config)

    f2_libs        = cparser.get('tc_normalization', 'ids')
    f3_files       = cparser.get('tc_normalization', 'files')
    f4_counts      = cparser.get('tc_normalization', 'counts')
    f5_normfactors = cparser.get('tc_normalization', 'normfactors')

    fo_tpm = run_tc_quantify(f1_pos,
                             f2_libs, f3_files, f4_counts, f5_normfactors,
                             f6_outdir)

    return fo_tpm

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=usage,
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-i', '--infile', dest='infile',
                        required=True,
                        help='''path to input file.
First column of tab file should be in format:
"chr<chr>:<start>..<stop>,<strand>"''')

    parser.add_argument('-o', '--outdir', dest='outdir',
                        required=True,
                        help='''specify out path/outdir''')

    parser.add_argument('-c', '--config', dest='config',
                        default='tcnorm.ini',
                        help='''path to tcnorm.ini configuration file''')

    ##get at the arguments
    args = parser.parse_args()

    ## do something..
    main(args.config, args.infile, args.outdir)
