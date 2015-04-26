#!/usr/bin/env python
# Author:  csiu
# Created: 2015-02-02
import argparse
from ConfigParser import SafeConfigParser
import sys
import os

from utils import get_value_from_keycolonvalue_list, ensure_dir, random_string
import features
import mirna_proximity
import correlation
import gff_unify_features
import promirna

usage = """Runs promi2

Example:
python2.7 promi2.py -i ../test/test.gff -o ../Testout-promi2
"""
def _read_params(f_param):
    params_dict = {}
    with open(f_param) as f:
        for l in f:
            k,v = l.strip().split(':')
            params_dict[k] = float(v)

    mu1 = params_dict['mu_promoter']
    mu2 = params_dict['mu_background']
    lambda1 = params_dict['lambda_promoter']
    lambda2 = params_dict['lambda_background']

    betas = [i for i in params_dict.keys() if i.startswith('beta')]
    betas.sort()

    betas = [params_dict[b] for b in betas]

    return (mu1, mu2, lambda1, lambda2, betas)

def _make_prediction(prior_prom, p_prom, p_back):
    if str(prior_prom).endswith('*'):
        note = '*'
    else:
        note = ''
    if p_prom >= p_back:
        prediction = 'prom'+note
    else:
        prediction = 'back'+note

    return prediction

def promi2(f_param, listoffeatures, infile, outfile):
    mu1, mu2, lambda1, lambda2, betas = _read_params(f_param)

    if len(betas) != len(listoffeatures)+1:
        sys.exit("ERROR: number of betas does not match number of features")

    with open(outfile, 'w') as out:
        with open(infile) as f:
            for line in f:
                line = line.strip()

                l = line.split('\t')
                x = float(l[5])

                _features = l[7].split(';')

                fvalues = []
                for lof in listoffeatures:
                    try:
                        fvalues.append(float(get_value_from_keycolonvalue_list(lof, _features)))
                    except ValueError:
                        fvalues.append(0)

                p_prom, p_back, prior_prom, prior_back = promirna.promirna(x, mu1, mu2, lambda1, lambda2,
                                                                           betas, fvalues)
                prediction = _make_prediction(prior_prom, p_prom, p_back)

                #line = '\t'.join([line,
                #                  ';'.join(['prior_prom:'+str(prior_prom), 'prior_back:'+str(prior_back),
                #                            'prob_prom:'+str(p_prom), 'prob_back:'+str(p_back)]),
                #                  prediction]) + '\n'
                line = line + '\t%s\t%s\t%s\t%s\t%s\n' % (prior_prom, prior_back, p_prom, p_back, prediction)
                out.write(line)
    return

def _cleanup_extra_positions(infile, outfile):
    ## cleanup of extra positions
    ## compare miRNA positions in PROX & CORR
    with open(outfile, 'w') as out:
        with open(infile) as f:
            for line in f:
                l = line.split('\t')
                descript = l[8].split('@')

                if (descript[1] != '') and (descript[2] != '\n'):
                    info_mprox = descript[1].split(';')
                    prox_start = get_value_from_keycolonvalue_list('mirna_start', info_mprox)
                    prox_stop  = get_value_from_keycolonvalue_list('mirna_stop', info_mprox)

                    info_corr = descript[2].split(';')
                    corr_start = get_value_from_keycolonvalue_list('mirna_start', info_corr)
                    corr_stop  = get_value_from_keycolonvalue_list('mirna_stop', info_corr)

                    if (prox_start == corr_start) and \
                           (prox_stop == prox_stop):
                        out.write(line)
                else:
                    out.write(line)
    return outfile

def main(f_config, gff_cage, is_gff, outdir):
    cparser = SafeConfigParser()
    cparser.read(f_config)

    in_bname = os.path.basename(gff_cage)

    if outdir == None:
        outdir = 'promi2_outdir_'+in_bname+'_'+random_string(6)
    ensure_dir(outdir, False)

    f_param        = cparser.get('promi2','params')
    listoffeatures = cparser.get('promi2','features')
    listoffeatures = listoffeatures.split(',')
    if 'corr' in listoffeatures:
        is_consider_corr = True
        corrmethod = cparser.get('correlation','corrmethod')
    else:
        is_consider_corr = False

    ## PART1: Feature extraction
    if not is_gff:
        ## feature extraction: cpg, cons, tata (features.py)
        outdir_seqfeatures = os.path.join(outdir, 'seqfeatures')
        ensure_dir(outdir_seqfeatures, False)

        gff_1kbfeatures = os.path.join(outdir_seqfeatures, 'features_1kbseq.gff')

        f_fasta      = cparser.get('genome','fasta')
        f_chromsizes = cparser.get('genome','chromsizes')
        d_phastcons  = cparser.get('cons','phastcons')
        TRAP         = cparser.get('tata','trap')
        f_psemmatrix = cparser.get('tata','psem')

        features.main(gff_cage, outdir_seqfeatures,
                      f_fasta, f_chromsizes, d_phastcons, TRAP, f_psemmatrix,
                      gff_1kbfeatures)

        ## feature extraction: mirna_proximity (mirna_proximity.py)
        outdir_mprox = os.path.join(outdir, 'mprox')
        ensure_dir(outdir_mprox, False)

        gff_mirnaprox = os.path.join(outdir_mprox, 'features_mirnaprox.gff')

        gff_mirna     = cparser.get('mirbase','gff2')

        mirna_proximity.main(gff_cage, gff_mirna, gff_mirnaprox)

        ## merge extracted features (gff_unify_features.py)
        gff_features = os.path.join(outdir, 'Features.1kb.mprox.'+in_bname)
        gff_unify_features.main(gff_1kbfeatures, gff_mirnaprox, 'mirna_prox', '0', gff_features)

        if is_consider_corr:
            ## merge extracted features (gff_unify_features.py) after compute correlation
            gff_features_corr = os.path.join(outdir,
                                             'Features.1kb.mprox.%s.%s' % (corrmethod, in_bname))

            outdir_corr = os.path.join(outdir, 'corr')

            m_mirna = cparser.get('correlation', 'srnaseqmatrix')
            m_tss   = cparser.get('correlation', 'cageseqmatrix')

            gff_corr = correlation.main(gff_mirna, m_mirna, m_tss, corrmethod, outdir_corr)
            gff_unify_features.main(gff_features, gff_corr, 'corr', '0', gff_features_corr)

            gff_allfeatures = gff_features_corr
        else:
            gff_allfeatures = gff_features
    else:
        gff_allfeatures = gff_cage
        with open(gff_allfeatures) as f:
            l = f.readline().split('\t')
            if not (':' in l[7]):
                sys.exit('ERROR: this is not a features.gff formatted file')

    ## PART2: extract parameters & run promirna
    f_prediction = os.path.join(outdir, 'Predictions.'+in_bname+'.txt')
    print 'COMPUTING: "%s"...' % f_prediction
    promi2(f_param, listoffeatures, gff_allfeatures, f_prediction)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=usage,
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-i', '--infile', dest='infile',
                        required=True,
                        help='''path to input gff input file.
Tab-separated columns should be like:
  1. chrom
  2. source
  3. feature
  4. start (+500)
  5. stop (-500)
  6. normalized tag count
  7. strand
  8. .
  9. info
''')
    parser.add_argument('-f', dest='is_gff',
                        action='store_true',
                        help='flag to specify that infile is already features.gff file')

    parser.add_argument('-c', '--config', dest='f_config',
                        default='config.ini',
                        help='path to config file; default="config.ini"')


    parser.add_argument('-o', '--outdir', dest='outdir',
                        help='specify output directory')

    ##get at the arguments
    args = parser.parse_args()

    ## do something..
    main(args.f_config, args.infile, args.is_gff, args.outdir)
