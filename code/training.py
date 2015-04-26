#!/usr/bin/env python
# Author:  csiu
# Created: 2015-04-19
import argparse
from ConfigParser import SafeConfigParser
import sys
import os

import math
from promirna import prior
from promirna import inverse_gaussian_pdf as pdf
from utils import lmean, get_value_from_keycolonvalue_list

usage = """Assumption
----------
- background set will contain '*back*' in 3rd column of gff file

EXAMPLE:
python2.7 training.py -i ../Test-tset/TrainingSet.gff
"""
def _trainset_to_formatgff(infile, outfile):
    ## this was the original format:
    ## CAGE TSS <chr> <start> <stop> <strand> <norm.count> <cpg> <cons> <tata> <mirna.prox> <label>
    with open(outfile, 'w') as out:
        with open(infile) as f:
            for l in f:
                source, feature, chrom, start, stop, strand, ncount, cpg, cons, tata, prox, label = l.strip().split('\t')

                info = ';'.join(['cpg:'+cpg,
                                 'cons:'+cons,
                                 'tata:'+tata,
                                 'mirna_prox:'+prox])

                newline = '\t'.join([chrom,
                                     source,
                                     feature,
                                     start,
                                     stop,
                                     ncount,
                                     strand,
                                     info,
                                     'label:'+label])
                out.write(newline + '\n')

def estimate_betas(f_trainingset, trainingfeatures):
    if 'mirna_prox' in trainingfeatures:
        mprox = trainingfeatures.index('mirna_prox')
        trainingfeatures.pop(mprox)
        add_mprox = True
    else:
        add_mprox = False

    ## preprocess to right format for Rscript
    f_intermediate = f_trainingset + '.intermediate.tmp'
    with open(f_intermediate, 'w') as out:
        with open(f_trainingset) as f:
            for l in f:
                l = l.strip().split('\t')
                features = l[7].split(';')

                fvalues = []
                for i in trainingfeatures:
                    try:
                        fvalues.append(float(
                            get_value_from_keycolonvalue_list(i, features)))
                    except ValueError:
                        fvalues.append(0)

                fvalues = [str(i) for i in fvalues]

                if 'back' in l[2].lower():
                    label = '0'
                else:
                    label = '1'

                out.write('\t'.join([label] + fvalues) +'\n')

    ## estimating the beta parameters
    ## note: mirna_proximity is not considered in the fitting of betas
    f_beta_tmp = f_trainingset + '.parameters_beta.tmp'
    os.system('R --slave --vanilla --args '+f_intermediate+' '+f_beta_tmp+\
              ' < external/choose_beta_params.R')
    betas = []
    with open(f_beta_tmp) as f:
        for l in f:
            betas.append(l.strip())

    ## estimate beta4 (for mirna_proximity)
    betas = [float(b) for b in betas]
    if add_mprox:
        beta_mprox = min(betas[1:])
        betas.insert(mprox+1, beta_mprox)

    os.system('rm %s %s' % (f_intermediate, f_beta_tmp))
    return betas

def estimate_initial_params(f_trainingset):
    x1 = []
    x2 = []
    with open(f_trainingset) as f:
        for l in f:
            l = l.strip().split('\t')
            label = l[2].lower()
            count = float(l[5])

            if 'back' in label:
                x2.append(count)
            else:
                x1.append(count)

        mu1 = lmean(x1)
        mu2 = lmean(x2)

        sum1 = 0
        for i in x1:
            sum1 += 1.0/i - 1.0/mu1
        lambda1 = 1.0/ (sum1/float(len(x1)))

        sum2 = 0
        for i in x2:
            sum2 += 1.0/i - 1.0/mu2
        lambda2 = 1.0/ (sum2/float(len(x2)))

        params = [mu1, mu2, lambda1, lambda2]
        return params

def initialize_all_params(f_trainingset, trainingfeatures):
    betas  = estimate_betas(f_trainingset, trainingfeatures)
    params = estimate_initial_params(f_trainingset)
    return betas, params

def read_data(f_trainingset, trainingfeatures):
    data = []
    with open(f_trainingset) as f:
        for l in f:
            l = l.split('\t')
            count = float(l[5])

            if count != 0:
                chrom  = l[0]
                start  = int(l[3])
                stop   = int(l[4])
                strand = l[6]

                features = l[7].split(';')

                fvalues = []
                for i in trainingfeatures:
                    try:
                        fvalues.append(float(
                            get_value_from_keycolonvalue_list(i, features)))
                    except ValueError:
                        fvalues.append(0)

                label  = l[2].lower()
                ##probability of promoter (z1) & background (z2)
                if 'back' in label:
                    z1 = 0.0
                    z2 = 1.0
                else:
                    z1 = 0.5
                    z2 = 0.5

                item = (chrom, start, stop, strand, count,
                        trainingfeatures, fvalues,
                        z1, z2)
                data.append(item)
    return data

def E_step(data, params, betas):
    ## unpack
    mu1, mu2, lambda1, lambda2 = params

    new_dataset = []
    for n in range(0,len(data)):
        chrom, start, stop, strand, x, fids, fvalues, z1_0, z2_0 = data[n]

        ## DISCRETIZING CpG content
        try:
            cpg = fids.index('cpg')
            if fvalues[cpg]<0.45:
                fvalues[cpg] = 0.1
            else:
                fvalues[cpg] = 0.9
        except ValueError:
            pass

        f1 = pdf(x, mu1, lambda1)
        f2 = pdf(x, mu2, lambda2)
        if f1==0.0 and f2==0.0:
            z1 = 1.0
            z2 = 0.0

        pi1 = prior(betas, fvalues)
        pi2 = 1 - pi1

        ## resetting z1, "prob_promoter"
        if z1_0==0 or z1_0==1.0:
            z1 = z1_0
        else:
            z1 = (f1*pi1) / (f1*pi1 + f2*pi2)

        ## resetting z2, "prob_background"
        if z2_0==0 or z2_0==1.0:
            z2 = z2_0
        else:
            z2 = (f2*pi2) / (f1*pi1 + f2*pi2)

        item = (x, z1, z2,
                fvalues,
                chrom, start, stop, strand,
                f1, f2, #FIXME:is this needed?
                pi1, pi2)
        new_dataset.append(item)

    return new_dataset

def M_step(new_dataset):
    ## re-estimate the parameters based on the marginal probability distributions
    sum_num_1 = 0
    sum_num_2 = 0
    sum_den_1 = 0
    sum_den_2 = 0
    for n in range(0, len(new_dataset)):
        x  = new_dataset[n][0]
        z1 = new_dataset[n][1]
        z2 = new_dataset[n][2]

        sum_num_1 += x*z1
        sum_den_1 += z1

        sum_num_2 += x*z2
        sum_den_2 += z2

    mu1_new = sum_num_1/sum_den_1
    mu2_new = sum_num_2/sum_den_2


    sum_num_lambda1 = 0
    sum_den_lambda1 = 0
    sum_num_lambda2 = 0
    sum_den_lambda2 = 0
    for n in range(0, len(new_dataset)):
        x  = new_dataset[n][0]
        z1 = new_dataset[n][1]
        z2 = new_dataset[n][2]

        sum_num_lambda1 += z1
        sum_den_lambda1 += z1*( pow((x-mu1_new),2) /x )

        sum_num_lambda2 += z2
        sum_den_lambda2 += z2*( pow((x-mu2_new),2) /x )

    lambda1_new = (sum_num_lambda1/sum_den_lambda1) * pow(mu1_new,2)
    lambda2_new = (sum_num_lambda2/sum_den_lambda2) * pow(mu2_new,2)

    return mu1_new, mu2_new, lambda1_new, lambda2_new

def calculate_likelihood(new_dataset, new_params):
    mu1, mu2, lambda1, lambda2 = new_params

    L_total = 0.0
    for n in range(0, len(new_dataset)):
        x   = new_dataset[n][0]
        pi1 = new_dataset[n][10]
        pi2 = new_dataset[n][11]
        f1 = pdf(x, mu1, lambda1)
        f2 = pdf(x, mu2, lambda2)

        L_total += math.log( (pi1*f1)+(pi2*f2) )
    return L_total

def _make_intermediate_from_EM(new_dataset, f_out):
    with open(f_out, 'w') as out:
        out.write('\t'.join(['chromosome', 'start', 'stop', 'strand',
                             'counts',
                             'prior_probability',
                             'prob_promoter', 'prob_background'])+'\n')
        for i in new_dataset:
            out.write('\t'.join([str(i) for i in [i[7], i[8], i[9], i[10],
                                                  i[0],
                                                  i[13],
                                                  i[1], i[2]]]) + '\n')
    return f_out

def training(f_trainingset, trainingfeatures,
             fo_params, get_intermediates=False):
    data = read_data(f_trainingset, trainingfeatures)

    betas, params  = initialize_all_params(f_trainingset, trainingfeatures)

    likelihood_old = -1000000.0

    itr = 1
    print 'Iteration: %s' % itr
    new_dataset = E_step(data, params, betas)
    new_params  = M_step(new_dataset)
    likelihood_new = calculate_likelihood(new_dataset, new_params)

    while (likelihood_new - likelihood_old) > 0.1:
        itr += 1
        print 'Iteration: %s' % itr
        likelihood_old = likelihood_new

        new_dataset = E_step(data, new_params, betas)
        new_params  = M_step(new_dataset)
        likelihood_new = calculate_likelihood(new_dataset, new_params)

    with open(fo_params, 'w') as out:
        mu1_new, mu2_new, lambda1_new, lambda2_new = new_params
        out.write('mu_promoter:%s\nmu_background:%s\n' % (mu1_new, mu2_new))
        out.write('lambda_promoter:%s\nlambda_background:%s\n' % (lambda1_new, lambda2_new))
        for i,b in enumerate(betas):
            out.write('beta%s:%s\n' % (i, b))

    if get_intermediates:
        print _make_intermediate_from_EM(new_dataset, fo_params+'.intermediate')

    return fo_params


def main(f_trainingset, f_config, fo_params, is_preprocess, is_get_intermediate):
    cparser = SafeConfigParser()
    cparser.read(f_config)
    trainingfeatures = cparser.get('training','trainingfeatures').split(',')

    if is_preprocess:
        with open(f_trainingset) as f:
            l = f.readline()
            if len(l.split('\t')) != 12:
                sys.exit("ERROR: preprocessing ('-p') was selected, " + \
                         "but '%s' is not in the original format... see '-h'" % f_trainingset)

        tmp = f_trainingset + '.gff'
        _trainset_to_formatgff(f_trainingset, tmp)
        f_trainingset = tmp
    else:
        ## simple check if in right format
        with open(f_trainingset) as f:
            l = f.readline()
            l = l.split('\t')
            if not 'cpg' in l[7]:
                sys.exit("ERROR: '%s' is not in the right format... see '-h'" % f_trainingset)

    f_param = training(f_trainingset, trainingfeatures, fo_params, is_get_intermediate)
    print f_param
    return f_param

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=usage,
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-i', '--infile', dest='infile',
                        required=True,
                        help='''Path to input file (e.g. the trainingset) in gff format
The tab-separated columns are as follows:
 1. chrom
 2. source
 3. feature (should contain "back" if is background sample
 4. start
 5. stop
 6. normalized count
 7. strand
 8. features (e.g. "cpg:#;cons:#;tata:#")
 9. other
''')

    parser.add_argument('-p', dest='preprocess',
                        action='store_true',
                        help='''flag to preprocess training file (from original format) to gff format.
The tab-separated columns in the original format are as follows:
 1. CAGE
 2. TSS
 3. chrom
 4. start
 5. stop
 6. strand
 7. normalized count
 8. cpg
 9. conservation score
 10. tata binding affinity
 11. mirna proximity schore
 12. label
''')

    parser.add_argument('-c', '--config', dest='f_config',
                        default='config.ini',
                        help='path to "config.ini"')

    parser.add_argument('-o', '--outfile', dest='outfile',
                        help='specify path to outfile')



    parser.add_argument('-g', dest='get_intermediate',
                        action='store_true',
                        help='flag to get intermediate from EM')

    ##get at the arguments
    args = parser.parse_args()

    if args.outfile == None:
        outfile = args.infile + '.finalparams'
    else:
        outfile = args.outfile

    ## do something..
    main(args.infile, args.f_config, outfile, args.preprocess, args.get_intermediate)
