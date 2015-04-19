#!/usr/bin/env python
# Author:  csiu
# Created: 2015-01-27

def inverse_gaussian_pdf(x, mu, lamb):
    import math
    import sys
    ## returns p(X_i | Z_k); ith region, kth class
    ## this is the 'likelihood' part
    _a = math.pow((lamb/(2*math.pi*math.pow(x,3))), 0.5)
    _b = math.exp((-lamb*math.pow((x-mu),2))/(2*math.pow(mu,2)*x))
    if _b == 0: _b = sys.float_info.min

    pdf = _a * _b
    return pdf

def prior(betas, features):
    import math

    features.insert(0, 'NA')
    if not (len(betas) == len(features)):
        sys.exit('Error: missing args')

    ## prior probability
    yi = betas[0]
    for i in range(1,len(betas)):
        yi += betas[i]*features[i]
    pi = 1.0/(1+math.exp(-yi))

    features.pop(0)
    return pi

def promirna(x,
             mu1, mu2, lambda1, lambda2,
             betas, features):
    ## class 1=promoter; 2=background

    ## prior probability
    pi1 = prior(betas, features)
    pi2 = 1.0 - pi1

    ## posterior probability
    _p_prom = (pi1*inverse_gaussian_pdf(x,mu1,lambda1))
    _p_back = (pi2*inverse_gaussian_pdf(x,mu2,lambda2))
    evidence = _p_prom + _p_back

    try:
        p_prom = _p_prom/evidence
        p_back = _p_back/evidence
    except ZeroDivisionError:
        p_prom = pi1
        p_back = pi2
        pi1 = '%s*' % pi1
        pi2 = '%s*' % pi2

    return (p_prom, p_back, pi1, pi2)
