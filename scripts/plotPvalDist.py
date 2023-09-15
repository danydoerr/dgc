#!/usr/bin/env python3

from sys import stdout, stderr, exit
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF

import os
if not os.environ.get('DISPLAY', None):
    import matplotlib; matplotlib.use('Agg')
import matplotlib.pylab as plt


if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('pval_distribution', type=str, 
            help='File containng computed p-values from calculatePvalue')
    args = parser.parse_args()

    data = plt.array([float(x) for x in open(args.pval_distribution)])

    h0975 = len(data[data > 0.975])
    print('number of p-values larger 0.975: %s (%s%%)' %(h0975,
            h0975/float(len(data))), file=stderr)
    h0025 = len(data[data < 0.025])
    print('number of p-values smaller 0.025: %s (%s%%)' %(h0025,
            h0025/float(len(data))), file=stderr)

    plt.figure()
    plt.hist(data, 20)

    plt.axvline(x=0.025, color='r', linestyle='--')
    plt.axvline(x=0.975, color='r', linestyle='--')
    plt.ylabel('count')
    plt.xlabel('p-value')
    plt.savefig(stdout, format='pdf')

