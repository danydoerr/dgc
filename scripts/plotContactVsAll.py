#!/usr/bin/env python3

from sys import stdout, stderr, exit
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF, \
        FileType
from os.path import basename
import logging
import csv

import os
if not os.environ.get('DISPLAY', None):
    import matplotlib; matplotlib.use('Agg')

import numpy as np
from matplotlib import pylab as plt

from hic import readHiCMapTRV, label2coords, PAT_COORD, assimilateMatrices
from identifyOutlier import readOutlier

LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)


def extractContactHist(mtrx, start, end, hist_size):

    tmtrx = mtrx.transpose()
    data = list()
    for i in xrange(start, end):
        diag = np.diagonal(tmtrx, i)
        diag = diag[diag.size/2:]
        data.extend(diag[diag != None])

    hist, bins = plt.histogram(data, hist_size)
    bin_centres = (bins[:-1] + bins[1:])/2

    return bin_centres, hist
    

if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('matrices', nargs='+', type=file, 
            help='matrices in TRV format')
    parser.add_argument('coordinate', type=str, 
            help='Genome coordinates of the regiono of interest, given in ' + \
                    'the format \'<chrid>:<start>-<end>\'. ')
    parser.add_argument('-s', '--hist_size', type=int, default=500,
            help='size of the histogram used in estimating Gaussian coefficients')
    parser.add_argument('-i', '--ignore', type=str, 
            help='ignore rows/columns specified in this file')
    args = parser.parse_args()


    # setup logging
    ch = logging.StreamHandler(stderr)
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(logging.Formatter('%(levelname)s\t%(asctime)s\t%(message)s'))
    LOG.addHandler(ch)

    ignoreRowNames = set()
    ignoreColNames = set()

    if args.ignore:
        ignoreRowNames, ignoreColNames  = readOutlier(open(args.ignore))


    hic_maps = list()
    for f in args.matrices:
        hic_maps.append(readHiCMapTRV(f, None, ignoreColNames, ignoreRowNames))

    ass_hic_maps = assimilateMatrices(hic_maps)
    m = PAT_COORD.match(args.coordinate)
    chrx, start, end = m.groups()

    coords = map(label2coords, hic_maps[0][1])

    i = coords.index((chrx, int(start)))
    j = coords.index((chrx, int(end)))

    plt.figure()
    for mtrx, rownames, colnames, axes_labels in ass_hic_maps:
        x, y = extractContactHist(mtrx, i, j+1, args.hist_size)
        plt.plot(x, y, label=' vs '.join(axes_labels))

    plt.legend()
    plt.savefig(stdout, format='pdf')
