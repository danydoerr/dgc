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

from hic import readHiCMapTRV, label2coords, PAT_COORD 
from testGCinCandidates import readBioGCs
from identifyOutlier import readOutlier

LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)


def calculateMeans(mtrx):

    res = np.zeros((mtrx.shape[0], ))
    for i in xrange(mtrx.shape[0]):
        data = mtrx[i, :][mtrx[i, :] != None]
        if data.size:
            res[i] = data.mean() 
    return res

def prepareData(means, max_steps=1000):

    data = list()
    total_lengths = dict()
    for m, labels, name in means:
        srtd_labels, A = zip(*sorted(zip(labels, xrange(len(labels)))))
        data.append((m[np.array(A)], srtd_labels, name))
        for i in xrange(len(srtd_labels)-1):
            chrx, end = srtd_labels[i]
            if chrx != srtd_labels[i+1][0]:
                if not total_lengths.has_key(chrx):
                    total_lengths[chrx] = 0
                total_lengths[chrx] = max(total_lengths[chrx], end)
            chrx, end = srtd_labels[-1]
            if not total_lengths.has_key(chrx):
                total_lengths[chrx] = 0
            total_lengths[chrx] = max(total_lengths[chrx], end)

    f = float(max_steps)/sum(total_lengths.values())

    start_pos = dict()
    p = 0
    for c, l in sorted(total_lengths.items()):
        start_pos[c] = p
        p += l

    return data, start_pos, f

def plotDistributions(means, start_pos, f, gene_clusters, out):

    plt.figure()
    
    for name, (chrx, start, end) in gene_clusters:
        plt.axvspan(xmin=(start_pos[chrx]+start)*f,
                xmax=(start_pos[chrx]+end)*f, alpha=0.5, label=name)

    for m, labels, name in means:
        x = np.array(map(lambda x: (start_pos[x[0]] + x[1])*f, labels))
        plt.plot(x[0::20], m[0::20], label=name)

    isFirst = True
    for start in start_pos.values():
        if isFirst:
            plt.axvline(x=start*f, linestyle='dotted', color='black',
                    label='chromosome boundary')
            isFirst = False
        else:
            plt.axvline(x=start*f, linestyle='dotted', color='black')

    plt.legend()
    plt.savefig(out, format='pdf')


if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('matrices', nargs='+', type=file, 
            help='matrices in TRV format')
    parser.add_argument('gene_clusters', type=file, 
            help='file containing coordinates of biological gene clusters')
    parser.add_argument('-i', '--ignore', nargs='+', type=file, 
            help='ignore rows/columns specified in this file')
    args = parser.parse_args()


    # setup logging
    ch = logging.StreamHandler(stderr)
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(logging.Formatter('%(levelname)s\t%(asctime)s\t%(message)s'))
    LOG.addHandler(ch)

    ignoreRowNames = [set() for _ in args.matrices]
    ignoreColNames = [set() for _ in args.matrices]

    if args.ignore:
        if len(args.ignore) == 1:
            # one ignore file for all
            ignore_rows, ignore_cols = readOutlier(args.ignore[0])
            for i in xrange(len(args.matrices)):
                ignoreRowNames[i] = ignore_rows
                ignoreColNames[i] = ignore_cols

        elif len(args.ignore) == len(args.matrices):
            # one ignore file for each
            for i in xrange(len(args.matrices)):
                ignore_rows, ignore_cols = readOutlier(args.ignore[i])
                ignoreRowNames[i] = ignore_rows
                ignoreColNames[i] = ignore_cols
        else:
            LOG.fatal('Number of ignore files is not 1, nor matches the ' + 
                    'number of matrices. Exiting')
            exit(1)


    means = list()
    for i, f in enumerate(args.matrices):
        name = basename(f.name).rsplit('.', 1)[0]
        mtrx, colnames, rownames, _ = readHiCMapTRV(f, None, ignoreColNames[i],
                ignoreRowNames[i])
        assert mtrx.shape[0] == mtrx.shape[1]
        means.append((calculateMeans(mtrx), map(label2coords, colnames), name))

    srtd_means, start_pos, f = prepareData(means)
    gene_clusters = readBioGCs(args.gene_clusters)
    plotDistributions(srtd_means, start_pos, f, gene_clusters, stdout)
    
