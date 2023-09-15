#!/usr/bin/env python3

from sys import stdout, stderr, exit
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF, \
        FileType
from itertools import repeat, combinations, chain, product, compress
from os.path import basename
import logging
import csv

import cooler as clr
import numpy as np

LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)


if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('-c', '--raw_counts', action='store_true',
            help = 'output matrix with raw counts, rather than balanced counts')
    parser.add_argument('-s', '--chromosome_sizes', type=FileType('w'),
            help = 'write chromosome sizes to specified file')
    parser.add_argument('cool_file', type=file,
            help='Hi-C matrix in COOLER format')
    args = parser.parse_args()

    # setup logging
    ch = logging.StreamHandler(stderr)
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(logging.Formatter('%(levelname)s\t%(asctime)s\t%(message)s'))
    LOG.addHandler(ch)

    LOG.info('loading COOLER %s' %args.cool_file.name)
    c = clr.Cooler(args.cool_file.name)

    if args.chromosome_sizes:
        LOG.info('writing chromosome sizes')
        for name, size in zip(c.chromnames, c.chromsizes):
            print('\t'.join((name, str(size))) args.chromosome_sizes.close(), file=args.chromosome_sizes)

    out = stdout
    LOG.info('reading %s matrix' %(args.raw_counts and 'raw count' or \
            'balanced'))
    mtrx = c.matrix(balance = not args.raw_counts)[:]
    bins = c.bins()[:]
    coords = zip(bins['chrom'].values, map(str, bins['start'].values))

    LOG.info('writing output in Juicer\'s short format')
    # <str1> <chr1> <pos1> <frag1> <str2> <chr2> <pos2> <frag2> <score>
    for (chrx, chry) in chain(zip(*repeat(c.chromnames, 2)),
            combinations(c.chromnames, 2)):
        sel_chrx = bins['chrom'].values == chrx
        sel_chry = bins['chrom'].values == chry
        coord1 = list(compress(coords, sel_chrx))
        coord2 = list(compress(coords, sel_chry))

        submtrx = mtrx[np.ix_(sel_chrx, sel_chry)]
        if chrx == chry:
            pairs = combinations(xrange(submtrx.shape[0]), 2)
        else:
            pairs = product(xrange(submtrx.shape[0]), xrange(submtrx.shape[1]))

        for i, j in pairs:
            if np.isfinite(submtrx[i,j]):
                print('\t'.join((
                    '0', coord1[i][0], coord1[i][1], '0',
                    '0', coord2[j][0], coord2[j][1], '1',
                    str(submtrx[i,j]))), file=out)

    LOG.info('DONE')
