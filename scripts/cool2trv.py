#!/usr/bin/env python3

from sys import stdout, stderr, exit
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF, \
        FileType
from os.path import basename
import logging
import csv

import cooler as clr
import numpy as np

from hic import writeMtrx as writeMtrxTRV

LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)


if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('-c', '--raw_counts', action='store_true',
            help='output matrix with raw counts, rather than balanced counts')
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

    LOG.info('reading bins')
    names = ['%s-%s' %(x.chrom.values[0], int(x.start)) for x in c.bins()]
    
    LOG.info('reading %s matrix' %(args.raw_counts and 'raw count' or \
            'balanced'))

    orig = c.matrix(balance=not args.raw_counts)[:]
    mtrx = np.array(orig, dtype=object)
    mtrx[np.logical_or(orig == 0, np.isnan(orig))] = None 

    LOG.info('writing output in TRV format')
    writeMtrxTRV(mtrx, names, names, stdout, xlabel='%s:cool2trv' %basename( \
            args.cool_file.name).rsplit('.', 1)[0])
    LOG.info('DONE')
