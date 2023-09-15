#!/usr/bin/env python

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF
from sys import stdout, stderr, exit
import logging
import math
import csv

LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)

from annotation2circos import readGFF3Data
from verifyWithDE import readDEData


if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('annotation', type=file,
            help='Annotation file in GFF format')
    parser.add_argument('expression_data', type=file,
            help='Expression data in tabular format')
    parser.add_argument('condition', type=str,
            help='column in expression table')
    parser.add_argument('-l', '--log', action='store_true',
            help='column in expression table')


    args = parser.parse_args()

    #
    # setup logging
    #
    ch = logging.StreamHandler(stderr)
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(logging.Formatter('%(levelname)s\t%(asctime)s\t%(message)s'))
    LOG.addHandler(ch)

    out = stdout

    LOG.info('reading gene annotations from %s' %args.annotation.name)
    genes = readGFF3Data(args.annotation, ('gene', ), ('protein_coding', ))

    LOG.info('reading expression data from %s' %args.expression_data.name)
    expression = readDEData(args.expression_data, (args.condition, ))

    gene2coords = dict(map(lambda x: (x[-1].upper(), x[:-1]), genes))
    LOG.info('writing expression data to stdout')
    for gene, expr in expression:
        locus = gene.upper()
        if gene2coords.has_key(locus):
            if args.log:
                if type(expr) is float and expr >= 1:
                    expr = math.log(expr)
                else:
                    expr = ''
            if type(expr) is float:
                chrx, start, end = gene2coords[locus][:3]
                print('%s %s %s %.5f locus=%s' %(chrx, start, end, expr, locus), file=out)
        else:
            LOG.warning('Unable to find gene %s in GFF file' %locus)
    LOG.info('DONE')
