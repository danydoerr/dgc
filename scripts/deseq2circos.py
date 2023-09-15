#!/usr/bin/env python

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF
from sys import stdout, stderr, exit
import logging
import math
import csv
import re

LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)

from annotation2circos import readGFF3Data
from verifyWithDE import readDEData


def readDeseq(data):

    res = list()
    isFirst = True
    for line in csv.reader(data, delimiter='\t'):
        if isFirst:
            isFirst = False
            continue
        res.append((line[0], float(line[2]), float(line[5])))
    return res


if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('annotation', type=file,
            help='Annotation file in GFF format')
    parser.add_argument('expression_data', type=file,
            help='Expression data in tabular format')

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
    genes = readGFF3Data(args.annotation, ('gene', 'lncRNA'),
            ignoreHypotheticalAndPseudo=False)

    LOG.info('reading expression data from %s' %args.expression_data.name)
    expression = readDeseq(args.expression_data)

    gene2coords = dict()
    for gene in genes:
        for ident in gene[-1]:
            ident = ident.upper()
            if gene2coords.has_key(ident):
                LOG.warning('gene identifier %s not unique' %ident)
            else:
                gene2coords[ident] = gene

    matched_locus = set()
    LOG.info('writing expression data to stdout')
    for gene, log2FoldChange, pvalue in expression:
        locus = gene.strip('"').upper()
        if gene2coords.has_key(locus):
            if locus in matched_locus:
                LOG.warning('gene %s matched more than once' %locus)
		continue
            chrx, start, end = gene2coords[locus][:3]
            print('%s %s %s %.5f locus=%s,pvalue=%s' %(chrx, start,
                    end, log2FoldChange, gene2coords[locus][-2], pvalue), file=out)
            matched_locus.add(locus)
        else:
            LOG.warning('Unable to find gene %s in GFF file' %locus)
    LOG.info('DONE')
