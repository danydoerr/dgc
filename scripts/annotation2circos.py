#!/usr/bin/env python3

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF
from sys import stdout, stderr, exit
import logging
import csv
import re

from plotSubmatrix import PAT_CHR, PAT_NAME, PAT_LOCUS, PAT_BIOTYPE

LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)

PAT_SYNONYMS = re.compile('^(?:.*;)?gene_synonym=([^;]+)(?:;.*)?$', re.I)
PAT_GENE = re.compile('^(?:.*;)?gene=([^;]+)(?:;.*)?$', re.I)

def readGFF3Data(data, feature_types, bio_types=(),
        ignoreHypotheticalAndPseudo=True):
    ''' parser for GFF3 annotation file '''

    chrMap = dict()
    anns = set()
    feature_types = set(feature_types)
    bio_types = set(map(lambda x: x.lower(), bio_types))
    for line in csv.reader(data, delimiter='\t'):
        if line[0].startswith('#'):
            continue
        if line[2] == 'region':
            m = PAT_CHR.match(line[8])
            if m:
                chrMap[line[0]] = m.group(1)
            else:
                m = PAT_NAME.match(line[8])
            if m and not chrMap.has_key(line[0]):
                chrMap[line[0]] = m.group(1)
        elif line[2] in feature_types and (not ignoreHypotheticalAndPseudo or \
                (line[8].find('pseudogene') < 0 and \
                line[8].find('hypothetical') < 0)):
            m = PAT_BIOTYPE.match(line[8])
            if not bio_types or (m and m.group(1).lower() in bio_types):
                locus = ''
                m = PAT_LOCUS.match(line[8])
                if m:
                    locus = m.group(1)

                identifiers = set([locus, line[0]])
                m = PAT_NAME.match(line[8])
                if m:
                    identifiers.add(m.group(1))
                m = PAT_GENE.match(line[8])
                if m:
                    identifiers.add(m.group(1).replace('%3B', '\t'))
#                m = PAT_SYNONYMS.match(line[8])
#                if m:
#                    identifiers.update(m.group(1).split(','))
                anns.add((line[0], int(line[3]), int(line[4]), line[6], line[2],
                    locus, tuple(identifiers)))
    res = list()
    for ann in anns:
        res.append((chrMap.get(ann[0], ann[0]), ) + ann[1:])

    # inefficient, but hey, data is not that large
    res.sort(key=lambda x: (x[0], x[1], -x[2]))
    return res

if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('annotation', type=file,
            help='Annotation file in GFF format')

    args = parser.parse_args()

    #
    # setup logging
    #
    ch = logging.StreamHandler(stderr)
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(logging.Formatter('%(levelname)s\t%(asctime)s\t%(message)s'))
    LOG.addHandler(ch)

    LOG.info('reading gene annotations from %s' %args.annotation.name)
    genes = readGFF3Data(args.annotation, ('gene', ), ('protein_coding', ))

    LOG.info('writing identified protein-coding genes to stdout')
    out = stdout
    for gene in genes:
        print('%s %s %s locus=%s' %(gene[0], gene[1], gene[2], gene[-2].upper(), file=out)
    LOG.info('DONE')
