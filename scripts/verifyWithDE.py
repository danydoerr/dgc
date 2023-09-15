#!/usr/bin/env python

from sys import stdout, stderr, exit
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF, \
        FileType
from bisect import bisect_left
from itertools import izip
import logging
import csv
import re

from testGCinCandidates import parseCandidates, printCluster

LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)

PAT_CHR = re.compile('^(?:.*;)?chromosome=([^;]+)(?:;.*)?$', re.I)
PAT_NAME = re.compile('^(?:.*;)?Name=([^;]+)(?:;.*)?$', re.I)
PAT_TAIR = re.compile('^(?:.*[,;])?TAIR:([^,;]+)', re.I)


def readGFF3Data(data, feature_types):
    ''' parser for GFF3 annotation file '''

    chrMap = dict()
    anns = list()
    feature_types = set(feature_types)
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
        elif line[2] in feature_types:
            name = line[0]
            m = PAT_NAME.match(line[8])
            if m:
                name = m.group(1)
            m = PAT_TAIR.match(line[8])
            if m:
                tairID = m.group(1)
                anns.append((line[0], int(line[3]), int(line[4]), tairID,
                    line[6], line[2], name))
    res = map(lambda x: (chrMap.get(x[0], x[0]), ) + x[1:], anns)
    # inefficient, but hey, data is not that large
    res.sort(key=lambda x: (x[0], x[1], -x[2]))

    return res 


def readDEData(data, tissue_types):
   
    res = list()

    tissue_idx = None
    isFirst = True
    for line in csv.reader(data, delimiter='\t'):
        if line[0].startswith('#'):
            continue
        if isFirst:
            isFirst = False
            header = map(lambda x: x.lower(), line)
            tissue_idx = list()
            for t in tissue_types:
                try:
                    tissue_idx.append(header.index(t.lower()))
                except ValueError, e:
                    LOG.fatal('tissue %s not found in DE file %s' %(t,
                        data.name))
                    exit(1)
            continue
        res.append((line[0], ) + tuple(line[i] and float(line[i]) for i in
            tissue_idx))

    return res


def identifyDEGenesInCandidates(candidates, annotation, deMap):

    res = list()

    for candidate in candidates:
        c = 0
        hits = list()
        for chrx, regions in candidate[-1]:
            i = 0
            for start, end in regions:
                i = bisect_left(annotation, (chrx, start, start) + (None,) *
                        (len(annotation[0])-3), i)
                if i > 0 and annotation[i-1] == chrx and annotation[i-1][2] < end:
                    i -= 1
                while i < len(annotation) and annotation[i][:2] < (chrx, end):
                    tairID = annotation[i][3]
                    if deMap.has_key(tairID) and '' not in deMap[tairID]:
                        c += 1
                        if deMap[tairID][0] < deMap[tairID][1]:
                            hits.append(i)
                    i += 1
        res.append((c, hits))
    return res


def printResults(candidates, annotation, deMap, hits, tissues, coverage, out):

    for candidate, (c, hit) in izip(candidates, hits):
        cov = float(len(hit))/(c or 1) * 100
        if hit and cov >= coverage:
            printCluster(candidate, out)
            print >> out, ('\n\t%.2f%% (%s/%s) genes in these regions are ' + \
                    'underexpressed in %s:\n') %(cov, len(hit), c, tissues[0])
            for pos in hit:
                tairID = annotation[pos][3]
                name = annotation[pos][-1]
                print >> out, '\t%.3f (%s/%s)\t%s\t%s' %(deMap[tairID][0] /
                        deMap[tairID][1], deMap[tairID][0], deMap[tairID][1],
                        tairID, name)
            out.write('\n')

if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('gc_candidates', type=file, 
            help='tab-separated list of p-value, GC candidates ')
    parser.add_argument('differential_expression_profile', type=file, 
            help='tab-separated list of differentially expressed genes ' + \
                    'in different tissues')
    parser.add_argument('annotation', type=file, 
            help='annotation file in GCF format')
    parser.add_argument('-c', '--coverage', type=float, default = 0,
            help='percentage (value between 0 and 100) of genes that must ' + \
                    'be underexpressed in region in order to be reported')
    parser.add_argument('-t', '--tissues', nargs = 2, default = ('leaf', 'root'),
            help='tissues to compare')
    args = parser.parse_args()

    # setup logging
    ch = logging.StreamHandler(stderr)
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(logging.Formatter('%(levelname)s\t%(asctime)s\t%(message)s'))
    LOG.addHandler(ch)

    LOG.info('reading GC candidate list %s' %args.gc_candidates.name)
    candidates = parseCandidates(args.gc_candidates)
    LOG.info('reading genome annotation in GFF format %s' %args.annotation.name)
    anns = readGFF3Data(args.annotation, ['gene'])
    LOG.info('reading differential expression profile ' + \
            args.differential_expression_profile.name)
    deData = readDEData(args.differential_expression_profile, args.tissues)
    deMap = dict((x[0], x[1:]) for x in deData)

    LOG.info(('identifying differentially expressed genes within GC ' + \
            'candidate regions %s') %args.gc_candidates.name)
    hits = identifyDEGenesInCandidates(candidates, anns, deMap) 
    LOG.info(('identified %s out of %s candidates with at least %s%% DE ' + \
            'genes. printing') %(len(filter(lambda x: x[1] and \
            float(len(x[1]))/(x[0] or 1) * 100 >= args.coverage, hits)), \
            len(hits), args.coverage))
    printResults(candidates, anns, deMap, hits, args.tissues, args.coverage,
            stdout)

    LOG.info('DONE')
