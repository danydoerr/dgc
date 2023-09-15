#!/usr/bin/env python3

from sys import stdout, stderr, exit, maxint
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF
from bisect import bisect_left
import csv
import re


PAT_BIOGC = re.compile('^(\w+):(\d+)-(\d+)$')


def readBioGCs(data):

    res = list()
    for name, _, coords in csv.reader(data, delimiter='\t'):
        m = PAT_BIOGC.match(coords)
        if not m:
            print(('unable to parse coords %s of GC %s - ' + \
                    'expected string matching %s') %(name, coords, \
                    PAT_BIOGC.pattern), file=stderr)
            exit(1)

        chrx, start, end = m.groups()
        res.append((name, (chrx, int(start), int(end))))

    return res


def parseCandidates(data):

    res = list()
    for i, (pval, val, coords) in enumerate(csv.reader(data, delimiter='\t')):
        candidate = (float(pval), float(val), list())
        for same_chr in coords.split(';'):
            chrx, regions = same_chr.split(':',1)
            candidate[-1].append((chrx, list())) 
            for block in regions.split(','):
                candidate[-1][-1][1].append(tuple(map(lambda x: x != 'end' and
                    int(x) or maxint, block.split('-'))))
            candidate[-1][-1][1].sort()
        res.append(candidate)
    return res


def checkOverlap(realClusters, candidates):

    res = [set() for _ in xrange(len(realClusters))]

    chr2cl = dict()
    for i, (_, (chrx, start, end)) in enumerate(realClusters):
        if not chr2cl.has_key(chrx):
            chr2cl[chrx] = list()
        chr2cl[chrx].append((start, end, i))

    for occ in chr2cl.values():
        occ.sort()

    for j, c in enumerate(candidates):
        for chrx, regions in c[-1]:
            gcs = chr2cl.get(chrx, list())

            i = 0
            for start, end, p in gcs:
                i = bisect_left(regions, (start, start), i)
                if i < len(regions) and regions[i][0] < end:
                    res[p].add(j)
    return res


def printCluster(cluster, out):

    pval, val, r = cluster
    out.write('%s\t%s\t' %(pval, val))
    for j, (chrx, regions) in enumerate(r):
        if j:
            out.write(';')
        out.write('%s:' %chrx)
        for i, (start, end) in enumerate(regions):
            if i:
                out.write(',')
            out.write('%s-%s' %(start, end))


def printHits(realClusters, candidates, overlaps, out):

    for i in xrange(len(realClusters)):
        if overlaps[i]:
            print('%s\t%s:%s-%s:' %((realClusters[i][0], ) + \
                    realClusters[i][1]), file=out)
            for j in overlaps[i]:
                out.write('\t')
                printCluster(candidates[j], out)
                out.write('\n')


if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('biological_geneclusters', type=file, 
            help='tab-separated list of name, genes, and coordinates of GCs')
    parser.add_argument('gc_candidates', type=file, 
            help='tab-separated list of p-value, GC candidates ')
    args = parser.parse_args()

    realClusters = readBioGCs(args.biological_geneclusters) 
    candidates = parseCandidates(args.gc_candidates)
    overlaps = checkOverlap(realClusters, candidates)
    printHits(realClusters, candidates, overlaps, stdout)

