#!/usr/bin/env python

from sys import stdout, stderr, exit
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF, \
        FileType

import numpy as np
from scipy.stats import norm

from hic import readHiCMapTRV, label2coords, readRMap, readBaitmapIDs

DEFAULT_THRESHOLD = 0.1


def solve(m1, m2, std1, std2):
    a = 1./(2.*std1**2) - 1./(2.*std2**2)
    b = m2/(std2**2) - m1/(std1**2)
    c = m1**2 /(2*std1**2) - m2**2 / (2*std2**2) - np.log(std2/std1)
    return np.roots([a,b,c])

def identifyOutlier(mtrx, threshold):

    rows = list()
    cols = list()

    _mtrx = mtrx[mtrx != None]
    if not _mtrx.size:
        return rows, cols

    x = _mtrx.mean()
    v = _mtrx.var()
    for i in xrange(mtrx.shape[0]): 
        row = mtrx[i, :]
        _row = row[row != None]
        # compute mean/variance only in case row is non-empty (i.e., row is not
        # composed of all 'None's)
        if _row.size:
            rx = _row.mean()
            rv = _row.var()
            area = 0
            if rv > 0:
                roots = solve(x, rx, v, rv)
                if x < rx:
                    area = norm.cdf(roots[0],rx,rv) + (1.-norm.cdf(roots[0], x, v))
                else:
                    area = norm.cdf(roots[0],x,v) + (1.-norm.cdf(roots[0], rx, rv))
            if area < threshold:
                 rows.append(i)

    for j in xrange(mtrx.shape[1]): 
        col = mtrx[:, j]
        _col = col[col != None]
        if _col.size:
            cx = _col.mean()
            cv = _col.var()
            area = 0
            if cv > 0:
                roots = solve(x, cx, v, cv)
                if x < rx:
                    area = norm.cdf(roots[0],cx,cv) + (1.-norm.cdf(roots[0], x, v))
                else:
                    area = norm.cdf(roots[0],x,v) + (1.-norm.cdf(roots[0], cx, cv))
            if area < threshold:
                cols.append(j)
    return rows, cols


def readOutlier(data):

    cols = list()
    rows = list()
    for line in data:
        if line.startswith('rows\t'):
            rows = line.strip().split('\t')[1:]
        elif line.startswith('columns\t'):
            cols = line.strip().split('\t')[1:]

    return map(label2coords, rows), map(label2coords, cols)

if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('hic_map', type=str, help='(Ca-)HI-C map')
    parser.add_argument('-t', '--threshold', type=float,
            default=DEFAULT_THRESHOLD,
            help='Outlier threshold: minimum overlap between outlier ' + \
                    'distribution and overall distribution')
    parser.add_argument('-d', '--restriction_digest', type=str, 
            help='restriction fragments file in RMAP format, required ' + \
                    'for CaHi-C maps')
    parser.add_argument('-b', '--baitmap', type=str, 
            help='baitmap from the capture part of CHI-C, required for ' + \
                    'CaHi-C maps')
    args = parser.parse_args()
    
    probedFragments = None
    if args.restriction_digest and args.baitmap:
        rMap = readRMap(open(args.restriction_digest))
        rMapDict = dict((frag[3], frag[:3]) for frag in rMap)
        probeIDs = readBaitmapIDs(open(args.baitmap))
        probedFragments = set('%s-%s'%rMapDict[x][:2] for x in probeIDs)


    mtrx, colnames, rownames, _ = readHiCMapTRV(open(args.hic_map),
            regions=probedFragments)
    ignore_rows, ignore_cols = identifyOutlier(mtrx, args.threshold)

    print('rows\t%s' %('\t'.join(map(lambda x: rownames[x], ignore_cols))))
    print('columns\t%s' %('\t'.join(map(lambda x: colnames[x], ignore_cols))))
