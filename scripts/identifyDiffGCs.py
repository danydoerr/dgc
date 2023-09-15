#!/usr/bin/env python3

from sys import stdout, stderr, exit, maxint
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF, \
        FileType
from itertools import chain, combinations, product
from functools import partial
from random import sample
import logging
import csv

import os
if not os.environ.get('DISPLAY', None):
    import matplotlib; matplotlib.use('Agg')

from matplotlib import pylab as plt
from scipy.cluster.hierarchy import linkage 
from scipy.misc import comb
from scipy.optimize import least_squares
from scipy.stats import norm as gaussian 
import networkx as nx
import numpy as np
import pandas as pd

from identifyOutlier import readOutlier
from hic import readHiCMapTRVAsPdist, pDistIndex, condenseCoords

MX_EMPIRICAL_DIST = 100
GAUSSIAN_HIST_RESOLUTION = 1000
LEAST_SQ_INLIER_MARGIN = 1
SUBSAMPLE_SIZE = 100000

GDATA_LEAVES = 'leaves'
GDATA_MTRXSUM = 'submatrix_sum'
GDATA_DELAYOUT = 'delay_output'

MX_VAL = maxint 

LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)


def __gaussian__(p, x, y=0):
    return gaussian.pdf(x, loc=p[0], scale=p[1]) - y


def fitGaussian(data, hist_size, inlier_margin):
    hist, bins = plt.histogram(data, hist_size, density=True)
    bin_centres = (bins[:-1] + bins[1:])/2
    p0 = [np.mean(data), np.sqrt(np.var(data))]
    try: 
        res = least_squares(__gaussian__, p0, loss='soft_l1',
                f_scale=inlier_margin, args=(bin_centres, hist))
    except ValueError, e:
        LOG.error('Fitting of Gaussian curve failed: %s' %e.message)
        return float('inf'), (np.mean(data), np.var(data)), (bin_centres, hist)
    # mean squared error as quality of fit
    mse = np.sum(np.square(res.fun-hist))/len(bin_centres)
    return mse, res.x, (bin_centres, hist)


def hier2graph(hier):

    G = nx.DiGraph()

    n = len(hier)+1
    G.add_nodes_from(xrange(n), dist=-np.inf)
    for i, line in enumerate(hier):
        G.add_node(i+n, dist=line[2])
        u = int(line[0])
        v = int(line[1])

#        if G.node[u]['dist'] == line[2] or G.node[v]['dist'] == line[2]:
#            G.add_edges_from((i+n, w) for w in chain(G.neighbors(u),
#                G.neighbors(v)))
#            G.remove_node(u)
#            G.remove_node(v)
#        else:
        G.add_edge(i+n, u)
        G.add_edge(i+n, v)
    return G


def depth(G):

    res = dict()
    queue = [(v, 0) for v in G.nodes() if G.in_degree(v) == 0]
    while queue:
        v, d = queue.pop()
        res[v] = d
        for u in G.neighbors(v):
            queue.append((u, d+1))
    return res


def minPvalClusters(G, subgraphValue, calculatePvalue, alpha, largestOnly=False):
    """
        Traverses hierarchy from leaves to root, thereby computing the pvalue of
        the 'submatrix value' from each induced subtree. 
        
        If largestOnly is set to False, all subgraphs are reported whose pvalues
        are below alpha. Otherwise, if the pvalue of a node is below alpha, it
        will output those subtrees whose pvalues are blow, UNLESS the node is
        included in a bigger subtree whose pvalue is below alpha. Then the
        bigger subtree will be returned (unless ... recursively).
        
        In the output, subtrees are presented by their leaves. 

        The function uses two functions as input:

            subgraphValue(G, v) -- G: Graph, v: node in G that induces subgraph

            calculatePvalue(n, val) -- n: number of leaves of the subgraph, val
            -- subgraph value of subgraph
    """

    res = list()
    v2d = depth(G)
    leaves = [v for v in G.nodes() if G.out_degree(v) == 0]

    mxd = max(map(v2d.get, leaves))
    queue = [set() for _ in xrange(mxd+1)]

    for v in leaves:
        queue[v2d[v]].add(v)

    for d in xrange(mxd, -1, -1):
        for v in queue[d]:
            # p denotes the predecessor of v; it is None if we are at the root
            p = None
            if d > 0:
                # by definition, there is only one predecessor
                p = next(G.predecessors(v))
                if largestOnly and not G.node[p].has_key(GDATA_DELAYOUT):
                    G.node[p][GDATA_DELAYOUT] = list()
                if not G.node[p].has_key(GDATA_LEAVES):
                    G.node[p][GDATA_LEAVES] = list()
                G.node[p][GDATA_LEAVES].extend(G.node[v].get(GDATA_LEAVES,
                    (v,)))
                # enqueue predecessors for the next level
                queue[v2d[p]].add(p)

            if G.out_degree(v):
                val, count = subgraphValue(G, v)
                pval = 1
                if count:
                    pval = calculatePvalue(count, val)
                if pval <= alpha:
                    if not largestOnly or d == 0:
                        res.append((pval, val, count, G.node[v][GDATA_LEAVES], v))
                    elif largestOnly:
                        G.node[p][GDATA_DELAYOUT].append((pval, val, count, v))
                elif largestOnly:
                    if d == 0:
                        res.extend(map(lambda x: x[:-1] +
                            (G.node[x[-1]][GDATA_LEAVES], x[-1]),
                            G.node[v].get(GDATA_DELAYOUT, ())))
                    else:
                        G.node[p][GDATA_DELAYOUT].extend(G.node[v].get(
                            GDATA_DELAYOUT, ()))
    return res


def computeMtrxSum(G, v, n, pDist):
    count = 0
    val = 0
    for u in G.successors(v):
        if G.node[u].has_key(GDATA_MTRXSUM):
            val += G.node[u][GDATA_MTRXSUM][0]
            count += G.node[u][GDATA_MTRXSUM][1]

    for u, w in combinations(G.successors(v), 2):
        for (c1, c2) in product(G.node[u].get(GDATA_LEAVES, (u, )),
                G.node[w].get(GDATA_LEAVES, (w, ))):
            # sorting necessary for calling pDistIndex
            if c1 > c2:
                c1, c2 = c2, c1
            k = pDistIndex(n, c1, c2)
            if pDist[k] < MX_VAL:
                val += pDist[k]
                count += 1
    # store in tree for next calculation
    G.node[v][GDATA_MTRXSUM] = (val, count)
    return val, count


def getCDF(G, dists, data, s):

    assert s >= 1

    cdf = None

    #
    # dists is an array of estimated gaussian distributions from sum counts of
    # submatrices, so that dists[s] is the distribution nof the sum of s+2 terms
    # 
    # if the corresponding distribution has not been estimated (dists[s] ==
    # None), it will be done so here
    # 
    # if s > len(dists), then the sum distribution will be extrapolated based on
    # the coefficients of the last distribution dists[-1]. The folowing
    # estimation is an exact estimation of the the distribution of sums of iid
    # normal random variables. It applies to the empirical data (which might not
    # be normal distributed) because of the central limit theorem.  

    if s-1 < len(dists):
        if dists[s-1] == None:
            empirical = [sum(sample(data, s)) for _ in xrange(SUBSAMPLE_SIZE)]
            mse, coeff, (pool_x, pool_y) = fitGaussian(empirical,
                    GAUSSIAN_HIST_RESOLUTION, LEAST_SQ_INLIER_MARGIN)
            cdf = lambda x: gaussian.cdf(x, loc=coeff[0], scale=coeff[1])
            dists[s-1] = [cdf, coeff[0], coeff[1]]
        else:
            cdf = dists[s-1][0]
    else:
        if dists[-1] == None:
            empirical = [sum(sample(data, len(dists)+1)) for _ in
                    xrange(SUBSAMPLE_SIZE)]
            mse, coeff, (pool_x, pool_y) = fitGaussian(empirical,
                    GAUSSIAN_HIST_RESOLUTION, LEAST_SQ_INLIER_MARGIN)
            cdf = lambda x: gaussian.cdf(x, loc=coeff[0], scale=coeff[1])
            dists[-1] = [cdf, coeff[0], coeff[1]]
        else:
            coeff = dists[-1][1:3]

        k = float(s)/(len(dists)+1)
        cdf = lambda x: gaussian.cdf(x, loc=k*coeff[0],
                scale=np.sqrt(k)*coeff[1])
            
    return cdf


def printClusters(clusters, labels, optimization, out):

    regions = [(a, b, map(lambda x: labels[x], c)) for a, b, _, c, _ in
            clusters]
    regions.sort(key=lambda x: x[-1][0])
   
    l2i = dict(zip(labels, xrange(len(labels))))

    for pval, val, r in regions:
        r.sort()
        coords = condenseCoords(r, labels, l2i)
        print('%s\t%s\t%s' %(pval, optimization == 'min' and val or -val, coords), file=out)


if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('diff_map', type=file, help='diff\'ed (Ca-)HI-C map')
    parser.add_argument('-a', '--alpha', type=float, default=0.05,
            help='P-value significance threshold')
    parser.add_argument('-b', '--bonferroni', action='store_true',
            help='Bonferroni multiple testing correction lowers alpha ' + \
                    '(significance threshold)')
    parser.add_argument('-m', '--max_dist_from_diag', type=int, default=1000000, 
            help='maximal distance in base pairs that is considered in ' + \
                'spatial folding of differential gene cluster analysis')
    parser.add_argument('-c', '--clustering_method', choices=('single', \
            'complete', 'average'), default='average', \
            help='type of hierarchical clustering')
    parser.add_argument('-i', '--ignore', type=str, 
            help='ignore rows/columns specified in this file')
    parser.add_argument('-p', '--plot_distribution', action='store',
            type=FileType('w'), help='plot distribution')
    parser.add_argument('-o', '--optimize', choices=('max', 'min'),
            default='max')
    args = parser.parse_args()

    # setup logging
    ch = logging.StreamHandler(stderr)
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(logging.Formatter('%(levelname)s\t%(asctime)s\t%(message)s'))
    LOG.addHandler(ch)

    ignoreRowNames = list()
    ignoreColNames = list()

    if args.ignore:
        LOG.info('reading %s containing row / columns that will be ignored' %(
            args.ignore))
        ignoreRowNames, ignoreColNames  = readOutlier(open(args.ignore))

    pDist, labels, _, _ = readHiCMapTRVAsPdist(args.diff_map,
            args.max_dist_from_diag, ignoreColNames, ignoreRowNames)
    if args.optimize == 'max':
        pDist[pDist != None] *= -1
    pDist[pDist == None] = MX_VAL
    # change dtype of array to float
    pDist = np.array(pDist, dtype=float)
    hier = linkage(pDist, args.clustering_method)
    G = hier2graph(hier)

    LOG.info('fitting a gaussian curve to Hi-C count data')
    data = pDist[pDist < MX_VAL]

    gaussian_dists = [None] * MX_EMPIRICAL_DIST
    mse, coeff, (test_x, test_y) = fitGaussian(data, GAUSSIAN_HIST_RESOLUTION,
            LEAST_SQ_INLIER_MARGIN)
    cdf = lambda x: gaussian.cdf(x, loc=coeff[0], scale=coeff[1])
    gaussian_dists[0] = [cdf, coeff[0], coeff[1]]

    LOG.info('fitted gaussian has a mean squared error of %s' %mse)
    if args.plot_distribution != None:
        LOG.info(('plotting empirical and fitted distribution to ' + \
                '%s') %args.plot_distribution.name)
        plt.plot(test_x, test_y, label = 'observed distribution')
        plt.plot(test_x, __gaussian__(coeff, test_x), 
                label = 'fitted normal distribution')
        plt.title('differential contact count distribution')
        plt.xlim([coeff[0]-5*coeff[1], coeff[0]+5*coeff[1]])
        plt.legend()
        #plt.show()
        plt.savefig(args.plot_distribution, format='pdf')

    alpha = args.alpha 
    if args.bonferroni:
        LOG.warning('correction for multiple testing is activated, ' + \
                'setting significance threshold from %s to %s' %(alpha, \
                alpha/float(G.number_of_nodes())))
        alpha /= float(G.number_of_nodes())

    cdfPartial = partial(getCDF, G=G, dists=gaussian_dists, data=data)
    clusters = minPvalClusters(G, partial(computeMtrxSum, n=len(labels), pDist=pDist),
            lambda s, val: cdfPartial(s=s)(val), alpha, largestOnly=True)
    printClusters(clusters, labels, args.optimize, stdout)

    LOG.info('DONE')
