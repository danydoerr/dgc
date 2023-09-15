#!/usr/bin/env python3

from sys import stdout, stderr, exit
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF, \
        FileType
from functools import partial
from itertools import chain, compress
import logging

import os
if not os.environ.get('DISPLAY', None):
    import matplotlib; matplotlib.use('Agg')

from scipy.cluster.hierarchy import linkage
from scipy.stats import norm as gaussian
import cooler as clr
import matplotlib.pylab as plt
import numpy as np

from identifyDiffGCs import fitGaussian, minPvalClusters, hier2graph, \
        computeMtrxSum, MX_VAL
from hic import condenseCoords, parseCoords, equiSegmentRegions

LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)


def getSamplePools(mtrx, sel_baited, sel_indices_row, sel_indices_col):
    """ this function assumes that the selected coordinates correspond to a
    single run in sel_indices_row and sel_indices_col, respectively """

    assert mtrx.shape[0] == mtrx.shape[1]

    pools_baited = list()
    pools_others = list()
    
    isSymmetric = np.all(sel_indices_row == sel_indices_col)
    start_diag, end_diag = getStartEndDiagonals(sel_indices_row,
            sel_indices_col, isSymmetric)

    n = sel_indices_row.size
    m = sel_indices_col.size

    # matrix of baited mask for debugging purposes
    #mask_mtrx = np.zeros((n, m))
    #mask_mtrx[:] = np.nan
    for k in xrange(start_diag, end_diag+1):
        x_mask = sel_baited[max(0, k):min(n+k, m)]
        y_mask = sel_baited[max(0, -k):n-max(0, k-m+n)]
        d = np.diag(mtrx, k=k)
        baited = d[np.logical_and(y_mask, x_mask)]
        pools_baited.append(baited[np.isfinite(baited)])
        others = d[np.logical_xor(y_mask, x_mask)]
        pools_others.append(others[np.isfinite(others)])
        #idx = np.arange(max(0, -k), n-max(0, k-m+n))
        #idy = np.arange(max(0, k), min(n+k, m))
        #d = np.zeros(x_mask.shape)
        #d[np.logical_and(y_mask, x_mask)] = 1
        #d[np.logical_xor(y_mask, x_mask)] = 2
        #mask_mtrx[idx, idy] = d

    return pools_baited, pools_others


def getStartEndDiagonals(sel_coord1, sel_coord2, upper_triangle_only=False):
    """ this function assumes that the selected coordinates correspond to a
    single run in sel_coord1 and sel_coord2, respectively """
    
    if upper_triangle_only:
        start_diag = 1
        end_diag = int(np.sum(sel_coord2)-1)
    else:
        m = int(np.sum(sel_coord1) + np.sum(sel_coord2))

        # i0\j0------j1
        #   |        |
        #   |        |
        #   |        |
        #   |        |
        #   i1-------o
        i1 = sel_coord1.size - np.argmax(np.flip(sel_coord1))
        j0 = np.argmax(sel_coord2)

        start_diag = j0-i1+1
        end_diag = start_diag + m-2

    return (start_diag, end_diag)


def getBaitedAndOthersCounts(sel_baited_row, sel_baited_col, isSymmetric=False):

    assert not isSymmetric or (sel_baited_row.size == sel_baited_col.size and
            np.all(sel_baited_row == sel_baited_col))


    # indexing scheme of diagonals
    #
    #   012---n
    #  -1012  |
    #   -1..2 |
    #   |  ..2|
    #   |   .12
    #   |   ..1
    #  -n-- -10

    n = sel_baited_row.size
    m = sel_baited_col.size

    if isSymmetric:
        start_diag = 1
    else:
        start_diag = 1-n

    end_diag = m-1

    baited_count = np.zeros(end_diag-start_diag+1, dtype=int)
    others_count = np.zeros(end_diag-start_diag+1, dtype=int)

    # matrix of baited mask for debugging purposes
    #mask_mtrx = np.zeros((n, m))
    for c, k in enumerate(xrange(start_diag, end_diag+1)):
        x_mask = sel_baited_col[max(0, k):min(n+k, m)]
        y_mask = sel_baited_row[max(0, -k):n-max(0, k-m+n)]
        #idx = np.arange(max(0, -k), n-max(0, k-m+n))
        #idy = np.arange(max(0, k), min(n+k, m))
        #d = np.zeros(x_mask.shape)
        #d[np.logical_and(y_mask, x_mask)] = 1
        #d[np.logical_xor(y_mask, x_mask)] = 2
        #mask_mtrx[idx, idy] = d

        baited_count[c] = np.sum(np.logical_and(y_mask, x_mask))
        others_count[c] = np.sum(np.logical_xor(y_mask, x_mask))
    return baited_count, others_count


def sampleSum(pools_baited, pools_others, baited_count, others_count, repeats):

    assert len(pools_baited) == len(pools_others) and baited_count.size == \
            others_count.size and len(pools_baited) == baited_count.size

    dist = np.empty(repeats)
    for k in xrange(repeats):
        instance = 0.
        for i in xrange(len(pools_others)):
            if pools_baited[i].size:
                instance += np.sum(np.random.choice(pools_baited[i],
                    baited_count[i]))
            if pools_others[i].size:
                instance += np.sum(np.random.choice(pools_others[i],
                    others_count[i]))
        dist[k] = instance

    #dist.sort()
    return dist


def doLinkage(mtrx, ltype):
    lmtrx = mtrx.copy()
    lmtrx[np.isnan(lmtrx)] = MX_VAL
    return linkage(lmtrx, ltype)



if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('hic_map', type=file, help='(Ca-)HI-C map in COOL format')
    parser.add_argument('coord1', type=str,
            help='coordinates of submatrix that is subject to p-value calculation')
    parser.add_argument('coord2', type=str, nargs='?',
            help='coordinates of submatrix that is subject to p-value calculation')
    parser.add_argument('-r', '--repeats', type=int, default=10000,
            help='Number of repeats used in sampling the ' + \
                    'empirical distribution')
    parser.add_argument('-p', '--plot_distribution', action='store',
            type=FileType('w'), help='plot distribution')
    parser.add_argument('-a', '--alpha', type=float, default=0.05,
             help='signficance threshold')
    args = parser.parse_args()

    # setup logging
    ch = logging.StreamHandler(stderr)
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(logging.Formatter('%(levelname)s\t%(asctime)s\t%(message)s'))
    LOG.addHandler(ch)


    LOG.info('loading (Ca-) Hi-C file %s' %args.hic_map.name)
    c = clr.Cooler(args.hic_map.name)
    mtrx = c.matrix()[:]
    bins = c.bins()[:]
    coords = zip(bins['chrom'].values, bins['start'].values, bins['end'].values)
    resolution = coords[0][2] - coords[0][1]

    if 'outlier' in bins.columns:
        sel = bins['outlier'].values  == False
        mtrx = mtrx[np.ix_(sel, sel)]
        coords = list(compress(coords, sel))
    else:
        sel = np.ones(bins.shape[0], dtype=bool)

    if 'baited' in bins.columns:
        LOG.info(('%s describes Capture Hi-C dataset - testing only ' + \
                'bait2bait / bait2other-end counts...') %args.hic_map.name)
        sel_baited = bins['baited'].values[sel]
    else:
        sel_baited = np.ones(sel.size, dtype=bool)


    LOG.info(('parsing coordinates and producing segmentation with ' + \
            'resolution %sbp') %(resolution))
    coord1 = map(lambda x: x[:2] + (x[2]+1, ), equiSegmentRegions( \
            (parseCoords(args.coord1), ), resolution))
    __coord1_set = set(coord1)
    sel_coord1 = np.array([x in __coord1_set for x in coords])

    if args.coord2:
        coord2 = map(lambda x: x[:2] + (x[2]+1, ), equiSegmentRegions( \
                (parseCoords(args.coord2), ), resolution))
        __coord2_set = set(coord2)
        sel_coord2 = np.array([x in __coord2_set for x in coords])
    else:
        coord2 = coord1
        sel_coord2 = sel_coord1

    isSymmetric = np.all(sel_coord1 == sel_coord2)
    submtrx = mtrx[np.ix_(sel_coord1, sel_coord2)]

    if isSymmetric:
        LOG.info('selected coordinates describe a symmetric submatrix, ' + \
                'therefore only using values in its upper triangle ' + \
                '(omitting the main diagonal), to revent bias in ' + \
                'computing p-values')
        m = submtrx.shape[0]
        linkageMtrx =  submtrx[np.triu_indices(m, k=1)]
        linkageLabels = list(compress(coords, sel_coord1))
    else:
        LOG.info('gene cluster matrix is assumed to be not symmetric, ' + \
                'computing p-values on the sum of *all* values')
        m = sum(submtrx.shape)
        X = np.empty((m, m))
        X[:] = np.nan
        X[:submtrx.shape[1],-submtrx.shape[0]:] = submtrx.transpose()
        linkageMtrx = X[np.triu_indices(m, k=1)]
        linkageLabels = list(chain(compress(coords, sel_coord1),
            compress(coords, sel_coord2)))

    if not linkageMtrx.size:
        LOG.fatal('linkage matrix is empty. exiting.')
        exit(1)

    baited_count, others_count = getBaitedAndOthersCounts( \
            sel_baited[sel_coord1], sel_baited[sel_coord2], isSymmetric = \
            isSymmetric)

    LOG.info('computing baited/ other-end pools for sampling')
    pools_baited, pools_others = getSamplePools(mtrx, sel_baited, sel_coord1,
        sel_coord2)
    LOG.info(('computing empirical distribution of sum of contact counts ' + \
            'based on sampling with %s repeats') %args.repeats)
    dist = sampleSum(pools_baited, pools_others, baited_count, others_count, \
            args.repeats)
    LOG.info('fitting a gaussian curve to sum of Hi-C counts of diagonals ' + \
            'from main matrix')
    mse, coeff, (dist_x, dist_y) = fitGaussian(dist, args.repeats/10,
            1)
    LOG.info('fitted gaussian has a mean squared error of %s' %mse)
    LOG.info('computing hierarchical clustering of submatrix')
    total = float(sum(baited_count) + sum(others_count))
    # calculate p-value for left-hand side (when smallest values are added
    # first)
    calculatePvalue = lambda s, val: gaussian.cdf(val, loc=s/total*coeff[0],
            scale=np.sqrt(s/total)*coeff[1])

    hier = doLinkage(linkageMtrx, 'single')
    G = hier2graph(hier)
    min_clusters = minPvalClusters(G, partial(computeMtrxSum, n=m,
            pDist=linkageMtrx), calculatePvalue, alpha=1, largestOnly=False)

    if any(map(lambda x: x[1] < 0, min_clusters)):
        min_clusters = filter(lambda x: x[1] < 0, min_clusters)

    min_cluster = min(min_clusters, key=lambda x: (x[0], x[2]))

    hier = doLinkage(linkageMtrx *-1, 'single')
    Gp= hier2graph(hier)

    # note: pDist=linkageMtrx is correct, because we are using the same
    # coefficients/distribution inferred from the diagonals of the Hi-C matrix,
    # except that we are now summing up from the big values to the low (the
    # right side of the distribution)
    clusters = minPvalClusters(Gp, partial(computeMtrxSum, n=m,
        pDist=linkageMtrx), lambda s, val: 1-calculatePvalue(s, val), alpha=1,
        largestOnly=False)
    max_clusters = list(clusters)
    if any(map(lambda x: x[1] > 0, max_clusters)):
        max_clusters = filter(lambda x: x[1] > 0, max_clusters)
    max_cluster = min(max_clusters, key=lambda x: (x[0], x[2]))

    root = [v for v in G.nodes() if G.in_degree(v) == 0][0]
    root_cluster = [c for c in clusters if c[-1] == root][0]

    if args.plot_distribution != None:
        # diagonal neighborhood of submatrix
        distPdf = lambda x: gaussian.pdf(x, loc=coeff[0], scale=coeff[1])
        plt.plot(dist_x, dist_y, label = 'sample pool' )
        plt.plot(dist_x, distPdf(dist_x), label = 'fitted gaussian')
        plt.title('distribution of Hi-C contact counts')
        plt.xlim([coeff[0]-3*coeff[1], coeff[0]+3*coeff[1]])
        plt.legend()
        plt.savefig(args.plot_distribution, format='pdf')

    labels = map(lambda x: x[:2], coords)
    l2i = dict(zip(labels, xrange(len(coords))))

    min_cluster_labels = sorted(set(map(lambda x: linkageLabels[x][:2],
        min_cluster[-2])))
    max_cluster_labels = sorted(set(map(lambda x: linkageLabels[x][:2],
        max_cluster[-2])))
    root_cluster_labels = sorted(set(map(lambda x: linkageLabels[x][:2],
        root_cluster[-2])))

    LOG.info('printing computed p-value')

    alpha = args.alpha/G.number_of_nodes()
    LOG.info('bonferroni-corrected significance threshold is %s' %alpha)
    print('high contact region\t%s\t%ssignificant\t%.5f\t%s\t%s\t%s' %(\
            max_cluster[0], max_cluster[0] > alpha and 'not ' or '',
            max_cluster[1], max_cluster[2], len(max_cluster[-2])/float(m),
            condenseCoords(max_cluster_labels, labels, l2i)))
    print('low contact region\t%s\t%ssignificant\t%.5f\t%s\t%s\t%s' %(\
            min_cluster[0], min_cluster[0] > alpha and 'not ' or '',
            min_cluster[1], min_cluster[2], len(min_cluster[-2])/float(m),
            condenseCoords(min_cluster_labels, labels, l2i)))
    print('full submatrix\t%s\t%ssignificant\t%.5f\t%s\t%s\t%s' %(\
            root_cluster[0], root_cluster[0] > args.alpha and 'not ' or '',
            root_cluster[1], root_cluster[2], len(root_cluster[-2])/float(m),
            condenseCoords(root_cluster_labels, labels, l2i)))

    LOG.info('DONE')
