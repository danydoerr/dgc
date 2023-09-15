#!/usr/bin/env python3

from sys import stdout, stderr, exit
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF, \
        FileType
from functools import partial
from itertools import chain, compress, product, izip, combinations
import logging

import os
if not os.environ.get('DISPLAY', None):
    import matplotlib; matplotlib.use('Agg')

from matplotlib import pylab as plt
import cooler as clr
import numpy as np

from hic import condenseCoords, parseCoords, equiSegmentRegions
from calculatePvalue import getSamplePools, getStartEndDiagonals

LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)


def sampleContactAreas(baited_mask, others_mask, pools_baited, pools_others,
        repeats, upper_triangle_only=False, doPositive=True,
        return_sample_representative=False):

    assert len(pools_baited) == len(pools_others) and baited_mask.shape == \
            others_mask.shape and not np.any(baited_mask & others_mask)

    n, m = baited_mask.shape
    start_diag, end_diag = getStartEndDiagonals(np.ones(n), np.ones(m),
            upper_triangle_only)

    representative = None
    dist = np.empty(repeats)
    for r in xrange(repeats):
        instance = np.zeros((n, m), dtype=float)
        instance.fill(np.nan)

        for i, k in enumerate(xrange(start_diag, end_diag+1)):
            baited_diag = np.diag(baited_mask, k)
            others_diag = np.diag(others_mask, k)
            idx = np.arange(max(0, -k), n-max(0, k-m+n))
            idy = np.arange(max(0, k), min(n+k, m))
            if pools_baited[i].size:
                instance[idx[baited_diag], idy[baited_diag]] = np.random.choice(
                        pools_baited[i], np.sum(baited_diag))
            if pools_others[i].size:
                instance[idx[others_diag], idy[others_diag]] = np.random.choice(
                        pools_others[i], np.sum(others_diag))

        if not r and return_sample_representative:
            representative = instance
        if doPositive:
            dist[r] = maxArea(instance, upper_triangle_only)
        else:
            dist[r] = -1 * maxArea(-1* instance, upper_triangle_only)
    dist.sort()

    if return_sample_representative:
        return dist, representative
    return dist


def explore(mtrx, i, j, visited, return_coords=False):

    if return_coords:
        res = list()
    stack = list()
    if not visited[i,j]:
        stack.append((i, j))
    val = 0
    while stack:
        i, j = stack.pop()
        if not visited[i,j] and mtrx[i,j] > 0:
            if return_coords:
                res.append((i, j))
            visited[i,j] = True
            val += mtrx[i,j]
            if return_coords:
                res.append((i, j))
            if i and not visited[i-1,j]:
                stack.append((i-1, j))
            if j+1  < mtrx.shape[1] and not visited[i, j+1]:
                stack.append((i, j+1))
            if i+1 < mtrx.shape[0] and not visited[i+1, j]:
                stack.append((i+1, j))
            if j and not visited[i, j-1]:
                stack.append((i, j-1))
    if return_coords:
        return val, res
    return val


def maxArea(mtrx, upper_triangle_only=False, return_coords=False,
        return_clusters=False):
    assert not upper_triangle_only or mtrx.shape[0] == mtrx.shape[1]

    visited = np.zeros(mtrx.shape, dtype = bool)
    if upper_triangle_only:
        pairs = combinations(xrange(mtrx.shape[0]), 2)
    else:
        pairs = product(xrange(mtrx.shape[0]), xrange(mtrx.shape[1]))

    mx_area = 0
    mx_pos = (-1, -1)
    if return_clusters:
        clusters = list()
    for i, j in pairs:
        if return_clusters:
            cur, coords = explore(mtrx, i, j, visited, return_coords=True)
            if coords:
                clusters.append((cur, coords))
        else:
            cur = explore(mtrx, i, j, visited)
        if cur > mx_area:
            mx_area = cur
            mx_pos = (i, j)
    if return_clusters:
        if return_coords:
            visited.fill(False)
            mx_area, coords = explore(mtrx, mx_pos[0], mx_pos[1], visited,
                    return_coords)
            return mx_area, coords, clusters
        else:
            return mx_area, clusters
    elif return_coords:
        visited.fill(False)
        return explore(mtrx, mx_pos[0], mx_pos[1], visited, return_coords)
    return mx_area


def doPlot(mtrx, sampled_mtrx, clusters, max_cluster_id, distribution, out,
        upper_triangle_only=False):

    nrows = 1
    ncols = 5
    if mtrx.shape[1] > mtrx.shape[0]:
        nrows, ncols = ncols, nrows
    fig = plt.figure(figsize=(ncols*5, nrows*5))

    # show positive cells of the matrix
    plt.subplot(nrows, ncols, 1)
    mask = np.array(mtrx > 0, dtype=float)
    if upper_triangle_only:
        mask[np.tril_indices(mtrx.shape[0], k=0, m=mtrx.shape[1])] = np.nan
    plt.imshow(mask, cmap=plt.cm.winter)
    plt.title('positive contact counts')

    # show identified connected areas
    plt.subplot(nrows, ncols, 2)
    clmtrx = np.zeros(mtrx.shape)
    clmtrx[:] = np.nan
    colors = np.random.permutation(len(clusters))
    for i, (_, coords) in enumerate(clusters):
        idx, idy = zip(*coords)
        clmtrx[idx, idy] = colors[i]
    plt.imshow(clmtrx, cmap=plt.cm.jet)
    plt.title('identified connected areas')

    # show max cluster
    plt.subplot(nrows, ncols, 3)
    maxmtrx = np.zeros(mtrx.shape)
    maxmtrx[:] = np.nan
    maxcoords = clusters[max_cluster_id][1]
    idx, idy = zip(*maxcoords)
    maxmtrx[idx, idy] = colors[i]
    plt.imshow(maxmtrx, cmap=plt.cm.jet)
    plt.title('area with largest volume')

    # show sampled matrix
    plt.subplot(nrows, ncols, 4)
    sampled_mask = np.array(sampled_mtrx > 0, dtype=float)
    if upper_triangle_only:
        sampled_mask[np.tril_indices(mtrx.shape[0], k=0, m=mtrx.shape[1])] = \
                np.nan
    plt.imshow(sampled_mask, cmap=plt.cm.winter)
    plt.title('positive contacts of sampled matrix')

    # show contact count distribution 
    plt.subplot(nrows, ncols, 5)
    plt.hist(distribution, 100, density=True)
    plt.axvline(clusters[max_cluster_id][0], label='observed cluster')
    plt.legend()
    plt.title('empirical dist. of contact volumes')

    plt.savefig(out, format='pdf')


def findNewestCluster(cluster):

    observed = set()
    observed.add(cluster)

    while cluster.new_cluster and cluster.new_cluster not in observed:
        cluster = cluster.new_cluster
        observed.add(cluster)
    if cluster.new_cluster:
        raise Exception, 'Circular path'
    return cluster


def constructMasks(sel_baited1, sel_baited2, upper_triangle_only):

    baited_mask = np.zeros((sel_baited1.size, sel_baited2.size), dtype=bool)
    others_mask = np.zeros((sel_baited1.size, sel_baited2.size), dtype=bool)
    for i in xrange(sel_baited1.size):
        baited_mask[i, :] = sel_baited1[i] & sel_baited2
        others_mask[i, :] = sel_baited1[i] ^ sel_baited2

    if upper_triangle_only:
        sel_lower_triangle = np.tril_indices(sel_baited1.size, k=0,
                m=sel_baited2.size)
        baited_mask[sel_lower_triangle] = False
        others_mask[sel_lower_triangle] = False

    return  baited_mask, others_mask


if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('hic_map', type=file, help='(Ca-)HI-C map in COOL format')
    parser.add_argument('coord1', type=str,
            help='coordinates of submatrix that is subject to p-value calculation')
    parser.add_argument('coord2', type=str, nargs='?',
            help='coordinates of submatrix that is subject to p-value calculation')
    parser.add_argument('-t', '--test', default='both', choices=('positive',
        'negative', 'both'), help='choose which values to test')
    parser.add_argument('-r', '--repeats', type=int, default=10000,
            help='Number of repeats used in sampling the ' + \
                    'empirical distribution')
    parser.add_argument('-p', '--plot', action='store',
            type=FileType('w'), help='plot identified areas, area with ' + \
                    'largest volume, contact count distribution, etc')
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
        m = (submtrx.shape[0] * (submtrx.shape[0]-1))/2
        LOG.info('selected coordinates describe a symmetric submatrix, ' + \
                'therefore only using values in its upper triangle ' + \
                '(omitting the main diagonal), to revent bias in ' + \
                'computing p-values')
    else:
        m = submtrx.shape[0] * submtrx.shape[1]
        LOG.info('gene cluster matrix is assumed to be not symmetric, ' + \
                'computing p-values on the sum of *all* values')

    LOG.info('computing baited/ other-end pools for sampling')
    pools_baited, pools_others = getSamplePools(mtrx, sel_baited, sel_coord1,
        sel_coord2)

    baited_mask, others_mask = constructMasks(sel_baited[sel_coord1],
            sel_baited[sel_coord2], isSymmetric)

    if args.test in ('positive' , 'both'):
        LOG.info(('computing empirical distribution of positive contact ' + \
                'areas based on sampling with %s repeats') %args.repeats)
        dist, sampled_mtrx = sampleContactAreas(baited_mask, others_mask,
                pools_baited, pools_others, args.repeats, isSymmetric,
                return_sample_representative=True)
        LOG.info('computing size of largest positive region within %s%s' %(
            args.coord1, args.coord2 and ', ' + args.coord2 or ''))
        submtrx_val, region, clusters = maxArea(submtrx, isSymmetric,
                return_coords=True, return_clusters=True)
        pval = 1-(np.searchsorted(dist, submtrx_val)/float(args.repeats))
        idx, idy = zip(*region)
        max_coord1 = set(coord1[x] for x in idx)
        max_coord2 = set(coord2[x] for x in idy)
        print('high contact region\t%s\t%ssignificant\t%.5f\t%s\t%s\t%s' %(\
                pval, pval > args.alpha and 'not ' or '', submtrx_val,
                len(region), len(region)/float(m),
                condenseCoords(sorted(max_coord1.union(max_coord2)), coords)))

        if args.plot:
            LOG.info('drawing figures and plot to %s' %args.plot.name)
            max_cluster_id = next(i for i, (x, _)  in enumerate(clusters) if x ==
                    submtrx_val)
            doPlot(submtrx, sampled_mtrx, clusters, max_cluster_id, dist,
                    args.plot, isSymmetric)

    if args.test in ('negative', 'both'):
        LOG.info(('computing empirical distribution of negative contact ' + \
                'areas based on sampling with %s repeats') %args.repeats)
        dist = sampleContactAreas(baited_mask, others_mask, pools_baited,
                pools_others, args.repeats, isSymmetric, doPositive=False)
        LOG.info('computing size of largest negative region within %s%s' %(
            args.coord1, args.coord2 and ', ' + args.coord2 or ''))
        submtrx_val, region = maxArea(-1 * submtrx, isSymmetric,
                return_coords=True)
        submtrx_val *= -1
        pval = np.searchsorted(dist, submtrx_val, side='right')/ \
                float(args.repeats)
        idx, idy = zip(*region)
        max_coord1 = set(coord1[x] for x in idx)
        max_coord2 = set(coord2[x] for x in idy)
        print('low contact region\t%s\t%ssignificant\t%.5f\t%s\t%s\t%s' %(\
                pval, pval > args.alpha and 'not ' or '', submtrx_val,
                len(region), len(region)/float(m),
                condenseCoords(sorted(max_coord1.union(max_coord2)), coords)))

    LOG.info('DONE')
