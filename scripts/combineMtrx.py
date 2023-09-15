#!/usr/bin/env python3

from sys import stdout, stderr, exit
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF, \
        FileType
from functools import partial
from itertools import izip, chain, combinations, compress, repeat
import logging

from hic import assimilateMatrices

from pandas import DataFrame, Series
import cooler as clr
import numpy as np

LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)

INTRA_CUTOFF = 1500000


def renormalize(mtrx, select_bait, coords, intra_cutoff, isdiff=False):
    """ this method assumes that the coordinates are consecutive, i.e., no
    coordinate is missing """

    resolution = coords[0][2] - coords[0][1]
    res = np.empty(mtrx.shape)
    idx = np.arange(mtrx.shape[0])
    res[idx, idx] = np.nan

    chr_mask = np.array([coords[i][0] == coords[i+1][0] for i in \
            xrange(len(coords)-1)])
    pos = np.argwhere(chr_mask == False)

    for i in xrange(1, mtrx.shape[0]):
        d = np.array(np.diag(mtrx, k=i))
        target = 1.-(i-1)/float(mtrx.shape[0])
        # perform intra-chromosomal max-distance cut-off
        if i*resolution > intra_cutoff:
            chr_mask[:] = False
            p = np.empty(0)

        nz_mask = (d != 0) & np.isfinite(d)
        # construct baited/non-baited mask
        bait_mask = select_bait[i:] & select_bait[:select_bait.size-i]
        #
        # four different cases:
        #   (baited, other-end) x (intra-chromosomal, inter-chormosomal)
        #
        mask = nz_mask & bait_mask & chr_mask[i-1:]
        baited_intra = d[mask]
        if baited_intra.size > 1:
            if isdiff:
                d[mask] = (baited_intra-baited_intra.mean()) / \
                        np.sqrt(baited_intra.var())
            else:
                m = np.sum(baited_intra)
                d[mask] *= target/m * baited_intra.size/float(d.size)

        mask = nz_mask & bait_mask & (chr_mask[i-1:] == False)
        baited_inter = d[mask]
        if baited_inter.size > 1:
            if isdiff:
                d[mask] = (baited_inter-baited_inter.mean()) / \
                        np.sqrt(baited_inter.var())
            else:
                m = np.sum(baited_inter)
                d[mask] *= target/m * baited_inter.size/float(d.size)

        mask = nz_mask & (bait_mask == False) & chr_mask[i-1:]
        othere_intra = d[mask]
        if othere_intra.size > 1:
            if isdiff:
                d[mask] = (othere_intra-othere_intra.mean()) / \
                        np.sqrt(othere_intra.var())
            else:
                m = np.sum(othere_intra)
                d[mask] *= target/m * othere_intra.size/float(d.size)

        mask = nz_mask & (bait_mask == False) & (chr_mask[i-1:] == False)
        othere_inter = d[mask]
        if othere_inter.size > 1:
            if isdiff:
                d[mask] = (othere_inter-othere_inter.mean()) / \
                        np.sqrt(othere_inter.var())
            else:
                m = np.sum(othere_inter)
                d[mask] *= target/m * othere_inter.size/float(d.size)

        res[idx[:-i], idx[i:]] = res[idx[i:], idx[:-i]] = d

        # update chromosome-mask
        for p in pos:
            if p+i < chr_mask.size:
                chr_mask[p+i] = False
    return res


def normavg(mtrcs_type1, mtrcs_type2):

    min_type1 = reduce(np.minimum, mtrcs_type1)
    max_type1 = reduce(np.maximum, mtrcs_type1)

    min_type2 = reduce(np.minimum, mtrcs_type2)
    max_type2 = reduce(np.maximum, mtrcs_type2)

    overlap = np.minimum(max_type1, max_type2)-np.maximum(min_type1, min_type2)
    above_zero = overlap > 0
    below_zero = overlap <= 0
    overlap[above_zero] /= np.minimum(
            max_type1[above_zero] - min_type1[above_zero],
            max_type2[above_zero] - min_type2[above_zero])
    overlap[above_zero] += 1.0
    overlap[below_zero] = 1.0

    zeros = (max_type1 == 0) & (max_type2 == 0)
    mtrx = (reduce(np.add, mtrcs_type1)/len(mtrcs_type1) - reduce(np.add,
        mtrcs_type2)/len(mtrcs_type2)) / overlap
    mtrx[zeros] = np.nan
    return mtrx


if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('-r', '--renormalize', action = 'store_true',
            help = 'renormalize diagonals for combined matrix')
    parser.add_argument('-1', '--hic_maps1', type=file, default=(), nargs='+',
            help='(Ca-) HI-C map(s) in COOLER format')
    parser.add_argument('-2', '--hic_maps2', type=file, default=(), nargs='+',
            help='(Ca-) HI-C map(s) in COOLER format')
    parser.add_argument('aggregator', type=str, choices = ('diff', 'min',
        'max', 'mean', 'normavg'),
        help='Aggregator function that maps two values from each ' + \
                'input map, respectively, to one value in the output matrix')
    parser.add_argument('out', type=FileType('w'),
            help='path to output file in COOLER format')
    args = parser.parse_args()

    #
    # setup logging
    #
    ch = logging.StreamHandler(stderr)
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(logging.Formatter('%(levelname)s\t%(asctime)s\t%(message)s'))
    LOG.addHandler(ch)


    if args.aggregator == 'normavg' and \
            (len(args.hic_maps1) < 1 or len(args.hic_maps2) < 1):
        print('chosen aggregator requires at least one matrix ' + \
                'of each type', file=stderr)
        exit(1)
    if args.aggregator == 'diff' and \
            (len(args.hic_maps1) != 1 or len(args.hic_maps2) != 1):
        print('chosen aggregator doesn\'t support more than two ' + \
                'matrices', file=stderr)
        exit(1)

    hic_maps = list()
    for f in chain(args.hic_maps1, args.hic_maps2):
        LOG.info('loading Hi-C map %s' %f.name)
        c = clr.Cooler(f.name)
        mtrx = c.matrix()[:]
        bins = c.bins()[:]
        coords = zip(bins['chrom'].values, bins['start'].values,
                bins['end'].values)
        hic_maps.append((mtrx, coords, list(coords), bins))

    LOG.info('make sure matrices have the same size...')
    hic_maps, common = assimilateMatrices(hic_maps, return_common=True)
    hic_maps1 = hic_maps[:len(args.hic_maps1)]
    hic_maps2 = hic_maps[len(args.hic_maps2):]
    coords = hic_maps[0][1]
    __one__ = lambda x, y, f: f(x[0], y[0])

    agg = dict((
        ('min', lambda x, y: reduce(np.minimum, chain(x, y))),
        ('max', lambda x, y: reduce(np.maximum, chain(x, y))),
        ('mean', lambda x, y: reduce(np.add, chain(x, y))/(len(x)+len(y))),
        ('diff', partial(__one__, f=lambda a, b: a-b)),
        ('normavg', normavg)))

    aggMtrx = apply(agg[args.aggregator], (map(lambda x: x[0], hic_maps1),
            map(lambda x: x[0], hic_maps2)))

    common_outlier = np.zeros(aggMtrx.shape[0], dtype=bool)
    hasOutlier = np.any(map(lambda x: 'outlier' in x[3].columns, hic_maps))
    if hasOutlier:
        LOG.info('merge outlier of all input files')
        for (_, _, _, bins), (c, _) in izip(hic_maps, common):
            if 'outlier' in bins.columns:
                common_outlier |= bins['outlier'].values[c]

    select_bait = np.ones(aggMtrx.shape[0], dtype=bool)
    hasBaited = np.all(map(lambda x: 'baited' in x[3].columns, hic_maps))
    if hasBaited:
        LOG.info('merge bait index of all input files')
        for (_, _, _, bins), (c, _) in izip(hic_maps, common):
            select_bait &= bins['baited'].values[c]

    if args.renormalize:
        LOG.info('renormalizing combined matrix')
        active = common_outlier == False
        select_bait = select_bait[active]
        coords = list(compress(coords, active))
        aggMtrx = renormalize(aggMtrx[np.ix_(active, active)], select_bait,
                coords, INTRA_CUTOFF, isdiff=args.aggregator in ('normavg',
                    'diff'))

    LOG.info('construct output cooler')
    bins = DataFrame(data=map(lambda x: x + (1.0, ), coords), columns =
            ('chrom', 'start', 'end', 'weight'))

    if hasOutlier and not args.renormalize:
        bins = bins.assign(outlier = Series(common_outlier))
    if hasBaited:
        bins = bins.assign(baited = select_bait)

    pixels = DataFrame(data=((x, y, aggMtrx[x,y]) for x,y in \
            chain(combinations(xrange(aggMtrx.shape[0]), 2), \
            zip(*repeat(np.arange(aggMtrx.shape[0]), 2))) if aggMtrx[x,y]), \
            columns=('bin1_id', 'bin2_id', 'count'))

    LOG.info('storing output cooler in %s' %args.out.name)
    clr.create_cooler(args.out.name, bins=bins, pixels=pixels, ordered=True,
            dtypes=dict((('bin1_id', int), ('bin2_id', int), ('count', float))))
    LOG.info('DONE')


