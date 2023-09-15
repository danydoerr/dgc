#!/usr/bin/env python3

from sys import stdout, stderr, exit
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF
from itertools import chain
import logging
import os
if not os.environ.get('DISPLAY', None):
    import matplotlib; matplotlib.use('Agg')
import matplotlib.pylab as plt

from calculatePvalue import readHiCMapTRV, readRMap, readBaitmapIDs


LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)


if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('restriction_digest', type=str, 
            help='restriction fragments file in RMAP format')
    parser.add_argument('baitmap', type=str, 
            help='baitmap from the capture part of CHI-C')
    parser.add_argument('hic_map', type=str, help='CHI-C map')
    parser.add_argument('genecluster_map', type=str, help='CHI-C map of gene cluster')
    args = parser.parse_args()

    # setup logging
    ch = logging.StreamHandler(stderr)
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))
    LOG.addHandler(ch)


    LOG.info('reading restriction site information file %s' %args.restriction_digest)
    rMap = readRMap(open(args.restriction_digest))
    rMapDict = dict((frag[3], frag[:3]) for frag in rMap)
    LOG.info('reading file contained probed fragments %s' %args.baitmap)
    probeIDs = readBaitmapIDs(open(args.baitmap))
    probedFragments = set('%s-%s'%rMapDict[x][:2] for x in probeIDs)

    LOG.info('reading CHi-C file %s' %args.hic_map)
    mtrx, _, _, _ = readHiCMapTRV(open(args.hic_map), probedFragments)

    LOG.info('reading genecluster file %s' %args.genecluster_map)
    gcMtrx, _, _, _ = readHiCMapTRV(open(args.genecluster_map), probedFragments)

    data = list(chain(*map(lambda x: filter(None, x),  mtrx)))
    gcData = list(chain(*map(lambda x: filter(None, x),  gcMtrx)))


    plt.figure()
    plt.hist([data, gcData], 100, density=True)
    plt.savefig(stdout, format='pdf')

