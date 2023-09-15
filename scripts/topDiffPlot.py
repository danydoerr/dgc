#!/usr/bin/env python3

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF, \
        FileType
from sys import stdout, stderr, exit, maxint
from os.path import dirname, join, isabs, basename, relpath
from itertools import izip, chain, repeat, compress
from cStringIO import StringIO
import logging
import csv
import re

import os, locale
locale.setlocale(locale.LC_ALL, 'en_US')
if not os.environ.get('DISPLAY', None):
    import matplotlib; matplotlib.use('Agg')

from matplotlib import pylab as plt
from scipy.stats import pearsonr
import numpy as np
import cooler as clr

from hic import assimilateMatrices, parseCoords


BREWER_COL = ['blues-%s-seq-3', 'greens-%s-seq-3', 'oranges-%s-seq-3', 'purples-%s-seq-3', \
        'reds-%s-seq-3', 'purd-%s-seq-3']

BREWER_COL_RANGE = range(3,11)

LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)


def condenseCoords(coords, merge=0, return_map=False):
    res = list()
    coord2condensed = dict()
    for chrx, start, end in coords:
        if not res or res[-1][0] != chrx or start-res[-1][2]-merge > 0:
            res.append([chrx, start, end])
        else:
            res[-1][2] = end
        if return_map:
            coord2condensed[(chrx, start, end)] = len(res)-1

    res = map(tuple, res)
    if return_map:
        return res, coord2condensed
    return res


def coarsenResolution(mtrx, coords, factor):
    """ coarsens the resolution along the x-axis of the matrix """
    assert mtrx.shape[1] == len(coords)
    res_mtrx = list()
    res_coords = list()
    i = 0
    while i < mtrx.shape[1]:
        j = i+1
        while j < (i/factor+1)*factor and j+1 < mtrx.shape[1] and coords[i][0] == \
                coords[j][0]:
            j += 1
        res_mtrx.append(np.median(mtrx[:, i:j], axis=1))
        res_coords.append(coords[i][:2] + (coords[j-1][2], ))
        i = j
    return np.array(res_mtrx).T, res_coords


def computeCorrelation(mtrx1, mtrx2):

    assert mtrx1.shape == mtrx2.shape

    res = np.empty((mtrx1.shape[0], 2))
    n = mtrx1.shape[1]
    for i in xrange(mtrx1.shape[0]):
        c, p = pearsonr(mtrx1[i, :], mtrx2[i, :])
        res[i, :] = (c, p)
    return res


def computeTopContacts(hic_maps, sel, coords, quantile, topx, resolution):

    if resolution > coords[0][2]-coords[0][1]:
        roughness_factor = resolution/(coords[0][2] - coords[0][1])
        mtrx1, coarse_coords = coarsenResolution(hic_maps[0][0][sel, :], coords,
                roughness_factor)
        mtrx2, _ = coarsenResolution(hic_maps[1][0][sel, :], coords,
                roughness_factor)
    else:
        coarse_coords = coords
        mtrx1 = hic_maps[0][0][sel, :]
        mtrx2 = hic_maps[1][0][sel, :]

    sel_coords = list(compress(coords, sel))

    corrs = computeCorrelation(mtrx1, mtrx2)
    A = np.flip(np.argsort(corrs[:, 1]))
    topx_changes = A[:topx]
    topx_changes = topx_changes[corrs[topx_changes, 1] > 0]
    diff_mtrx = hic_maps[0][0][sel, :] - hic_maps[1][0][sel, :]
    t2, t1 = np.quantile(diff_mtrx, (quantile/2, 1-quantile/2))
    top_contact_mtrx1 = diff_mtrx[topx_changes] > t1
    top_contact_mtrx2 = diff_mtrx[topx_changes] < t2

    top_contacts1 = list()
    top_contacts2 = list()
    for i, x in enumerate(topx_changes):
        top_contacts1.extend(zip(
            repeat(sel_coords[x]),
            compress(coarse_coords, top_contact_mtrx1[i]),
            -diff_mtrx[x, top_contact_mtrx1[i]]))
        top_contacts2.extend(zip(
            repeat(sel_coords[x]),
            compress(coarse_coords, top_contact_mtrx2[i]),
            diff_mtrx[x, top_contact_mtrx2[i]]))

    return top_contacts1, top_contacts2


def writeLinks(top_changes, links_out, circos_out0, circos_out1, radius,
        baited_coords=None, coord2region=None):

    if baited_coords == None:
        baited_coords = ()
    for c, changes in enumerate(top_changes):
        # sort by intensity
        links = sorted(changes, key=lambda x: (x[-1], ) + x[:-1])
        if not links:
            continue
        #
        # map intensity to 5 categories
        #
        intensities = np.array(map(lambda x: x[2], links))
        # six cuts make 5 categories
        thresholds = np.quantile(intensities, np.linspace(0, 1, 5))
        #category = np.ones(len(links), dtype=int)
        category = np.empty(intensities.size, dtype=int)
        pt = 0
        for i, t in enumerate(thresholds):
            category[np.logical_and(intensities > pt, intensities <= t)] = 5-i
            pt = t

        for k, (ci, cj, intensity) in enumerate(links):
            isBaited = (ci in baited_coords) and cj in (baited_coords)
            print(('%s %s %s %s %s %s condition=%s,intensity=%s'+\
                    ',bait_to_bait=%s%s') %( ci[0].lower(), ci[1], ci[2], \
                    cj[0].lower(), cj[1], cj[2], c, category[k], isBaited and \
                    1 or 0, coord2region != None and ',region=%s' % \
                    coord2region.get(ci, -1) or ''), file=links_out)

    for i, circos_out in enumerate((circos_out0, circos_out1)):
        print('<links>\n<link>', file=circos_out)
        print('file             = %s' %relpath(links_out.name,
                dirname(circos_out.name)), file=circos_out)
        print('radius          = %sr' %radius, file=circos_out)
        print('bezier_radius    = 0r', file=circos_out)
        print('ribbon           = yes', file=circos_out)
        print('flat             = yes', file=circos_out)
        print('show             = yes', file=circos_out)
        print('<rules>\n<rule>', file=circos_out)
        print('condition        = 1', file=circos_out)
        if coord2region:
            print('color            = eval(\'region\' . ' + \
                    'var(region))', file=circos_out) # . \'_a\' . var(intensity))'
        else:
            print('color            = eval(\'condition\' . ' + \
                'var(condition) . \'_a\' . var(intensity))', file=circos_out)
        print('flow             = continue', file=circos_out)
        print('</rule>\n<rule>', file=circos_out)
        print('condition        = !(var(condition) eq "%s")' %i, file=circos_out)
        print('show             = no', file=circos_out)
        print('</rule>\n</rules>', file=circos_out)
        print('</link>\n</links>', file=circos_out)


def writeKaryotype(chromsizes, karyotype_out, circos_out0, circos_out1):
    for c, s in chromsizes:
        # chr - ID LABEL START END COLOR
        print('chr - %s %s 0 %s ideogram.chr%s' %(c, c, s, c), file=karyotype_out)

    for circos_out in (circos_out0, circos_out1):
        print('karyotype = %s' %relpath(karyotype_out.name, dirname(circos_out0.name)), file=circos_out)
        print('chromosomes_units= 1000000', file=circos_out)
        print('chromosomes_display_default = no', file=circos_out)
        print('chromosomes      = %s' %(';'.join(map(lambda x:
            x[0], chromsizes))), file=circos_out)


def writeZoomRegions(coords, merge_window, chromsizes, scale, out):

    chrom2size = dict(chromsizes)
    print >> out, '<zooms>'
    i = 0
    while i < len(coords):
        if not i or coords[i-1][0] != coords[i][0] or \
                coords[i][1]-coords[i-1][2] > merge_window:
            if i:
                print >> out, 'end              = %sb' %min(coords[i-1][2] + \
                        merge_window/2, chrom2size[coords[i-1][0]])
                print >> out, '</zoom>'
            print >> out, '<zoom>'
            print >> out, 'smooth_distance  = 0.1r'
            print >> out, 'smooth_steps     = 3'
            print >> out, 'scale            = %s' %scale
            print >> out, 'chr              = %s' %coords[i][0]
            print >> out, 'start            = %sb' %max(0, coords[i][1] -
                    merge_window/2)
        i += 1
    print >> out, 'end   = %sb' %min(coords[-1][2] + merge_window/2, \
            chrom2size[coords[-1][0]])
    print >> out, ' </zoom>'
    print >> out, '</zoom>'


def writeExpressionPlot(expression_file, min_val, max_val, outer_radius,
        inner_radius, out, alpha=0.05):

    print >> out, '<plot>'
    print >> out, 'type             = histogram'
    print >> out, 'file             = %s' %relpath(expression_file.name,
            dirname(out.name))
    print >> out, 'r1               = %sr' %outer_radius
    print >> out, 'r0               = %sr' %inner_radius
    print >> out, 'max              = %s' %max_val
    print >> out, 'min              = %s' %min_val
    print >> out, 'stroke_type      = bin'
    print >> out, 'thickness        = 1'
    print >> out, 'extend_bin       = no'
    print >> out, 'max_gap          = 100000'
#    print >> out, 'color            = vdgrey'
    print >> out, '<backgrounds>\n<background>'
    #print >> out, 'color = vvlgrey'
    print >> out, 'color = white'
    print >> out, '</background>\n</backgrounds>'
    print >> out, '<axes>\n<axis>'
    print >> out, 'spacing          = 0.1r'
    print >> out, 'color            = lgrey'
    print >> out, 'thickness        = 2'
    print >> out, '</axis>\n</axes>'
    print >> out, '<rules>\n<rule>'
#    print >> out, '#condition  = var(value) < 0'
#    print >> out, 'condition  = 1'
#    print >> out, '#fill_color = lred'
#    print >> out, ('fill_color = eval(sprintf("spectral-9-div-%%d",' + \
#            'remap_int(var(value),%s,%s,1,9)))') %(min_val, max_val)
    print >> out, 'condition        = var(value) < 0 && var(pvalue) <= %s' %alpha
    print >> out, 'color            = spectral-9-div-1'
    print >> out, 'fill_color       = spectral-9-div-1'
    print >> out, '</rule>\n<rule>'
    print >> out, 'condition        = var(pvalue) > %s' %alpha
    print >> out, 'color            = spectral-9-div-3'
    print >> out, 'fill_color       = spectral-9-div-3'
    print >> out, 'show             = yes'
    print >> out, '</rule>\n<rule>'
    print >> out, 'condition        = var(value) > 0 && var(pvalue) <= %s' %alpha
    print >> out, 'color            = spectral-9-div-7'
    print >> out, 'fill_color       = spectral-9-div-7'
    print >> out, '</rule>\n</rules>'
    print >> out, '</plot>'


def writeHighlights(highlight_regions, highlight_out, circos_out0,
        circos_out1, outer_radius, inner_radius, color_prefix,
        enumerateColors=True, z=1):

    # assumes coords are sorted
    for circos_out in (circos_out0, circos_out1):
        print >> circos_out, '<highlight>'
        print >> circos_out, 'stroke_thickness = 0'
        print >> circos_out, 'show             = yes'
        print >> circos_out, 'z                = %s' %z
        print >> circos_out, 'file             = %s' %relpath(
                highlight_out.name, dirname(circos_out0.name))
        print >> circos_out, 'r0               = %sr' %outer_radius
        print >> circos_out, 'r1               = %sr' %inner_radius
        print >> circos_out, '</highlight>'

    for i, (chrx, start, end) in enumerate(highlight_regions):
        print >> highlight_out, '%s %s %s fill_color=%s%s' %(chrx, start,
                end, color_prefix, enumerateColors and str(i) or '')


def writeColors(out, chrs, highlights=None, region2id=None):

    print >> out, '<colors>'
    print >> out, 'condition0       = black'
    print >> out, 'condition1       = black'

    c_max = len(BREWER_COL_RANGE) * len(BREWER_COL)
    if len(chrs) > c_max:
        print >> stderr, '!! Too many chromosomes. Not enough colors ' + \
                'available to color all links! Exiting'
        exit(1)

    for x, c in enumerate(chrs):
        r = x / len(BREWER_COL)
        i = x % len(BREWER_COL)
        print >> out, 'ideogram.chr%s = %s' %(x, BREWER_COL[i]
                %BREWER_COL_RANGE[r])
        print >> out, 'band.chr%s = vvlgrey' %x
        # BREWER_COL[i] %BREWER_COL_RANGE[(r+3)%7]

    if highlights:
        for i in xrange(len(highlights)):
            print >> out, 'highlight%s = lgrey' %i

    if region2id:
        N = len(region2id)
        for rid, col in zip(sorted(region2id.values()), np.linspace(0, 360, N)):
            print >> out, 'region%s = hue%03d' %(rid, col)
            print >> out, 'band%s = hue%03d' %(rid, col)

    print >> out, '</colors>'


def writeAnnotationPlot(annotation_file, outer_radius, inner_radius, out, z=1):

    print >> out, '<plot>'
    print >> out, 'type             = tile'
    print >> out, 'file             = %s' %relpath(annotation_file.name,
            dirname(out.name))
    print >> out, 'r1               = %sr' %outer_radius
    print >> out, 'r0               = %sr' %inner_radius
    print >> out, 'layers           = 4'
    print >> out, 'margin           = 0.001u'
    print >> out, 'thickness        = 15'
    print >> out, 'padding          = 5'
    print >> out, 'layers_overflow  = hide'
    print >> out, 'orientation      = out'
    print >> out, 'stroke_thickness = 0'
    print >> out, 'stroke_color     = grey'
    print >> out, 'color            = grey'
    print >> out, 'z                = %s' %z
    print >> out, '<backgrounds>\n<background>'
    print >> out, 'color            = vvlgrey'
    print >> out, '</backgrounds>\n</background>'
    print >> out, '</plot>'

def writeHighlightGenes(genes, annotation_file, outer_radius, inner_radius,
        highlight_out, circos_out0, circos_out1, z=1):

    genes = set(genes)
    PAT_LOCUS=re.compile('.*locus=([^,]+)')

    for line in csv.reader(open(annotation_file.name), delimiter=' '):
        m = PAT_LOCUS.match(line[-1])
        if m and m.group(1) in genes:
            print >> highlight_out, ' '.join(line[:3] + [m.group(1)])

    for circos_out in (circos_out0, circos_out1):
        print >> circos_out, '<plot>'
        print >> circos_out, 'type             = text'
        print >> circos_out, 'file             = %s' %relpath(
                highlight_out.name, dirname(circos_out.name))
        print >> circos_out, 'r1               = %sr' %outer_radius
        print >> circos_out, 'r0               = %sr' %inner_radius
        print >> circos_out, 'show_links       = no'
        print >> circos_out, 'label_size       = 18p'
        print >> circos_out, 'label_font       = condensed'
        print >> circos_out, 'z                = %s' %z
        print >> circos_out, '<backgrounds>\n<background>'
        print >> circos_out, 'color            = white'
        print >> circos_out, '</backgrounds>\n</background>'
        print >> circos_out, '</plot>'

def writeGeneralConf(image_path, out, zoom_coord=None):

    print >> out, '<ideogram>'
    print >> out, '<spacing>'
    print >> out, 'default          = 0.005r'
    print >> out, '</spacing>'
    print >> out, 'radius           = 0.90r'
    print >> out, 'thickness        = 0p'
    print >> out, 'stroke_thickness = 0'
    print >> out, 'show_bands       = yes'
    print >> out, 'fill_bands       = yes'
    print >> out, 'band_stroke_thickness = 0'
    print >> out, 'band_transparency= 1'
    print >> out, 'show_label       = yes'
    print >> out, 'label_font       = default'
    print >> out, 'label_radius     = 1.1r'
    print >> out, 'label_size       = 60'
    print >> out, 'label_parallel   = no'
    print >> out, '</ideogram>'
    print >> out, 'show_ticks       = yes'
    print >> out, 'show_tick_labels = yes'
    print >> out, '<ticks>'
    print >> out, 'radius           = 1r'
    print >> out, 'color            = black'
    print >> out, 'thickness        = 2p'
    print >> out, 'multiplier       = 1e-6'
    print >> out, 'label_separation = 10p'
    print >> out, 'format           = %dMb'
    print >> out, '<tick>'
    print >> out, 'spacing          = 0.1u'
    print >> out, 'size             = 5p'
    print >> out, '</tick>'
    print >> out, '<tick>'
    print >> out, 'spacing          = 1u'
    print >> out, 'size             = 10p'
    print >> out, '</tick>'
    print >> out, '<tick>'
    print >> out, 'spacing          = 5u'
    print >> out, 'size             = 15p'
    print >> out, 'show_label       = yes'
    print >> out, 'label_size       = 30p'
    print >> out, 'label_offset     = 10p'
    print >> out, '</tick>'
    if zoom_coord != None:
        order = int(10**np.floor(np.log10(zoom_coord[2]-zoom_coord[1])-1))
        start = (zoom_coord[1]/order) *  order
        print >> out, '<tick>'
        print >> out, 'position = %s' %','.join(map(str, xrange(start,
            zoom_coord[2], order)))
        print >> out, 'size             = 10p'
        print >> out, 'show_label       = no'
        print >> out, '</tick>'

        ending = ['bp', 'kb', 'Mb']
        canon_order = int(np.log10(max(start, 1)))/3
        for i in xrange(zoom_coord[1]/(order*10) *  order*10,
                zoom_coord[2]+1, order*10):
            print >> out, '<tick>'
            print >> out, 'chromosomes_display_default = no'
            print >> out, 'chromosomes      = %s' %zoom_coord[0]
            print >> out, 'position         = %s' %i
            print >> out, 'size             = 15p'
            print >> out, 'show_label       = yes'
            print >> out, 'label            = %.1f%s'%(i/(10.**(canon_order*3)),
                    ending[canon_order])
            print >> out, 'label_size       = 30p'
            print >> out, 'label_offset     = 10p'
            print >> out, 'force_display    = yes'
            print >> out, '</tick>'

    print >> out, '</ticks>'
    print >> out, '<image>\nfile  = %s' %image_path
    print >> out, 'dir   = .\npng   = yes\nradius         = 1500p'
    print >> out, 'angle_offset      = -90'
    print >> out, 'auto_alpha_colors = yes\nauto_alpha_steps  = 5'
    print >> out, 'background = white\n</image>'
    print >> out, '<<include etc/colors_fonts_patterns.conf>>'
    print >> out, '<<include etc/colors.brewer.conf>>'
    print >> out, '<<include %s>>' %relpath('etc/housekeeping.conf',
            dirname(out.name))
    print >> out, 'data_out_of_range* = warning'


if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('cool_file', type=file, nargs=2,
            help='(Ca-) HI-C map(s) in COOLER format')
    parser.add_argument('annotation_file', type=file,
            help='genome annotation file in circos format')
    parser.add_argument('topx', type=int,
            help='visualize top x changes')
    parser.add_argument('-o', '--out_dir', type=str, default='.',
            help = 'output directory')
    parser.add_argument('-e', '--expression', type=file,
            help='gene epression file in circos format')
    parser.add_argument('-q', '--quantile', type=float, default=0.05,
            help = 'threshold for significant interaction')
    parser.add_argument('-g', '--highlight_genes', type=str, nargs='*',
            help='do not label all genes in the region, but only those' + \
                    'specified')
    parser.add_argument('-i', '--highlight_regions', type=str, nargs='*',
            help='highlight all regions specified here in condensed format')
    parser.add_argument('-c', '--chromosomes_only', type=str, nargs='+',
            help = 'only display specified chromosomes')
    parser.add_argument('-r', '--resolution', type=int, default=0, help =
            'resolution of the differential analysis')
    parser.add_argument('-z', '--zoom', type=str, help =
            'give a 1/3 zoom to a specific region')

    args = parser.parse_args()

    #
    # setup logging
    #
    ch = logging.StreamHandler(stderr)
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(logging.Formatter('%(levelname)s\t%(asctime)s\t%(message)s'))
    LOG.addHandler(ch)

    #
    # read Hi-C maps
    #
    hic_maps = list()
    for f in args.cool_file:
        LOG.info('loading (Ca-) Hi-C map %s' %f.name)
        c = clr.Cooler(f.name)
        if args.chromosomes_only:
            unknowns = set(args.chromosomes_only).difference(c.chromnames)
            if unknowns:
                LOG.fatal('unknown chromosomes %s, exiting' %(
                    ', '.join(unknowns)))
                exit(1)

        mtrx = c.matrix()[:]
        mtrx[np.isnan(mtrx)] = 0
        bins = c.bins()[:]
        coords = zip(bins['chrom'].values, bins['start'].values,
                bins['end'].values)
        hic_maps.append((mtrx, coords, list(coords), bins, c))

    chromnames = hic_maps[0][-1].chromnames
    chromsizes = hic_maps[0][-1].chromsizes
    if args.chromosomes_only:
        chromnames = filter(lambda x: x in args.chromosomes_only, chromnames)
        chromsizes = list(compress(chromsizes, map(lambda x: x in chromnames,
            hic_maps[0][-1].chromnames)))

    LOG.info('make sure matrices have the same size...')
    hic_maps, common = assimilateMatrices(hic_maps,
            return_common=True)
    coords = hic_maps[0][1]

    sel = np.ones(len(coords), dtype = 'bool')
    isBaited = False
    for i, hic_map in enumerate(hic_maps):
        bins = hic_map[3]
        if 'baited' in bins.columns:
            sel = np.logical_and(sel, bins['baited'][common[i][1]])
            isBaited = True

    LOG.info('compute correlation of the count distribution of each genomic ' +\
            'segment between the two Hi-C maps')
    if not args.resolution:
        args.resolution = coords[0][2]-coords[0][1]

    top_contacts1, top_contacts2 = computeTopContacts(hic_maps, sel, coords,
            args.quantile, args.topx, args.resolution)

    zoom_coord = None
    if args.zoom:
        zoom_coord = parseCoords(args.zoom)
        intersects = lambda x, y: x[0] == y[0] and ((x[1] < y[1] and y[1] <
            x[2]) or (y[1] < x[1] and x[1] < y[2]))
        top_contacts1 = filter(lambda x: intersects(x[0], zoom_coord) or
                intersects(x[1], zoom_coord), top_contacts1)
        top_contacts2 = filter(lambda x: intersects(x[0], zoom_coord) or
                intersects(x[1], zoom_coord), top_contacts2)

    #
    # write configuration
    #
    file_prefix0 = join(args.out_dir, basename( \
            args.cool_file[0].name).rsplit('.cool', 1)[0])
    file_prefix1 = join(args.out_dir, basename( \
            args.cool_file[1].name).rsplit('.cool', 1)[0])

    # open file handles for circos configuraiton
    circos_out0 = open('%s.circos.conf' %file_prefix0 , 'w')
    circos_out1 = open( '%s.circos.conf' %file_prefix1, 'w')

    # write karyotypes
    karyotype_out = open(join(args.out_dir, 'karyotype.txt'), 'w')
    writeKaryotype(zip(chromnames, chromsizes), karyotype_out, circos_out0,
            circos_out1)
    writeGeneralConf('%s.png' %file_prefix0, circos_out0, zoom_coord)
    writeGeneralConf('%s.png' %file_prefix1, circos_out1, zoom_coord)

    outer_band_r = 0.935
    inner_radius = args.expression and 0.65 or outer_band_r

    # mapping of regions (for coloring links) to IDs
    baited_coords = list()
    if isBaited:
        baited_coords = list(compress(coords, sel))
        baited_coords_condensed, coord2region = condenseCoords(baited_coords,
                merge=50000, return_map=True)
        region2id = dict(zip(baited_coords_condensed, xrange(len(
            baited_coords_condensed))))
    else:
        coord2region = region2id = dict(zip(sorted(set(chain(*(zip(
            *top_contacts1)[:2] + zip(*top_contacts2)[:2])))),
            xrange(len(coords))))

    # write links
    links_out = open(join(args.out_dir, 'links.txt'), 'w')
    writeLinks((top_contacts1, top_contacts2), links_out, circos_out0,
            circos_out1, inner_radius, baited_coords=set(baited_coords),
            coord2region=coord2region)

    writeColors(circos_out0, chromnames, highlights=args.highlight_regions,
            region2id=region2id)
    writeColors(circos_out1, chromnames, highlights=args.highlight_regions,
            region2id=region2id)

    print >> circos_out0, '<highlights>'
    print >> circos_out1, '<highlights>'

    highlight_out = open(join(args.out_dir, 'chromosomes.txt'), 'w')
    writeHighlights(map(lambda x: (x[0], 0, x[1]), zip(chromnames, chromsizes)),
            highlight_out, circos_out0, circos_out1, 1.0, outer_band_r, 'band.chr', 1)
    highlight_out.close()

    if args.highlight_regions:
        highlight_regions = map(parseCoords, args.highlight_regions)
        highlight_out = open(join(args.out_dir, 'highlights.txt'), 'w')
        writeHighlights(highlight_regions, highlight_out, circos_out0,
                circos_out1, 1.0, inner_radius, 'highlight', 2)
        highlight_out.close()

    if isBaited:
        highlight_out = open(join(args.out_dir, 'baited_regions.txt'), 'w')
        writeHighlights(baited_coords_condensed, highlight_out, circos_out0,
                circos_out1, 1.0, inner_radius, 'region', 3)
        highlight_out.close()
        highlight_out.close()


    print >> circos_out0, '</highlights>'
    print >> circos_out1, '</highlights>'

    if zoom_coord != None:
        N = sum(chromsizes)
        zoom_factor = N/((zoom_coord[2]-zoom_coord[1]) * 2.)
        writeZoomRegions([zoom_coord], 0, zip(chromnames, chromsizes),
                zoom_factor, circos_out0)
        writeZoomRegions([zoom_coord], 0, zip(chromnames, chromsizes),
                zoom_factor, circos_out1)

    print >> circos_out0, '<plots>'
    print >> circos_out1, '<plots>'
    writeAnnotationPlot(args.annotation_file, outer_band_r+0.05, outer_band_r,
            circos_out0, z=4)
    writeAnnotationPlot(args.annotation_file, outer_band_r+0.05, outer_band_r,
            circos_out1, z=4)
    if args.highlight_genes:
        highlight_out = open(join(args.out_dir, 'highlight_genes.txt'), 'w')
        writeHighlightGenes(args.highlight_genes, args.annotation_file,
                outer_band_r, outer_band_r-0.1, highlight_out, circos_out0,
                circos_out1, z=5)
        highlight_out.close()
    if args.expression:
        expression_data = np.array(tuple(float(row[3]) for row in
            csv.reader(args.expression, delimiter=' ')))
        min_val, max_val = np.quantile(expression_data, (0.01, 0.99))
        writeExpressionPlot(args.expression, min_val, max_val, outer_band_r -
                0.1, inner_radius+.02, circos_out0)
        writeExpressionPlot(args.expression, min_val, max_val, outer_band_r -
                0.1, inner_radius+.02, circos_out1)
    print >> circos_out0, '</plots>'
    print >> circos_out1, '</plots>'


