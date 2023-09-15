#!/usr/bin/env python3

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF
from sys import stdout, stderr, exit
from itertools import imap, product
import logging
import csv
import re

import os, locale
locale.setlocale(locale.LC_ALL, 'en_US')
if not os.environ.get('DISPLAY', None):
    import matplotlib; matplotlib.use('Agg')

from matplotlib import pylab as plt
from matplotlib import colors, patches
import numpy as np
import cooler as clr


from hic import condenseCoords, parseCondensedCoords, parseCoords, \
        equiSegmentRegions


LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)

PAT_CHR = re.compile('^(?:.*;)?chromosome=([^;]+)(?:;.*)?$', re.I)
PAT_NAME = re.compile('^(?:.*;)?Name=([^;]+)(?:;.*)?$', re.I)
PAT_LOCUS = re.compile('^(?:.*;)?locus_tag=([^;]+)(?:;.*)?$', re.I)
PAT_BIOTYPE = re.compile('^(?:.*;)?\w+_biotype=([^;]+)(?:;.*)?$', re.I)

HIGHLIGHT_COLOR = (158./255, 0, 150./255)
#HIGHLIGHT_COLOR = (140./255, 26./255, 1)
HIGHLIGHT_WIDTH = 3


def readGFF3Data(data, chrx, start, end, feature_types, bio_types=()):
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
        elif line[2] in feature_types:
            m = PAT_BIOTYPE.match(line[8])
            if not bio_types or (m and m.group(1).lower() in bio_types):
                names = [line[0]]
                m = PAT_NAME.match(line[8])
                if m:
                    names.append(m.group(1))
                m = PAT_LOCUS.match(line[8])
                if m:
                    names.append(m.group(1))
                anns.add((line[0], int(line[3]), int(line[4]), line[6], line[2],
                    tuple(names)))

    res = list()
    for ann in anns:
        chry = chrMap.get(ann[0], ann[0])
        if chry == chrx and ann[1] >= start and ann[2] <= end:
            res.append((chry, ) + ann[1:])

    # inefficient, but hey, data is not that large
    res.sort(key=lambda x: (x[0], x[1], -x[2]))

    return res


def drawCoordsX(axX, xoffset, yoffset, start, end):

    dim = int(np.log10(end-start))
    major_ticks = 10**dim
    minor_ticks = 10**(dim-1)

    line_size = 0.05
    # small ticks
    for p in xrange(int(np.ceil(start/float(minor_ticks))) * minor_ticks, end+1,
            minor_ticks):
        axX.plot([p-xoffset, p-xoffset], [yoffset, yoffset+line_size/2], 'k-',
                linewidth=1)
    # middle ticks
    for p in xrange(int((np.ceil(start/float(major_ticks))+ 0.5) *major_ticks),
            end+1, major_ticks):
        axX.plot([p-xoffset, p-xoffset], [yoffset, yoffset+line_size], 'k-',
                linewidth=1)

    for p in xrange(int(np.ceil(start/float(major_ticks))) * major_ticks, end+1,
            major_ticks):
        axX.plot([p-xoffset, p-xoffset], [yoffset, yoffset+line_size], 'k-',
                linewidth=2)
        axX.text(p-xoffset, yoffset+line_size + 0.02, locale.format("%d", p,
            grouping=True), fontsize=6, horizontalalignment='center')


def drawCoordsY(axY, xoffset, yoffset, start, end):

    dim = int(np.log10(end-start))
    major_ticks = 10**dim
    minor_ticks = 10**(dim-1)

    line_size = 0.1

    text_offset = 0.1
    xoffset += line_size+text_offset
    # small ticks
    for p in xrange(int(np.ceil(start/float(minor_ticks))) * minor_ticks, end+1,
            minor_ticks):
        axY.plot([xoffset-line_size/2, xoffset], [yoffset-p, yoffset-p], 'k-',
                linewidth=1)
    # middle ticks
    for p in xrange(int((np.ceil(start/float(major_ticks))+ 0.5) *major_ticks),
            end+1, major_ticks):
        axY.plot([xoffset-line_size, xoffset], [yoffset-p, yoffset-p], 'k-',
                linewidth=1)

    for p in xrange(int(np.ceil(start/float(major_ticks))) * major_ticks, end+1,
            major_ticks):
        axY.plot([xoffset-line_size, xoffset], [yoffset-p, yoffset-p], 'k-',
                linewidth=2)
        axY.text(xoffset-line_size-text_offset, yoffset-p, locale.format("%d", p,
            grouping=True), fontsize=6, rotation=90, verticalalignment='center')


def drawGenesX(axX, xoffset, yoffset, colannotations, f, max_stacks=15):

    lastpos = [float('-inf') for _ in xrange(max_stacks)]

    scaling = 0.006
    r = f.canvas.get_renderer()
    height = None
    for _, start, end, orient, atype, name in colannotations:
        if atype == 'gene':
            curpos = (start+end)/2-xoffset
            i = 0

            if orient == '+':
                axX.add_patch(patches.FancyArrow(start-xoffset, yoffset, end-start,
                    0, fill=False, overhang=1, width=0,
                    length_includes_head=True, head_length=500, head_width=0.1,
                    edgecolor='black'))
            else:
                axX.add_patch(patches.FancyArrow(end-xoffset, yoffset, start-end,
                    0, fill=False, overhang=1, width=0,
                    length_includes_head=True, head_length=500, head_width=0.1,
                    edgecolor='black'))

            if name:
                yoffsett = yoffset + 0.1
                t = axX.text(curpos, yoffsett, name, fontsize=6,
                        horizontalalignment='center')
                bb = t.get_window_extent(renderer=r)
                if height == None:
                    height = bb.height
                while i < len(lastpos) and lastpos[i] > bb.x0:
                    i += 1

                if i < len(lastpos):
                    yoffsett+= i*1.2*height*scaling
                    t.set_position([curpos, yoffsett])
                    lastpos[i] = bb.x1
                else:
                    # sent to nirvana
                    t.set_position([curpos, 100])

        elif atype == 'exon':
            continue
            axX.add_patch(patches.Rectangle((start-xoffset, yoffset), end-start,
                yoffset + 0.2, edgecolor='none', facecolor =
                colors.to_rgba('blue', 0.3)))
        else:
            print(('Warning: unknown annotation type %s, unable,' + \
                    ' to visualize') %atype, file=stderr)


def drawGenesY(axY, xoffset, yoffset, maxleny, rowannotations, f):

    r = f.canvas.get_renderer()
    lastpos = float('inf')
    for _, start, end, orient, atype, name in rowannotations:
        if atype == 'gene':
            if name:
                curpos = (start+end)/2-yoffset
                t = axY.text(xoffset+0.01, maxleny-curpos, name, fontsize=6,
                        verticalalignment='center')
                bb = t.get_window_extent(renderer=r)
                if lastpos < bb.y1:
                    # don't show / sent to nirvana
                    t.set_position((0.001, -10000))
                else:
                    lastpos = bb.y0

            if orient == '+':
                axY.add_patch(patches.FancyArrow(0.81,
                    maxleny-(start-yoffset), xoffset, start-end, fill=False,
                    overhang=1, width=0, length_includes_head=True,
                    head_length=500, head_width=0.2, edgecolor='black'))
            else:
                axY.add_patch(patches.FancyArrow(0.81,
                    maxleny-(end-yoffset), xoffset,
                    end-start, fill=False, overhang=1, width=0,
                    length_includes_head=True, head_length=500, head_width=0.2,
                    edgecolor='black'))
        elif atype == 'exon':
            continue
            axY.add_patch(patches.Rectangle((start-yoffset, 0.21), end-start,
                0.41, edgecolor='none', facecolo=colors.to_rgba('blue', 0.3)))
        else:
            LOG.info('unknown annotation type %s, unable, to visualize' %atype)


def drawMatrix(ax, offsetx, offsety, maxleny, colorFunc, mtrx,
        colcoords, rowcoords, outlier_cols, outlier_rows):
    for j in xrange(mtrx.shape[0]):
        for i in xrange(mtrx.shape[1]):
            if np.isfinite(mtrx[j][i]):
                ax.add_patch(patches.Rectangle((colcoords[i][1]-offsetx,
                    maxleny-(rowcoords[j][2]-offsety)),
                    colcoords[i][2]-colcoords[i][1],
                    rowcoords[j][2]-rowcoords[j][1], edgecolor='none',
                    facecolor=colorFunc(mtrx[j][i]), fill=True))
            if outlier_rows[j] or outlier_cols[i]:
                ax.add_patch(patches.Rectangle((colcoords[i][1]-offsetx,
                    maxleny-(rowcoords[j][2]-offsety)),
                    colcoords[i][2]-colcoords[i][1],
                    rowcoords[j][2]-rowcoords[j][1], edgecolor='none',
                    facecolor=(1, 0, 0, 0.5), fill=True))


def constructColorFunc(minval, maxval, colormap):
    return lambda x: colormap.colors[int(round(min(max((x-minval)/(maxval-minval), 0),
        1) * 255))]


def highlightRegions(ax, regions, colcoords, rowcoords):

    y_regions = list()
    x_regions = list()

    offsetx = colcoords[0][1]
    offsety = rowcoords[0][1]
    half = 0 #(rowcoords[1][1]-rowcoords[0][1])/2.
    maxleny = rowcoords[-1][2]-rowcoords[0][1]+1-half
    maxlenx = colcoords[-1][2]-colcoords[0][1]+1


    for chrx, coords in regions:
        for start, end in coords:
            if (chrx, start) in colcoords:
                x_regions.append((chrx, start, end))
            elif (chrx, start) in rowcoords:
                y_regions.append((chrx, start, end))

    for _, start, end in x_regions:
        ax.add_patch(patches.Rectangle((start-offsetx, half), end-start,
            maxleny, edgecolor=HIGHLIGHT_COLOR, linewidth=HIGHLIGHT_WIDTH,
            fill=False))

    for _, start, end in y_regions:
        ax.add_patch(patches.Rectangle((0, maxleny-(end-offsety)),
            maxlenx, end-start, edgecolor=HIGHLIGHT_COLOR,
            linewidth=HIGHLIGHT_WIDTH, fill=False))


def doPlot(mtrx, colcoords, rowcoords, colannotations, rowannotations,
        color_min, color_max, outlier_cols, outlier_rows, axes_labels,
        highlight_regions=None):

    colorfunc = constructColorFunc(color_min, color_max,
            plt.get_cmap('viridis'))

    startx = colcoords[0][1]
    endx = colcoords[-1][2]
    starty = rowcoords[0][1]
    endy = rowcoords[-1][2]
    ax = plt.figure(1, figsize = (8, 8))

    mtrxSize = 0.75
    geneHeight = 0.15
    left_start = 0
    if axes_labels:
        left_start = 0.025
        mtrxSize -= left_start
        axYlabel = plt.axes([0, 0, left_start, mtrxSize], frameon=False)
        axYlabel.set_ylim([0, 1])
        axYlabel.set_xlim([0, 1])
        axYlabel.set_xticks([])
        axYlabel.set_yticks([])
        axXlabel = plt.axes([geneHeight+left_start, 0.975, mtrxSize,
            left_start], frameon=False)
        axXlabel.set_ylim([0, 1])
        axXlabel.set_xlim([0, 1])
        axXlabel.set_xticks([])
        axXlabel.set_yticks([])
        axXlabel.text(0.5, 0.2, axes_labels[0], fontsize=12, fontweight='bold',
                horizontalalignment='center')
        axYlabel.text(0.1, 0.5, axes_labels[1], fontsize=12, fontweight='bold',
                verticalalignment='center', rotation='vertical')

    # [left, bottom, width, height]
    axGenesX = plt.axes([geneHeight+left_start, mtrxSize, mtrxSize,
        1-mtrxSize-left_start], frameon=False)
    axGenesY = plt.axes([left_start, 0, geneHeight, mtrxSize], frameon=False)
    axMtrx = plt.axes([geneHeight+left_start, 0, mtrxSize, mtrxSize], frameon=False)
    # XXX disable for Figure 2 plots
    axColor = plt.axes([0.95, 0, 0.1, mtrxSize], frameon=False)
    axColor.imshow(plt.transpose(( plt.linspace(1, 0, 256),)),
            plt.get_cmap('viridis'),
        aspect='auto', extent=[0,1,color_min, color_max] )
    axColor.set_xticks([])
    # XXX end
    axMtrx.set_xticks([])
    axMtrx.set_yticks([])
    axMtrx.set_xlim([0, endx-startx+1])
    axMtrx.set_ylim([0, endy-starty+1])
    axGenesX.set_xticks([])
    axGenesX.set_yticks([])
    axGenesY.set_xticks([])
    axGenesY.set_yticks([])
    axGenesX.set_xlim([0, endx-startx+1])
    axGenesX.set_ylim([0, 1])
    axGenesY.set_xlim([0, 1.2])
    axGenesY.set_ylim([0, endy-starty+1])
    #plt.colorbar(im, cax=axMtrx)
    drawMatrix(axMtrx, startx, starty, endy-starty+1, colorfunc,
            mtrx, colcoords, rowcoords, outlier_cols, outlier_rows)
    axMtrx.set_xticks([])
    axMtrx.set_yticks([])


    drawCoordsX(axGenesX, startx, 0.025, startx, endx)
    drawGenesX(axGenesX, startx, 0.21, colannotations, ax, 12)
    drawCoordsY(axGenesY, 0.95, endy+1, starty, endy)
    drawGenesY(axGenesY, 0, starty, endy-starty+1, rowannotations, ax)
    #plt.axis('off')
    if highlight_regions:
        highlightRegions(axMtrx, highlight_regions, colcoords, rowcoords)


if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('annotation', type=file,
            help='Annotation file in GFF format')
    parser.add_argument('hic_map', type=file, help='(Ca-)HI-C map in COOL format')
    parser.add_argument('coord1', type=str,
            help='coordinates of submatrix that is subject to p-value calculation')
    parser.add_argument('coord2', type=str, nargs='?',
            help='coordinates of submatrix that is subject to p-value calculation')
    parser.add_argument('-l', '--labels', type=str, nargs=2, help='axes labels')
    parser.add_argument('-r', '--highlight_regions', type=str,
            help='mark specified regions by a rectangle')
    parser.add_argument('-a', '--minval', default=None, type=float,
            help='fixed minimum value for color coding; if set to None, ' + \
                    'the minimum value of the matrix is used')
    parser.add_argument('-z', '--maxval', default=None, type=float,
            help='fixed maximum value for color coding; if set to None, ' + \
                    'the maximum value of the matrix is used')
    parser.add_argument('-g', '--highlight_genes', type=str, nargs='*',
            help='do not label all genes in the region, but only those' + \
                    'specified')
    parser.add_argument('-i', '--ignore_zeros', action = 'store_true',
            help='do not plot cells with zero count')
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


    if 'outlier' in bins.columns:
        outlier1 = bins['outlier'].values[sel_coord1]
        outlier2 = bins['outlier'].values[sel_coord2]
    else:
        outlier1 = np.zeros(len(coord1), dtype=bool)
        outlier2 = np.zeros(len(coord2), dtype=bool)

    isSymmetric = np.all(sel_coord1 == sel_coord2)
    submtrx = mtrx[np.ix_(sel_coord2,sel_coord1)]

    if args.ignore_zeros:
        LOG.info('clear zeros from submatrix')
        submtrx[submtrx == 0] = np.nan

    minval, maxval = np.quantile(mtrx[np.isfinite(mtrx)], (0.05, 0.95))
    if args.minval != None:
        minval = args.minval
    if args.maxval != None:
        maxval = args.maxval

    # TODO change to single reading of annotation file
    colannotations = readGFF3Data(open(args.annotation.name), coord1[0][0],
            coord1[0][1], coord1[-1][2], ('gene', ), ('protein_coding', ))
    #, 'exon'))
    rowannotations = readGFF3Data(open(args.annotation.name), coord2[0][0],
            coord2[0][1], coord2[-1][2], ('gene', ), ('protein_coding', ))
    #, 'exon'))

    if args.highlight_genes:
        __labels__ = set(map(lambda x: x.upper(), args.highlight_genes))
        colannotations_new = list()
        for chrx, start, end, orient, atype, names in colannotations:
            s = [x for x in names if x.upper() in __labels__]
            colannotations_new.append((chrx, start, end, orient, atype, s and
                s[0] or ''))

        colannotations = colannotations_new
        rowannotations_new = list()
        for chrx, start, end, orient, atype, names in rowannotations:
            s = [x for x in names if x.upper() in __labels__]
            rowannotations_new.append((chrx, start, end, orient, atype, s and
                s[0] or ''))
        rowannotations = rowannotations_new
    else:
        for i, el in enumerate(colannotations):
            name = ''
            if len(colannotations) < 200:
                name = el[-1] > 1 and el[-1][1] or el[-1][0]
            colannotations[i] = el[:-1] + (name, )
        for i, el in enumerate(rowannotations):
            name = ''
            if len(rowannotations) < 200:
                name = el[-1] > 1 and el[-1][1] or el[-1][0]
            rowannotations[i] = el[:-1] + (name, )

    axes_labels = args.labels or ('', '')

    highlight_regions = None
    if args.highlight_regions:
        highlight_regions = parseCondensedCoords(args.highlight_regions)

#    colannotations = list()
#    rowannotations = list()

    doPlot(submtrx, coord1, coord2, colannotations, rowannotations, minval,
            maxval, outlier1, outlier2, axes_labels,
            highlight_regions=highlight_regions)
    #plt.show()

    plt.savefig(stdout, format='pdf')
