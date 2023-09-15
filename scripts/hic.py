#!/usr/bin/env python2

from sys import stdout, stderr, exit, maxint
from bisect import bisect_left
from itertools import izip, compress
import csv
import re

from scipy.misc import comb
import numpy as np

PAT_LABEL = re.compile('^(.+)-(\d+)$')
PAT_COORD = re.compile('^([^:]+):(\d+)-(\d+)')


def equiSegmentRegions(rMap, bin_size, overlap = 0, ignore_regions=None):
    """ Extracts the larger perimeters of the regions that are captured by the
    set of reads that are specified in the rMap. Allows for small gaps within
    these regions. The overlap parameter allows to asign restriction fragments
    to nearby (non-overlapping) equi-sized segments that are within the
    perimeter of the overlap. If the overlap is negative, then each equi-sized
    segment must intersect with a restriction fragment by at least the overlap
    """
   
    assert bin_size > 0 #and overlap >= 0
    assert all(rMap[i][:3] < rMap[i+1][:3] for i in xrange(len(rMap)-1))

    # fill with dummy as first element
    res = [(None, 0, 0)]
    for i, row in enumerate(rMap):
        if ignore_regions is not None and not ignore_regions[i]:
            continue
        chrx, start, end = row[:3]
        boundary= (int(max(0, start-overlap))/bin_size) * bin_size
        overlap_relapse = min(int(start+overlap-1)/bin_size - int(max(0, \
                start-overlap))/bin_size + 1, len(res)) 
        while overlap_relapse:
            if res[-overlap_relapse][0] == chrx and res[-overlap_relapse][1] \
                    == boundary:
                res[-overlap_relapse][3][1] = i
                boundary += bin_size
            overlap_relapse -= 1
        while end + overlap - boundary > 0:
            res.append((chrx, boundary, boundary+bin_size-1, [i, i]))
            boundary += bin_size

    # boundaries may go beyond chromosome sizes, therefore, the list is iterated
    # again, removing those segments that start after chromosome ends
    #
    # at the same time, test whether the segmentation is sorted to ensure that
    # the result is correct
    res = res[1:]
    chr_ends = [(row[0], row[2]) for i, row in enumerate(rMap) if i ==
            len(rMap)-1 or rMap[i+1][0] != row[0]]
    i = j = 0
    while j < len(res):
        if j and res[j-1][:3] > res[j][:3]:
            raise Exception, ('Segmentation is not sorted, %s is before ' + \
                    '%s') %( str(res[j-1]), str(res[j]))
        # if chromosome has switched in res, advance in chr_ends
        if chr_ends[i][0] != res[j][0]:
            i += 1
        # either delete entry if it exceeds chromosome end, or advance in res
        if res[j][1] > chr_ends[i][1]:
            del res[j]
        else:
            j += 1

    return res




def label2coords(label):
    m = PAT_LABEL.match(label)
    if not m:
        raise Exception, ('Unknown label format "%s", expected ' + \
                'string matching pattern %s') %(label, PAT_LABEL.pattern)
    chrx, p = m.groups()
    return (chrx, int(p))
   

def parseCoords(coords):
    m = PAT_COORD.match(coords)
    if not m:
        raise Exception, ('Unknown coordinate format "%s", expected ' + \
                'string matching pattern %s') %(coords, PAT_COORD.pattern)
    chrx, start, end = m.groups()
    
    return (chrx, int(start), int(end))


def condenseCoords(coords,labels, l2i=None):

    assert all(coords[i] <=coords[i+1] for i in xrange(len(coords)-1))

    if l2i == None:
        l2i = dict(zip(labels, xrange(len(labels))))

    res = '%s:%s' %coords[0][:2]
    for i in xrange(1, len(coords)):
        next_label = labels[l2i[coords[i-1]]+1]
        if coords[i][0] != coords[i-1][0]:
            if next_label[0] == coords[i-1][0]:
                res += '-%s' %next_label[1]
            else:
                res += '-end'
            res += ';%s:%s' %coords[i][:2]
        elif coords[i] != next_label:
            res += '-%s,%s' %(next_label[1], coords[i][1])
    if l2i[coords[-1]]+1 >= len(labels) or labels[l2i[coords[-1]]+1][0] != coords[-1][0]:
        res += '-end'
    else:
        res += '-%s' %labels[l2i[coords[-1]]+1][1]

    return res


def parseCondensedCoords(condensed):

    res = list()

    for x in condensed.split(';'):
        chrx, coords = x.rsplit(':', 1)
        res.append((chrx, list()))
        for item in coords.split(','):
            start, end = item.split('-')
            start = int(start)
            if end != 'end':
                end = int(end)
            res[-1][1].append((start, end))

    return res


def writeMtrx(mtrx, colnames, rownames, out, xlabel=None, ylabel=None):
    print('%s\t%s\t%s' %(xlabel, ylabel, '\t'.join(colnames)), file=out)
    for i in xrange(mtrx.shape[0]):
        out.write('\t'.join((rownames[i], rownames[i])))
        for j in xrange(mtrx.shape[1]):
            out.write('\t')
            if mtrx[i][j] != None and mtrx[i][j] != 0:
                out.write(str(mtrx[i][j]))
        out.write('\n')

    i = len(mtrx)
    while i < len(rownames):
        print('\t'.join((rownames[i], rownames[i])), file=out)
        i += 1

def extractRegion(data, chr1, start1, end1, chr2, start2, end2):
    
    res = list()
    bounds = None 
    collabels = None
    rowlabels = None
    prev_line = None
    for c, line in enumerate(csv.reader(data, delimiter='\t')):
        if not c:
            coords = map(label2coords, line[2:])
            i = -1
            while i+1 < len(coords) and coords[i+1] <= (chr1, start1):
                i += 1
            j = i
            while j < len(coords) and coords[j] <= (chr1, end1):
                j += 1
            collabels = line[i+2:j+3]
            bounds = (i+2, j+2)
            continue
        coord = label2coords(line[0])
        if coord > (chr2, end2):
            if prev_line:
                rowlabels = [prev_line[0]]
                res.append(prev_line[bounds[0]:bounds[1]])
                prev_line = None
            rowlabels.append(line[0])
            break
        elif coord > (chr2, start2):
            if prev_line:
                rowlabels = [prev_line[0]]
                res.append(prev_line[bounds[0]:bounds[1]])
                prev_line = None
            rowlabels.append(line[0])
            res.append(line[bounds[0]:bounds[1]])
        else:
            prev_line = line
    
    return collabels, rowlabels, np.array(res)


def readRMap(data):

    res = list()

    for i, line in enumerate(csv.reader(data, delimiter='\t')):
        if len(line) < 4:
            print('skipping line no. %s because it has less ' + \
                    'than 4 columns') % i, file=stderr)
            continue
        chrx, start, end, numericID = line[:4]
        res.append((chrx, int(start), int(end), int(numericID)))

    res.sort()
    return res


def readBaitmapIDs(data):

    res = list()
    for line in csv.reader(data, delimiter='\t'):
        res.append(int(line[3]))
    return res


def readHiCMapTRV(data, regions=None, ignoreCols=set(), ignoreRows=set()):
   
    res = list()
    colnames = None
    rownames = None
    axes_labels = None

    select_indices = None
    for c, line in enumerate(csv.reader(data, delimiter='\t')):
        if not c:
            axes_labels = line[:2] 
            if regions != None:
                select_indices = [i for i, x in enumerate(line[2:]) if x in
                        regions and label2coords(x) not in ignoreCols]
            else:
                select_indices = [i for i, x in enumerate(line[2:]) if
                        label2coords(x) not in ignoreCols]
            # initialize col/rownames 
            colnames = line[2:]
            rownames = list()
        else:
            row = [None] * (len(line)-2)
            res.append(row)
            rownames.append(line[0])
            if label2coords(line[0]) in ignoreRows:
                continue
            else:
                sj = 0
                for j, x in enumerate(line[2:]):
                    while sj < len(select_indices) and select_indices[sj] < j:
                        sj += 1
                    val = None
                    if sj < len(select_indices) and j == select_indices[sj] \
                            and x:
                        val = float(x)
                    row[j] = val
    if res and not res[-1]:
        del res[-1]

    return np.array(res), colnames, rownames, axes_labels


def readHiCMapTRVAsPdist(data, max_dist, ignore_cols=set(), ignore_rows=set()):
    """ returns an numpy object representing the upper triangle matrix in the
    format described by pdist:
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html

        NOTE: values are None if cell is empty or cells are further away than
        max_dist (determined by column/row labeling)
    """
  
    def __f__(x):
        a, b = x.rsplit('-', 1)
        return (a, int(b))
    
    __d__ = lambda x, y: x[0] != y[0] and np.inf or y[1]-x[1]

    res = None
    ycoords = None
    axes_labels = None
    mx = -np.inf
    i = 0
    for c, line in enumerate(csv.reader(data, delimiter='\t')):
        if not c:
            axes_labels = line[:2]
            res = np.empty(int(comb(len(line)-2, 2)), dtype=object)
            ycoords = map(__f__, line[2:])
        else:
            if label2coords(line[0]) not in ignore_rows:
                p = c
                while p < len(ycoords) and __d__(ycoords[c-1], ycoords[p]) \
                        <= max_dist:
                    if line[p+2] != '' and ycoords[p] not in ignore_cols:
                        val = float(line[p+2])
                        mx = max(val, mx)
                        res[i+p-c] = val 
                    p += 1
            i += len(line)-c-2
        
    return res, ycoords, mx, axes_labels


def pDistIndex(n, i, j):
    assert i <= j

    return int(comb(n, 2)-comb(n-i, 2)+j-i-1)


def triu_condensed(mtrx, k=0):

    assert len(mtrx.shape) == 2 and mtrx.shape[0] == mtrx.shape[1]

    res = np.empty(int(comb(mtrx.shape[0]+1-k, 2)), dtype=mtrx.dtype)
    c = 0
    for i in xrange(mtrx.shape[0]-k):
        res[c:c+mtrx.shape[0]-k-i] = mtrx[i,i+k:]
        c += mtrx.shape[0]-k-i
    return res


def getDiagonalValues(mtrx, select_indices_col, select_indices_row,
        diag_offset, ignore_cols=None, ignore_rows = None, excludeNones=False):
  
    pool_baited = list()
    pool_others = list()

    si = 0
    start = max(0, -diag_offset)
    # min over: 
    # (1) boundary of rectengular matrix if diag_offset is negative and
    #   |columns| < |rows|
    # (1) boundary of rectengular matrix if diag_offset is positive and 
    #   |columns| < |rows|
    # (2) boundary if diag_offset is negative (then it will go to the last row)
    # (3) boundary if diag_offset is positive (then it will hit the column
    #   border at some earlier row)
    end = min(start+mtrx.shape[1], start+mtrx.shape[1]-diag_offset,
            mtrx.shape[0], mtrx.shape[0]-diag_offset)
    for i in xrange(start, end):
        if ignore_rows and i in ignore_rows:
            continue

        j = i+diag_offset
        if ignore_cols and j in ignore_cols:
            continue
        while si < len(select_indices_row) and select_indices_row[si] < i:
            si += 1
        sj = bisect_left(select_indices_col, j)

        i_selected = si < len(select_indices_row) and i == \
                select_indices_row[si]
        j_selected = sj < len(select_indices_col) and j == \
                select_indices_col[sj]

        if i_selected and j_selected:
            if not excludeNones or mtrx[i,j] != None:
                pool_baited.append(mtrx[i,j])
        elif i_selected or j_selected:
            if not excludeNones or mtrx[i,j] != None:
                pool_others.append(mtrx[i,j])
             
    return pool_baited, pool_others


def removeDisjointColumns(mtrx, colnames1, colnames2):
    c2j = dict(izip(colnames1, xrange(len(colnames1))))
    delc = set(colnames1).difference(colnames2)
    keep = [True] * len(colnames1)

    for c, j in enumerate(sorted(map(c2j.get, delc))):
        del colnames1[j-c]
        keep[j] = False

    return mtrx[:,keep]

def removeDisjointRows(mtrx, rownames1, rownames2):
    r2i = dict(izip(rownames1, xrange(len(rownames1))))
    delr = set(rownames1).difference(rownames2)
   
    keep = [True] * len(rownames1)
    for c, i in enumerate(sorted(map(r2i.get, delr))):
        del rownames1[i-c]
        keep[i] = False

    return mtrx[keep, :]


def assimilateMatrices(mtrxs, return_common=False):

    mtrxs = list(mtrxs)

    assert len(set(map(len, mtrxs))) == 1 and len(mtrxs[0]) >= 3

    res = list()

    colnames = None
    rownames = None

    for m in mtrxs: 
        mtrx, fcolnames, frownames = m[:3]
        if colnames == None:
            colnames = set(fcolnames)
        else:
            colnames.intersection_update(fcolnames)

        if rownames == None:
            rownames = set(frownames)
        else:
            rownames.intersection_update(frownames)
    common = list()
    for m in mtrxs:
        mtrx, fcolnames, frownames = m[:3]
        sel_col = np.array([x in colnames for x in fcolnames])
        sel_row = np.array([x in rownames for x in frownames])
        common.append((sel_col, sel_row))
        res.append((mtrx[np.ix_(sel_col, sel_row)], list(compress(fcolnames,
            sel_col)), list(compress(frownames, sel_row))) + tuple(m[3:]))

    if return_common:
        return res, common
    return res


if __name__ == '__main__':
    print('this is a module')
