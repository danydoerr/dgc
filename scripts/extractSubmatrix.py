#!/usr/bin/env python

from sys import stdout, stderr, exit
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF
from os.path import basename

from hic import PAT_COORD, extractRegion, writeMtrx


if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('-l', '--axes_labels', type=str, nargs='+',
            help='Axes labels: if one is specified, both axes are labeled ' + \
                    'equally, if two are specified, the first corresponds ' + \
                    'to the x-axis, the second to the y-axis')
    parser.add_argument('matrix', type=str, help='Matrix in TRV format')
    parser.add_argument('coordinate1', type=str, 
            help='Genome coordinates of the submatrix, given in the ' + \
                    'format \'<chrid>:<start>-<end>\'. ')
    parser.add_argument('coordinate2', nargs='?', type=str, default=None,
            help='Genome coordinates of the submatrix, given in the ' + \
                    'format \'<chrid>:<start>-<end>\'. If only one ' + \
                    'coordinate is given, the matrix is assumed to be ' + \
                    'the symmetric matrix of this region.')
    args = parser.parse_args()

    m = PAT_COORD.match(args.coordinate1)

    if not m:
        print('Unable to parse genome coordinates %s. Exiting.' %args.coordinate1, file=stderr)
        exit(1)

    chr1, start1, end1 = m.groups()
    chr2, start2, end2 = chr1, start1, end1

    if args.coordinate2:
        m = PAT_COORD.match(args.coordinate2)
        if not m:
            print('Unable to parse genome coordinates %s. Exiting.' %args.coordinate2, file=stderr)
            exit(1)
        chr2, start2, end2 = m.groups()
     
    collabels, rowlabels, extractedMtrx = extractRegion(open(args.matrix),
            chr1, int(start1), int(end1), chr2, int(start2), int(end2))

    xlabel = basename(args.matrix)
    ylabel = ''

    if args.axes_labels:
        xlabel = ylabel = args.axes_labels[0]
        if len(args.axes_labels) > 1:
            ylabel = args.axes_labels[1]

    writeMtrx(extractedMtrx, collabels, rowlabels, stdout, xlabel, ylabel)

