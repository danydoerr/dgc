#!/usr/bin/env python3

from sys import stdout, stderr, exit
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF, \
        FileType
from itertools import chain
from os.path import basename
import logging
import csv

import os
if not os.environ.get('DISPLAY', None):
    import matplotlib; matplotlib.use('Agg')

import numpy as np
from matplotlib import pylab as plt, patches, lines
import matplotlib

matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.unicode'] = True
matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']

LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)

GREEN = 'C1' #(0.2,0.91,0.39)
RED = 'C0' #(1.0,0.27,0.14)

def readStats(data, ttype):

    for line in csv.reader(data, delimiter='\t'):
        if line[0].startswith(ttype):
            return float(line[3]), line[2] == 'significant' and 1 or 0, \
                    float(line[1])


def plotSumCount(data, labels, out, ylim=None):

    f = plt.figure()
    p = plt.bar(xrange(data.shape[0]), data[:, 0], tick_label=labels, color =
            RED)

    for i, r in enumerate(p.patches):
        if type(r) != patches.Rectangle:
            continue
        height = r.get_height()
        va = 'bottom'
        if height < 0:
            va= 'top'
        v = data[i,2]
        b = np.floor(np.log10(v))
        plt.text(r.get_x()+r.get_width()/2.0, height, (np.isinf(b) or \
                np.isnan(b)) and r'$0$' or (r'$%.2f%s$' %( v * 10**-b, int(b) \
                and r'\cdot 10^{%.0f}'%b or '')), ha='center', va=va, \
                fontsize=12, color= data[i, 1] and GREEN or 'gray')

    plt.tick_params(axis='x', which='major', labelsize=14)
    plt.xticks(rotation=90)

    yrange = np.concatenate((data[:, 0], np.array((0, ))))
    if ylim:
        yrange = np.concatenate((yrange, np.array(ylim)))
    
    ticks = np.linspace(np.min(yrange), np.max(yrange), 5)
    #plt.yscale('symlog')
    plt.yticks(ticks)
    plt.axhline(0, color='k')
    plt.tight_layout(pad=2)
    plt.ylabel('sum-count', fontsize=14)
    plt.xlabel('gene clusters', fontsize=14)

    n = lines.Line2D([], [], color='gray', label=r'\textbf{$p$-value}')
    if args.type == 'low' or args.type == 'high': 
        s = lines.Line2D([], [], color=GREEN, label=
                r'\textbf{significant ($\alpha_{\textbf{corr}} \geq 7.69\cdot 10^{-5}$)}')
    else:
        s = lines.Line2D([], [], color=GREEN, label=
                r'\textbf{significant ($\alpha = 0.001$)}')

    plt.legend(handles=[n, s], loc = 'lower right')

    plt.savefig(out, format='pdf')


if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('stats', nargs='+', type=file, 
            help='statistics files')
    parser.add_argument('-t', '--type', type=str, choices=('full', 'high', 'low'),
            default='full', help='statistics files')
    parser.add_argument('-y', '--ylim', type=float, nargs=2,
            help='limits of y-axis (unless they are exceeded by data)')
    args = parser.parse_args()

    # setup logging
    ch = logging.StreamHandler(stderr)
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(logging.Formatter('%(levelname)s\t%(asctime)s\t%(message)s'))
    LOG.addHandler(ch)

    # two columns: sumcount, signficant 0/1
    data = np.zeros((len(args.stats), 3))
    labels = list()

    for i, f in enumerate(args.stats):
        data[i,:] = readStats(f, args.type)
        labels.append(basename(f.name).rsplit('.', 1)[0].replace('_', ' '))

    if args.ylim:
        plotSumCount(data, labels, stdout, ylim=args.ylim)
    else:
        plotSumCount(data, labels, stdout)

