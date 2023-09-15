#!/usr/bin/env python3

from sys import stdout, stderr, exit
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF
from os.path import basename, join
import csv

from hic import parseCoords, parseCondensedCoords

def readStats(data):

    res = list()

    isHeader = True
    for line in csv.reader(data, delimiter='\t'):
        if isHeader:
            isHeader = False
            # rudimentary check whether table is in proper format
            assert len(line) == 7
            continue

        name = line[0]
        center_coords = parseCoords(line[1])
        percent_contacts = float(line[5])
        interact_regions = parseCondensedCoords(line[6])
        res.append((name, center_coords, percent_contacts, interact_regions))

    return res


def regions2bed(name, description, chrx, browserPos, regions, score, out):

    if chrx.isdigit():
        chrx = 'chr%s' %chrx

    print('browser position %s:%s-%s' %((chrx, ) + tuple(browserPos)), file=out)
    print('browser hide all', file=out)
    print('track name="%s" description="%s" visibility=1 useScore=1' %(name, 
            description), file=out)
    for (start, end) in regions:
        print('\t'.join((chrx, str(start), str(end), '%s-%s' %(start,
            end), str(int(score)))),file=out)


def __intersect__(int1, int2):
    
    start1, end1 = int1
    start2, end2 = int2

    if end1 == 'end':
        end1 = float('inf')
    if end2 == 'end':
        end2 = float('inf')

    return (start1 <= start2 and end1 >= start2) or (start2 <= start1 and
            end2 >= start1)



if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('stat_table', type=file, 
            help='file containng summary statistic of interaction regions')
    parser.add_argument('out_dir', type=str, 
            help='output directory')
    args = parser.parse_args()

    fileDesc = basename(args.stat_table.name).rsplit('.', 1)[0].replace('_', 
            ' ')

    stats = readStats(args.stat_table)

    for name, center_coords, percent_contacts, interact_regions in stats:
        bedName = '%s - %s:%s-%s' %((name, ) + tuple(center_coords))
        #bedName = name
        bedDescription = ' '.join((name, fileDesc))
        
        _, regions = [x for x in interact_regions if x[0] ==
                center_coords[0]][0]
        regions = [r  for r in regions if __intersect__(r, center_coords[1:])]
        fName = '%s.bed' %('_'.join((name, '%s:%s-%s'%tuple(center_coords))))

        out = open(join(args.out_dir, fName), 'w')
        regions2bed(bedName, bedDescription, center_coords[0],
                center_coords[1:], regions, percent_contacts * 1000, out)
        out.close()

        
