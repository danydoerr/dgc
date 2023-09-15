#!/usr/bin/env python3

from collections.abc import Mapping
from itertools import chain, product, combinations
from os.path import isfile, join
from sys import stdout, stderr, exit
import configparser as cp
import csv
import re

import cooler as clr

PAT_COORD = re.compile('^([^:]+):(\d+)-(\d+)')
REGIONS_OF_INTEREST_FILE            = 'roi.csv'
COMPARISONS_OF_INTEREST_FILE        = 'coi.csv'
DATASET_CONFIG_FILE                 = 'datset.cfg'
ANNOT_FILE_KEY                      = 'gcf_annotation_file'
EXPR_FILE_KEY                       = 'deseq2_expression_file'

ROI_HEADER = ('genomic region name', 'categories', 'genomic coordinates',
        'highlighted genes')
COI_HEADER = ('genomic region name/category 1',
        'genomic region name/category 2', 'pools')


class Manager(Mapping):

    def __init__(self,
            config_file             = DATASET_CONFIG_FILE,
            roi_file                = REGIONS_OF_INTEREST_FILE,
            coi_file                = COMPARISONS_OF_INTEREST_FILE,
            create_if_not_existing  = False):

        #
        # read config
        #

        self.config = ConfigParser(defaults={
            ANNOT_FILE_KEY: '',
            EXPR_FILE_KEY: ''}
            )
        self.config.readfp(open(config_file))
        self.expression_file = self.config.get(cp.DEFAULTSECT, EXPR_FILE_KEY)

        #
        # check interest files
        #
        for fname, header in ((roi_file, ROI_HEADER), (coi_file, COI_HEADER)):
            if not isfile(fname):
                if create_if_not_existing:
                    f = open(fname, 'w')
                    f.write('\t'.join(map(lambda x: x.capitalize(), header)))
                    f.write('\n')
                    f.close()
                    print('file %s not found, created new file at given path' %fname)
                else:
                    raise OSError('file %s not found -- wrong path?' %fname)

        #
        # set up poolerimental data structures
        #
        self.rois = self.parseRegionsOfInterest(open(roi_file))

        self.chromosomes = frozenset(map(lambda x: PAT_COORD.match(x[2]).group(
            1), self.rois))

        self.categories = frozenset(chain(*map(lambda x: x[1], self.rois)))
        self.category_members = self.createCategory2IdsMap(self.categories,
                self.rois)
        self.name2id = dict(map(lambda x: (x[1][0], x[0]), \
                enumerate(self.rois)))
        self.filename2id = dict(map(lambda x: (x[1][0].replace(' ', \
                '_').lower(), x[0]), enumerate(self.rois)))
        self.conditions = frozenset(self.config.sections())

        pools = frozenset((x[0] == x[1] and x[0] or '-'.join(x) for x in
            product(self.conditions, self.conditions)))

        cois = self.parseComparisonsOfInterest(open(coi_file), self.categories,
                pools)
        self.cois = dict(map(lambda x: (x, set()), pools))

        single_pools = set()
        paired_pools = set()
        for c1, c2, pools in cois:
            for pool in pools:
                self.cois[pool].add((c1, c2))
                if pool in self.conditions:
                    single_pools.add(pool)
                else:
                    paired_pools.add(pool)
        self.single_pools = frozenset(single_pools)
        self.paired_pools = frozenset(paired_pools)
        # make cois immutable
        for k, v in self.cois.items():
            self.cois[k] = frozenset(v)


    def parseRegionsOfInterest(self, data):
        res = list()

        isFirst = True
        for line in csv.reader(data, delimiter='\t'):
            if line[0].strip().startswith('#'):
                continue
            if isFirst:
                if tuple(map(lambda x: x.lower(), line))[:len(ROI_HEADER)] != \
                        ROI_HEADER:
                    raise AssertionError(('poolected header with column ' + \
                            'names "%s", but got "%s" instead') %(
                                ', '.join(ROI_HEADER), ', '.join(line)))
                isFirst = False
                continue

            name, categories, coords, genes = line[:4]

            coords = coords.strip()
            if not PAT_COORD.match(coords):
                raise AssertionError(('poolected coordinates in format %s, ' + \
                        'but got %s instead') %(PAT_COORD.pattern, coords))

            if not name.strip():
                name = coords

            categories = frozenset(chain(filter(None, map(lambda x:
                x.strip().lower(), categories.split(';'))), (name.lower(),
                    coords.lower())))

            genes = tuple(map(lambda x: x.strip(), genes.split(';')))
            res.append((name, categories, coords, genes))

        return res


    def parseComparisonsOfInterest(self, data, categories, pools):
        res = list()

        isFirst = True
        for line in csv.reader(data, delimiter='\t'):
            if line[0].strip().startswith('#'):
                continue
            if isFirst:
                if tuple(map(lambda x: x.lower(), line))[:len(COI_HEADER)] != \
                        COI_HEADER:
                    raise AssertionError(('expected header with column ' + \
                            'names "%s", but got "%s" instead') %(
                            ', '.join(COI_HEADER), ', '.join(line)))
                isFirst = False
                continue
            category1, category2, pools = line[:3]

            category1 = category1.strip().lower()
            category2 = category2.strip().lower()
            pools = frozenset(map(lambda x: x.strip().lower(), pools.split(';')))

            if pools.difference(pools):
                raise AssertionError('pool(s) %s is (are) unknown' %(
                    ', '.join(pools.difference(pools))))

            for c in (category1, category2):
                if c not in categories:
                    raise AssertionError('category %s is unknown' %c)
            res.append((category1, category2, pools))
        return res


    def createCategory2IdsMap(self, categories, rois):

        res = dict((c, set()) for c in categories)

        for i, (_, cats, _, _) in enumerate(rois):
            for cat in cats:
                res[cat].add(i)
        return res


    def samples(self, condition=None):
        assert condition == None or condition in self.conditions
        conditions = condition and (condition, ) or self.conditions
        for c in conditions:
            for replicate in self.config.options(c, no_defaults=True):
                if replicate != ANNOT_FILE_KEY:
                    yield (c, replicate)


    def originalFiles(self, condition=None):

        assert condition == None or condition in self.conditions

        conditions = condition and (condition, ) or self.conditions
        for c in conditions:
            for replicate in self.config.options(c, no_defaults=True):
                if replicate != ANNOT_FILE_KEY:
                    yield self.config.get(c, replicate)


    def originalFile(self, condition, replicate):

        assert self.config.has_option(condition, replicate, no_defaults=True)
        return self.config.get(condition, replicate)


    def comparisonFiles(self, pool=None, coi=None):

        # check and setup pool collection
        assert pool == None or pool in self.cois
        pools = pool != None and (pool, ) or self.cois.keys()

        # check comparison of interest paremeter
        assert coi == None or all(tuple(coi) in self.cois[p] for p in pools)

        for pool in pools:
            cois = coi != None and (coi, ) or self.cois[pool]
            for cat1, cat2 in cois:
                ids1 = self.category_members[cat1]
                ids2 = self.category_members[cat2]
                cat_dir = cat1 == cat2 and cat1 or '-'.join((cat1, cat2))
                if (len(ids1) ==  1 and cat1 == \
                        self.rois[next(iter(ids1))][0].lower()) and \
                        (len(ids2) == 1 and cat2 == \
                        self.rois[next(iter(ids2))][0].lower()):
                    cat_dir = ''
                # the comparison within a category is always symmetric
                if cat1 == cat2:
                    pairs = chain(combinations(ids1, 2), zip(ids1, ids1))
                else:
                    pairs = product(ids1, ids2)

                for id1, id2 in pairs:
                    fname = id1 == id2 and self.rois[id1][0] or \
                            '-'.join((self.rois[id1][0], self.rois[id2][0]))
                    yield join(pool, cat_dir, self.canonicalFname(fname))


    def canonicalFname(self, fname):
        return fname.replace(' ', '_').lower()


    def topDiffFiles(self, pool_strategy, topx, interaction_threshold,
            sub_dir='', zoom_category=None):

        zooms = zoom_category != None and map(lambda x: self.canonicalFname( \
                self.rois[x][0]), self.category_members.get(zoom_category, ())) \
                or ('', )

        for zoom in zooms:
            for x in self.paired_pools:
                c1, c2 = x.split('-', 1)
                for tx, it in product(topx, interaction_threshold):
                    yield join(x, sub_dir, 'top%s_peak%s%s' %(tx, it, zoom and
                        '_%s' %zoom), '%s:%s' %( c1, pool_strategy))
                    yield join(x, sub_dir, 'top%s_peak%s%s' %(tx, it, zoom and
                        '_%s' %zoom), '%s:%s' %( c2, pool_strategy))


    def isProperCategory(self, category):

        return category in self.category_members and \
                (self.rois[next(iter(self.category_members[category])) \
                ][0].lower() != category)


    def getResolution(self, cooler_path):

        c = clr.Cooler(cooler_path)
        bins = c.bins()[:]
        _, row = next(bins.iterrows())
        return row['end']-row['start']


    def getSingleColorRangeParameter(self):
        minval = maxval = ''
        if self.config.has_option(cp.DEFAULTSECT, 'single_color_min'):
            minval = '-a%s' %self.config.get(cp.DEFAULTSECT, 'single_color_min')
        if self.config.has_option(cp.DEFAULTSECT, 'single_color_max'):
            maxval = '-z%s' %self.config.get(cp.DEFAULTSECT, 'single_color_max')
        return ' '.join((minval, maxval))


    def getPairwiseColorRangeParameter(self):
        minval = maxval = ''
        if self.config.has_option(cp.DEFAULTSECT, 'pairwise_color_min'):
            minval = '-a%s' %self.config.get(cp.DEFAULTSECT, 'pairwise_color_min')
        if self.config.has_option(cp.DEFAULTSECT, 'pairwise_color_max'):
            maxval = '-z%s' %self.config.get(cp.DEFAULTSECT, 'pairwise_color_max')
        return ' '.join((minval, maxval))


    def getAnnotationFile(self, pool=None):
        if pool == None:
            return self.config.get(cp.DEFAULTSECT, ANNOT_FILE_KEY)
        return self.config.get(pool, ANNOT_FILE_KEY)


    def __getitem__(self, key):
        return self.cois.__getitem__(key)


    def __iter__(self):
        return self.cois.__iter__()


    def __len__(self):
        return self.cois.__len__()


class ConfigParser(cp.ConfigParser):
    """Can get options() without defaults
    """
    def options(self, section, no_defaults=False, **kwargs):
        if no_defaults:
            try:
                return list(self._sections[section].keys())
            except KeyError:
                raise NoSectionError(section)
        else:
            return super().options(section, **kwargs)


    def has_option(self, section, option, no_defaults=False, **kwargs):
        """Check for the existence of a given option in a given section."""
        if no_defaults:
            if not section or section == cp.DEFAULTSECT or section not in \
                    self._sections:
                return False
            else:
                option = self.optionxform(option)
                return option in self._sections[section]
        else:
            return super().has_option(section, option, **kwargs)


if __name__ == '__main__':
    print('THIS IS A MODULE')
