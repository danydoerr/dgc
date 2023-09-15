configfile: 'config.yaml'

from itertools import chain, repeat
from importlib import util as imputil
from os.path import basename

# where is the software directory?
SCRIPT_DIR      = config.get('script_dir', 'scripts')
# load experiment manager
spec = imputil.spec_from_file_location('experiment', '%s/experiment.py' % \
        SCRIPT_DIR)
exp = imputil.module_from_spec(spec)
spec.loader.exec_module(exp)

#
# read configuration and setup local variables
#

BEDS_DIR        = config.get('beds_dir', 'beds')
CLUSTER_TYPES   = config.get('cluster_types', ['single'])
COI_FILE        = config.get('coi_file', exp.COMPARISONS_OF_INTEREST_FILE)
CIRCOS_BIN      = config.get('circos_bin', 'circos')
DATASET_FILE    = config.get('dataset_file', exp.DATASET_CONFIG_FILE)
JUICER_TOOLS    = config.get('juicer_tools_path', '/path/to/juicer_tools.jar')
ORIG_MAPS_DIR   = config.get('orig_maps_dir', 'original')
PLOT_DIR        = config.get('submtrx_plots_dir', 'plots')
POOLED_MAPS_DIR = config.get('pooled_maps_dir', 'pooled')
POOL_STRAT_CON  = config.get('pooling_strategy_conditions', 'normavg')
POOL_STRAT_REP  = config.get('pooling_strategy_replicates', 'min')
PVALUE_ALPHA    = config.get('significance_threshold_calc_pval', 0.05)
ROI_FILE        = config.get('roi_file', exp.REGIONS_OF_INTEREST_FILE)
SAMPLE_SIZE     = config.get('empirical_dist_samples', 1000)
STATS_DIR       = config.get('stats_dir', 'stats')
TOP_DIFF_TOPX   = config.get('top_diff_plot_topx', [200])
TOP_DIFF_IT     = config.get('top_diff_plot_interation_threshold', [0.01])

#
# setup experiment manager
#
MANAGER = exp.Manager(config_file=DATASET_FILE, roi_file=ROI_FILE,
        coi_file=COI_FILE)


rule all:
    input:
#        expand('%s/{condition}:%s.cool' %(POOLED_MAPS_DIR, POOL_STRAT_REP),
#                condition=MANAGER.single_pools),
#        expand('%s/{condition_pair}:%s.cool' %(POOLED_MAPS_DIR, POOL_STRAT_CON),
#                condition_pair=MANAGER.paired_pools),
#        expand('%s/{out}.txt' %STATS_DIR,
#                out=MANAGER.comparisonFiles()),
        expand('%s/{out}.pdf' %PLOT_DIR,
                out=MANAGER.comparisonFiles()),
#        expand('%s/{out}.csv' %STATS_DIR, out=chain(*(map(lambda x:
#                '/'.join((x[0], x[1][0] == x[1][1] and x[1][0] or
#                '-'.join(x[1]))), zip(repeat(p), filter(lambda x:
#                all(map(MANAGER.isProperCategory, x)), MANAGER.cois[p]))) for p
#                in chain(MANAGER.single_pools, MANAGER.paired_pools)))),
#        expand('%s/{out}.shapiro.csv' %STATS_DIR, out=chain(*(map(lambda x:
#                '/'.join((x[0], x[1][0] == x[1][1] and x[1][0] or
#                '-'.join(x[1]))), zip(repeat(p), filter(lambda x:
#                all(map(MANAGER.isProperCategory, x)), MANAGER.cois[p]))) for p
#                in MANAGER.paired_pools))),
#        expand('%s/{file}.png' %PLOT_DIR, file=MANAGER.topDiffFiles(
#                POOL_STRAT_REP, TOP_DIFF_TOPX, TOP_DIFF_IT)),
#        expand('%s/{file}.png' %PLOT_DIR, file=MANAGER.topDiffFiles(
#                POOL_STRAT_REP, TOP_DIFF_TOPX, TOP_DIFF_IT,
#                zoom_category='genecluster'))



rule link_hic_maps:
    input:
        hic_maps = MANAGER.originalFiles()
    output:
        expand('%s/{hic_map}.cool' %ORIG_MAPS_DIR, hic_map=map(lambda x:
                ':'.join(x), MANAGER.samples()))
    run:
        from os import symlink, path
        for (c, r) in MANAGER.samples():
            src = MANAGER.originalFile(c, r)
            symlink(path.relpath(src, ORIG_MAPS_DIR), path.join(ORIG_MAPS_DIR,
                    '%s:%s.cool' %(c, r)))

rule cooler_to_juicer:
    input:
        maps = '{subdir}/{hic_map}.cool'
    output:
        sf = temp('{subdir}/{hic_map,[^/]+}.sf'),
        chr_size = temp('{subdir}/{hic_map}.chromsizes')
    log:
        'logs/{subdir}/cool2juicer_{hic_map}.log'
    shell:
        '%s/cool2juicer.py -s {output.chr_size} {input} ' %SCRIPT_DIR +
        '> {output.sf} 2> {log}'


rule juicer_to_hic:
    input:
        sf = '{subdir}/{hic_map}.sf',
        chr_size = '{subdir}/{hic_map}.chromsizes'
    params:
        resolution = lambda wildcards: MANAGER.getResolution('%s/%s.cool' %(
                wildcards.subdir, wildcards.hic_map))
    output:
        '{subdir}/{hic_map,[^/]+}.hic'
    log:
        'logs/{subdir}/juicer_tools_pre_{hic_map}.log'
    shell:
        'java -jar %s pre -n -r {params.resolution} {input.sf} ' % JUICER_TOOLS +
        '{output} {input.chr_size} > {log}'


rule combine_conditions_single:
    input:
        lambda wildcards: expand('%s/%s:{rep}.cool' %(ORIG_MAPS_DIR,
                wildcards.condition), rep=map(lambda x: x[1],
                MANAGER.samples(wildcards.condition))),
    output:
        hic_map = '%s/{condition,[^:-]+}:{opt,[^:]+}.cool' %POOLED_MAPS_DIR,
    log:
        'logs/combineMtrx_{condition}_{opt}.log'
    shell:
        '%s/combineMtrx.py -1 {input} -- {wildcards.opt} ' %SCRIPT_DIR +
        '{output} 2> {log}'

rule combine_conditions_paired:
    input:
        condition1 = lambda wildcards: expand('%s/%s:{rep}.cool' %(
                ORIG_MAPS_DIR, wildcards.condition1), rep=map(lambda x: x[1],
                MANAGER.samples(wildcards.condition1))),
        condition2 = lambda wildcards: expand('%s/%s:{rep}.cool' %(
                ORIG_MAPS_DIR, wildcards.condition2), rep=map(lambda x: x[1],
                MANAGER.samples(wildcards.condition2))),
    output:
        '%s/{condition1,[^:-]+}-{condition2,[^:-]+}:{opt,[^:]+}.cool' % \
                POOLED_MAPS_DIR,
    log:
        'logs/combineMtrx_{condition1}_{condition2}_{opt}.log'
    shell:
        '%s/combineMtrx.py -1 {input.condition1} -2 ' %SCRIPT_DIR +
        '{input.condition2} -- {wildcards.opt} {output} 2> {log}'


rule plot_submatrix_single:
    input:
        hic_map = lambda wildcards: '%s/%s:%s.cool' %(POOLED_MAPS_DIR,
                wildcards.hic_map, wildcards.hic_map in MANAGER.single_pools and
                POOL_STRAT_REP or POOL_STRAT_CON),
        annotation = lambda wildcards: MANAGER.getAnnotationFile(
                wildcards.hic_map)
    params:
        coords = lambda wildcards: MANAGER.rois[MANAGER.filename2id[ \
                wildcards.gc]][2],
        highlight = lambda wildcards: MANAGER.rois[MANAGER.filename2id[ \
                wildcards.gc]][3],
        label = lambda wildcards: MANAGER.rois[MANAGER.filename2id[\
                wildcards.gc]][0],
        ignore_zeros = lambda wildcards: wildcards.hic_map in
                MANAGER.single_pools and '-i' or '',
        color_range = MANAGER.getSingleColorRangeParameter()
    output:
        '%s/{hic_map,[^/]+}{category,/?[^/]+/}{gc,[^/-]+}.pdf' %PLOT_DIR,
    log:
        'logs{category}plotSubmatrix_{hic_map}_{gc}.log'
    shell:
        '%s/plotSubmatrix.py -l \'{params.label}\' ' %SCRIPT_DIR +
        '\'{params.label}\' {params.color_range} {params.ignore_zeros} '
        '-g {params.highlight} -- {input.annotation} {input.hic_map} '
        '"{params.coords}" > {output} 2> {log}'


rule plot_submatrix_pairwise:
    input:
        hic_map = lambda wildcards: '%s/%s:%s.cool' %(POOLED_MAPS_DIR,
                wildcards.hic_map, wildcards.hic_map in MANAGER.single_pools and
                POOL_STRAT_REP or POOL_STRAT_CON),
        annotation = MANAGER.getAnnotationFile(),
    params:
        coords1 = lambda wildcards: MANAGER.rois[MANAGER.filename2id[\
                wildcards.gc1]][2],
        coords2 = lambda wildcards: MANAGER.rois[MANAGER.filename2id[\
                wildcards.gc2]][2],
        highlight1 = lambda wildcards: MANAGER.rois[MANAGER.filename2id[\
                wildcards.gc1]][3],
        highlight2 = lambda wildcards: MANAGER.rois[MANAGER.filename2id[\
                wildcards.gc2]][3],
        label1 = lambda wildcards: MANAGER.rois[MANAGER.filename2id[\
                wildcards.gc1]][0],
        label2 = lambda wildcards: MANAGER.rois[MANAGER.filename2id[\
                wildcards.gc2]][0],
        ignore_zeros = lambda wildcards: wildcards.hic_map in
                MANAGER.single_pools and '-i' or '',
        color_range = MANAGER.getPairwiseColorRangeParameter()
    output:
        '%s/{hic_map,[^/]+}{category,.*/}{gc1,[^./-]+}-{gc2,[^./-]+}.pdf' %PLOT_DIR,
    log:
        'logs{category}plotSubmatrix_{hic_map}_{gc1}-{gc2}.log'
    shell:
        '%s/plotSubmatrix.py -l \'{params.label1}\' ' %SCRIPT_DIR +
        '\'{params.label2}\' {params.color_range} {params.ignore_zeros} '
        '-g {params.highlight1} {params.highlight2} -- {input.annotation} '
        '{input.hic_map} "{params.coords1}" "{params.coords2}" > {output} '
        '2> {log}'


#rule plot_centromeric_contact_counts_vs_all:
#    input:
#        hic_maps = expand('%s/{rep}:{comb}.trv' %ORIG_MAPS_DIR, rep=REPLICATES,
#                comb=DATASET.sections()),
#        outliers = OUTLIER_FILE
#    output:
#        'contact_dist_{centromere,[^:+]:\d+-\d+}.pdf'
#    shell:
#        '%s/plotContactVsAll.py -i {input.outliers} ' %SCRIPT_DIR +
#        '{input.hic_maps} {wildcards.centromere} > {output}'


rule calc_pvalues_single:
    input:
        hic_map = lambda wildcards: '%s/%s:%s.cool' %(
                POOLED_MAPS_DIR, wildcards.hic_map, wildcards.hic_map in
                MANAGER.single_pools and POOL_STRAT_REP or POOL_STRAT_CON)
    params:
        repeats = SAMPLE_SIZE,
        alpha = PVALUE_ALPHA,
        coords = lambda wildcards: MANAGER.rois[MANAGER.filename2id[\
                wildcards.gc]][2],
    output:
        pval = '%s/{hic_map,[^/]+}{category,.*/}{gc,[^./-]+}.txt' %STATS_DIR,
        plot = '%s/{hic_map}{category}{gc}.pdf' %STATS_DIR
    log:
        'logs{category}calculatePvalue_{hic_map}_{gc}.log'
    shell:
        '%s/calculatePvalue.py -r {params.repeats} ' %SCRIPT_DIR +
        '-a {params.alpha} -p {output.plot} {input.hic_map} {params.coords} '
        '> {output.pval} 2> {log} || true'


rule calc_pvalues_pairwise:
    input:
        hic_map = lambda wildcards: '%s/%s:%s.cool' %(
                POOLED_MAPS_DIR, wildcards.hic_map, wildcards.hic_map in
                MANAGER.single_pools and POOL_STRAT_REP or POOL_STRAT_CON)
    params:
        repeats = SAMPLE_SIZE,
        alpha = PVALUE_ALPHA,
        coords1 = lambda wildcards: MANAGER.rois[MANAGER.filename2id[\
                wildcards.gc1]][2],
        coords2 = lambda wildcards: MANAGER.rois[MANAGER.filename2id[\
                wildcards.gc2]][2],
    output:
        pval = '%s/{hic_map,[^/]+}{category,.*/}{gc1,[^./-]+}-{gc2,[^./-]+}.txt' %STATS_DIR,
        plot = '%s/{hic_map}{category}{gc1}-{gc2}.pdf' %STATS_DIR
    log:
        'logs{category}calculatePvalue_{hic_map}_{gc1}-{gc2}.log'
    shell:
        '%s/calculatePvalue.py -r {params.repeats} ' %SCRIPT_DIR +
        '-a {params.alpha} -p {output.plot} {input.hic_map} {params.coords1} '
        '{params.coords2} > {output.pval} 2> {log} || true'


rule calc_pvalues_shapiro_single:
    input:
        hic_map = lambda wildcards: '%s/%s:%s.cool' %(
                POOLED_MAPS_DIR, wildcards.hic_map, wildcards.hic_map in
                MANAGER.single_pools and POOL_STRAT_REP or POOL_STRAT_CON)
    params:
        repeats = SAMPLE_SIZE,
        alpha = PVALUE_ALPHA,
        coords = lambda wildcards: MANAGER.rois[MANAGER.filename2id[\
                wildcards.gc]][2],
    output:
        pval = '%s/{hic_map,[^/]+}{category,.*/}{gc,[^./-]+}.shapiro.txt' %STATS_DIR,
        plot = '%s/{hic_map}{category}{gc}.shapiro.pdf' %STATS_DIR
    log:
        'logs{category}calculateShapiroPvalue_{hic_map}_{gc}.log'
    shell:
        '%s/calculateShapiroPvalue.py -r {params.repeats} ' %SCRIPT_DIR +
        '-p {output.plot} -a {params.alpha} {input.hic_map} {params.coords} '
        '> {output.pval} 2> {log} || true'


rule calc_pvalues_shapiro_pairwise:
    input:
        hic_map = lambda wildcards: '%s/%s:%s.cool' %(
                POOLED_MAPS_DIR, wildcards.hic_map, wildcards.hic_map in
                MANAGER.single_pools and POOL_STRAT_REP or POOL_STRAT_CON)
    params:
        repeats = SAMPLE_SIZE,
        alpha = PVALUE_ALPHA,
        coords1 = lambda wildcards: MANAGER.rois[MANAGER.filename2id[\
                wildcards.gc1]][2],
        coords2 = lambda wildcards: MANAGER.rois[MANAGER.filename2id[\
                wildcards.gc2]][2],
    output:
        pval = '%s/{hic_map,[^/]+}{category,.*/}{gc1,[^/-]+}-{gc2,[^/-]+}.shapiro.txt' %STATS_DIR,
        plot = '%s/{hic_map}{category}{gc1}-{gc2}.shapiro.pdf' %STATS_DIR
    log:
        'logs{category}calculateShapiroPvalue_{hic_map}_{gc1}-{gc2}.log'
    shell:
        '%s/calculateShapiroPvalue.py -r {params.repeats} ' %SCRIPT_DIR +
        '-p {output.plot} -a {params.alpha} {input.hic_map} {params.coords1} '
        '{params.coords2} > {output.pval} 2> {log} || true'


rule interaction_summary_table:
    input:
        lambda wildcards: map(lambda x: '%s/%s.txt' %(STATS_DIR, x),
                MANAGER.comparisonFiles(pool=wildcards.pool,
                coi=wildcards.coi.find('-') == -1 and (wildcards.coi,
                wildcards.coi) or wildcards.coi.split('-')))
    output:
        '%s/{pool,[^/]+}/{coi,[^.]+}.csv' %STATS_DIR,
    shell:
        'echo -e \'compared region(s)\\tp-value\\t'
        'sumcount\\t#high-int. contacts\\t%high-int contacts\\t'
        'specific interacting regions\' > {output};'
        'grep high {input} | grep -v \'not significant\' | sed \'s/\\.txt:'
        'high contact region//\' | sed \'s/^.*\\///\' | cut -f1,2,4- >> {output} || true'


rule interaction_summary_table_shapiro:
    input:
        lambda wildcards: map(lambda x: '%s/%s.shapiro.txt' %(STATS_DIR, x),
                MANAGER.comparisonFiles(pool=wildcards.pool,
                coi=wildcards.coi.find('-') == -1 and (wildcards.coi,
                wildcards.coi) or wildcards.coi.split('-')))
    output:
        '%s/{pool,[^/]+}/{coi,[^.]+}.shapiro.csv' %STATS_DIR,
    shell:
        'echo -e \'compared region(s)\\tp-value\\t'
        'sumcount\\t#high-int. contacts\\t%high-int contacts\\t'
        'specific interacting regions\' > {output};'
        'grep high {input} | grep -v \'not significant\' | sed \'s/\\.txt:'
        'high contact region//\' | sed \'s/^.*\\///\' | cut -f1-2,4- '
        '>> {output} || true'


rule circos_annotation:
    input:
        MANAGER.getAnnotationFile()
    output:
        '%s/annotation.txt' %PLOT_DIR
    log:
        'logs/annotation2circos_%s.log' %basename(MANAGER.getAnnotationFile())
    shell:
        '%s/annotation2circos.py {input} > {output} 2> {log}' %SCRIPT_DIR


rule circos_expression:
    input:
        annot = MANAGER.getAnnotationFile(),
        expr = MANAGER.expression_file
    output:
        '%s/expression.txt' %PLOT_DIR
    log:
        'logs/deseq2circos_%s.log' %basename(MANAGER.expression_file)
    shell:
        '%s/deseq2circos.py {input.annot} {input.expr} > {output} ' %SCRIPT_DIR +
        '2> {log}'


rule plot_spatial_diff_circos:
    input:
        hic_map1 = '%s/{hic_map1}:{strategy}.cool' %POOLED_MAPS_DIR,
        hic_map2 = '%s/{hic_map2}:{strategy}.cool' %POOLED_MAPS_DIR,
        annotation = '%s/annotation.txt' %PLOT_DIR
    params:
        outdir = '%s/{hic_map1}-{hic_map2}/top{topx}' %PLOT_DIR +
                '_peak{interaction_threshold}',
        centros = list(map(lambda x: MANAGER.rois[x][2],
                MANAGER.category_members.get('chromocenter', ()))),
        chroms = MANAGER.chromosomes,
        it = '{interaction_threshold}',
        tx = '{topx}'
    output:
        circos1 = '%s/{hic_map1,[^/]+}-{hic_map2,[^/]+}/' %PLOT_DIR +
                'top{topx}_peak{interaction_threshold,[0-9.]+}/'
                '{hic_map1}:{strategy,[^_]+}.circos.conf',
        circos2 = '%s/{hic_map1,[^/]+}-{hic_map2,[^/]+}/' %PLOT_DIR +
                'top{topx}_peak{interaction_threshold,[0-9.]+}/'
                '{hic_map2}:{strategy,[^_]+}.circos.conf',
    log:
        'logs/topDiffPlot_{hic_map1}-{hic_map2}:{strategy}_top{topx}_'
        'peak{interaction_threshold}.log'
    shell:
        '%s/topDiffPlot.py -q {params.it} -c {params.chroms} ' %SCRIPT_DIR +
        '-i {params.centros} -o {params.outdir} -- {input.hic_map1} '
        '{input.hic_map2} {input.annotation} {params.tx} 2> {log}'


rule plot_spatial_diff_circos_zoom:
    input:
        hic_map1 = '%s/{hic_map1}:{strategy}.cool' %POOLED_MAPS_DIR,
        hic_map2 = '%s/{hic_map2}:{strategy}.cool' %POOLED_MAPS_DIR,
        annotation = '%s/annotation.txt' %PLOT_DIR,
        expression = '%s/expression.txt' %PLOT_DIR
    params:
        outdir = '%s/{hic_map1}-{hic_map2}/top{topx}' %PLOT_DIR +
                '_peak{interaction_threshold}_{gc}',
        centros = list(map(lambda x: MANAGER.rois[x][2],
                MANAGER.category_members.get('chromocenter', ()))),
        chroms = MANAGER.chromosomes,
        it = '{interaction_threshold}',
        tx = '{topx}',
        gc_coords = lambda wildcards: MANAGER.rois[MANAGER.filename2id[\
                wildcards.gc]][2],
        gc_genes = lambda wildcards: MANAGER.rois[MANAGER.filename2id[\
                wildcards.gc]][3]
    output:
        circos1 = '%s/{hic_map1,[^/]+}-{hic_map2,[^/]+}/' %PLOT_DIR +
                'top{topx}_peak{interaction_threshold,[0-9.]+}_{gc}/'
                '{hic_map1}:{strategy}.circos.conf',
        circos2 = '%s/{hic_map1,[^/]+}-{hic_map2,[^/]+}/' %PLOT_DIR +
                'top{topx}_peak{interaction_threshold,[0-9.]+}_{gc}/'
                '{hic_map2}:{strategy}.circos.conf',
    log:
        'logs/topDiffPlot_{gc}_{hic_map1}-{hic_map2}:{strategy}_top{topx}_'
        'peak{interaction_threshold}.log'
    shell:
        '%s/topDiffPlot.py -q {params.it} -c {params.chroms} ' %SCRIPT_DIR +
        '-i {params.centros} -z {params.gc_coords} -o {params.outdir} -g '
        '{params.gc_genes} -e {input.expression} -- {input.hic_map1} '
        '{input.hic_map2} {input.annotation} {params.tx} 2> {log}'


rule do_circos:
    input:
        circos = '%s/{prefix}/top{topx}_peak{it}/{circos_file}.circos.conf' %PLOT_DIR,
        annotation = '%s/annotation.txt' %PLOT_DIR,
        expression = '%s/expression.txt' %PLOT_DIR
    output:
        '%s/{prefix,[^/]+}/top{topx,[0-9.]+}_peak{it}/{circos_file}.png' %PLOT_DIR
    log:
        'logs/circos_{prefix}_top{topx}_peak{it}_{circos_file}.log'
    shell:
        '%s -conf {input.circos} > {log}' %CIRCOS_BIN

#rule centromeric_interaction_beds:
#    input:
#        '%s/{hic_map}_high_interaction.csv' %STATS_DIR
#    output:
#        temp('%s/{hic_map}/dummy' %BEDS_DIR)
#    shell:
#        'touch {output};'
#        '%s/stats_to_bed.py {input} %s/{wildcards.hic_map}' %(SCRIPT_DIR, BEDS_DIR)


#rule plot_genecluster_summary:
#    input:
#        expand('%s/{{hic_map}}/{gc}.txt' %STATS_DIR, gc=[x for x in GCS if not
#                x.lower().startswith('positive') and not
#                x.lower().startswith('negative')]),
#    params:
#        type='{type}',
#        minbound = -0.003,
#        maxbound = 0.003
#    output:
#        '%s/{hic_map}_genecluster_{type}_interaction.pdf' %STATS_DIR
#    shell:
#        '%s/plotGeneclusterSummary.py -t {params.type} ' %SCRIPT_DIR +
#        '-y {params.minbound} {params.maxbound} {input} > {output}'
#
#
#rule plot_chromocenter_summary:
#    input:
#        expand('%s/{{hic_map}}/{{gc}}_{centromere}.txt' %STATS_DIR,
#                centromere=CENTROMERES),
#    params:
#        type='{type}',
#        minbound = -0.04,
#        maxbound = 0.04
#    output:
#        '%s/{hic_map}_{type}_interaction_{gc}.pdf' %STATS_DIR
#    shell:
#        '%s/plotChromocenterSummary.py -t {params.type} ' %SCRIPT_DIR +
#        '-y {params.minbound} {params.maxbound} {input} > {output}'


