import pandas as pd
import os
from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
S3 = S3RemoteProvider()

## S3 root bucket
SROOT=config['s3root']

## vg docker container
if 'vg_docker' not in config:
    config['vg_docker'] = 'quay.io/vgteam/vg:v1.29.0'

## samples to use, based on the 'dataset' label
SAMPS = ""
if 'dataset' not in config:
    config['dataset'] = 'hpp60'
if config['dataset'] == 'lc2019_12ont':
    SAMPS = 'CHM13 GM24143 GM24149 GM24385 HG00733 HG01109 HG01243 HG02055 HG02080 HG02723 HG03098 HG03492'.split()
if config['dataset'] == 'lc2019_4ont':
    SAMPS = 'CHM13 HG00733 GM24385 HG01109'.split()
if config['dataset'] == 'hpp60':
    SAMPS = 'HG00436 HG00437 HG00619 HG00620 HG00671 HG00672 HG00734 HG00739 HG00740 HG01047 HG01069 HG01070 HG01104 HG01105 HG01121 HG01122 HG01173 HG01174 HG01356 HG01357 HG01889 HG01890 HG01926 HG01927 HG01950 HG01951 HG01976 HG01977 HG02146 HG02147 HG02255 HG02256 HG02484 HG02485 HG02557 HG02558 HG02570 HG02571 HG02620 HG02621 HG02628 HG02629 HG02715 HG02716 HG02884 HG02885 HG03451 HG03452 HG03514 HG03515 HG03538 HG03539 HG03577 HG03578 HG01256 HG01257 HG01359 HG01360 HG003 HG004 CHM13'.split()

## chromosome currently analyzed
CHR = 'chr20'
if 'chr' in config:
    CHR = config['chr']

## reads
READS=['HG002-glenn_giabfeb26']
READS_LR=['HG002-SequelII-merged_15kb_20kb-pbmm2']

## List experiments (pangenomes) in the form {method}/{params}/{dataset}.{method}.{params}
##  and {method}/{params}/map/{dataset}.{method}.{params}.{read}.{mapper}
## Note: because some pangenomes use embedded paths and don't, we need this loop and
##       can't use 'expand' easily to deal with different {mapper} mode used on different pangenomes
EMB_PATH_METHS = ['cactus', 'seqwish', 'vg']   ## methods where we should use embedded paths for mapping (for the 'auto' option)
EXPS = []  ## {method}/{params}/{dataset}.{method}.{params}
EXPS_MAP = [] ## {method}/{params}/map/{dataset}.{method}.{params}.{read}.{mapper}
EXPS_MAP_LR = [] ## {method}/{params}/map/{dataset}.{method}.{params}.{read}.{mapper}
EXPS_SUBGRAPHS = [] ## {method}/{params}/subgraphs/{dataset}.{method}.{params}
EXP_LABELS = [] ## 
dataset_s = config['dataset'].split()
exp_s = config['exp'].split()
for ii in range(len(exp_s)):
    exp = exp_s[ii]
    if len(dataset_s) == len(exp_s):
        DATASET = '{}-hg38'.format(dataset_s[ii])
        EXP_LABELS.append(DATASET + '_' + exp)
    else:
        DATASET = '{}-hg38'.format(dataset_s[0])
        EXP_LABELS.append(exp)
    ## Made from 'exp' config that contains {method}/{params}, and the DATASET value computed above
    EXPS.append(exp + '/' + DATASET + '.' + exp.replace('/', '.'))
    if exp != 'vg/linear':
        # no need to visualize the linear reference "graph"
        EXPS_SUBGRAPHS.append(exp + '/subgraphs/' + DATASET + '.' + exp.replace('/', '.'))
    exp_meth = exp.split('/')[0]
    ## eventually match the MAPPER to pangenome ('auto' mode) or non-default mapping mode
    MAPPER=config['mapper']
    if MAPPER == 'giraffe':
        if config['cover_n'] == 'paths' or (config['gbwt_auto'] and exp_meth in EMB_PATH_METHS):
            ## use embedded paths
            MAPPER = 'giraffe{}k{}wPaths'.format(config['min_k'], config['min_w'])
        else:
            ## use greedy algo to make a path cover from the pangenome
            MAPPER = 'giraffe{}k{}w{}N'.format(config['min_k'], config['min_w'], config['cover_n'])    
    for read in READS:
        EXPS_MAP.append(exp + '/map/' + DATASET + '.' + exp.replace('/', '.') + '.' + read + '.' + MAPPER)
    for read in READS_LR:
        EXPS_MAP_LR.append(exp + '/map/' + DATASET + '.' + exp.replace('/', '.') + '.' + read + '.' + config['mapper_lr'])

## Regions to visualize
regs = pd.read_csv(config['viz_bed'], sep='\t', header=None)
REGIONS = []
for ii in range(len(regs)):
    if regs[0][ii] == CHR:
        REGIONS.append('hg38_' + CHR + '-' + str(regs[1][ii]) + '-' + str(regs[2][ii]))

## if the input graph didn't follow the naming conventions, define "links" in this file
## CSV: s3 path, new s3 path (that follows naming conventions)
in_ln = {}
if 's3_input_links' in config:
    s3ln = pd.read_csv(config['s3_input_links'], header=None)
    for ii in range(len(s3ln)):
        in_ln[s3ln[1][ii]] = s3ln[0][ii]

rule eval_all:
    input:
        rmd='evaluation-report.Rmd',
        graph=S3.remote(expand('{sroot}/{exp}.vgstats.txt', sroot=SROOT, exp=EXPS)),
        degree=S3.remote(expand('{sroot}/{exp}.vgstats-degree.tsv', sroot=SROOT, exp=EXPS)),
        diststats=S3.remote(expand('{sroot}/{exp}.dist-stats.tsv.gz', sroot=SROOT, exp=EXPS)),
        map=S3.remote(expand('{sroot}/{exp}.stats.txt', sroot=SROOT, exp=EXPS_MAP)),
        maplr=S3.remote(expand('{sroot}/{exp}.stats.txt', sroot=SROOT, exp=EXPS_MAP_LR)),
        decon=expand('{exp}.decon.altsplit.vcf', exp=EXPS)
    output: S3.remote(SROOT + '/' + config['html_out'])
    params: exps=EXP_LABELS
    shell:
        """
        echo {params.exps} > files.info
        echo graph {input.graph} >> files.info
        echo degree {input.degree} >> files.info
        echo diststats {input.diststats} >> files.info
        echo map {input.map} >> files.info
        echo maplr {input.maplr} >> files.info
        echo decon {input.decon} >> files.info
        Rscript -e 'rmarkdown::render("evaluation-report.Rmd")'
        cp evaluation-report.html {output}
        """

rule eval_srmap_vgstats:
    input:
        rmd='evaluation-report.Rmd',
        graph=S3.remote(expand('{sroot}/{exp}.vgstats.txt', sroot=SROOT, exp=EXPS)),
        degree=S3.remote(expand('{sroot}/{exp}.vgstats-degree.tsv', sroot=SROOT, exp=EXPS)),
        map=S3.remote(expand('{sroot}/{exp}.stats.txt', sroot=SROOT, exp=EXPS_MAP))
    output: S3.remote(SROOT + '/' + config['html_out'])
    params: exps=EXP_LABELS
    shell:
        """
        echo {params.exps} > files.info
        echo graph {input.graph} >> files.info
        echo degree {input.degree} >> files.info
        echo map {input.map} >> files.info
        Rscript -e 'rmarkdown::render("evaluation-report.Rmd")'
        cp evaluation-report.html {output}
        """

rule eval_srmap:
    input:
        rmd='evaluation-report.Rmd',
        map=S3.remote(expand('{sroot}/{exp}.stats.txt', sroot=SROOT, exp=EXPS_MAP))
    output: S3.remote(SROOT + '/' + config['html_out'])
    params: exps=EXP_LABELS
    shell:
        """
        echo {params.exps} > files.info
        echo map {input.map} >> files.info
        Rscript -e 'rmarkdown::render("evaluation-report.Rmd")'
        cp evaluation-report.html {output}
        """

rule viz:
    input:
        rmd='visualization-report.Rmd',
        bed=config['viz_bed'],
        vgview=S3.remote(expand('{sroot}/{exp}.{region}_c50.dot.svg', sroot=SROOT, exp=EXPS_SUBGRAPHS, region=REGIONS)),
        odgiviz=S3.remote(expand('{sroot}/{exp}.{region}_c50.odgi.png', sroot=SROOT, exp=EXPS_SUBGRAPHS, region=REGIONS))
    output: S3.remote(SROOT + '/' + config['viz_html_out'])
    run:
        outf = open('files.info', 'w')
        outf.write('vgview {}\n'.format(' '.join(input.vgview)))
        outf.write('odgiviz {}\n'.format(' '.join(input.odgiviz)))
        outf.close()
        shell("Rscript -e 'rmarkdown::render(\"visualization-report.Rmd\")' files.info {input.bed}")
        shell("cp visualization-report.html {output}")

rule prep_viz:
    input:
        r=S3.remote(SROOT + '/visualization-app-prepdata.R'),
        bed=S3.remote(SROOT + '/' + config['viz_bed']),
        gfas=S3.remote(expand('{sroot}/{exp}.{region}_c50.gfa', sroot=SROOT, exp=EXPS_SUBGRAPHS, region=REGIONS))
    output: S3.remote(SROOT + '/viz-app-data.RData')
    params: exps=config['exp'].split()
    shell:
        """
        echo gfa {input.gfas} > files.info
        Rscript {input.r} files.info {input.bed} {output}
        """    

rule main_dev:
    input:
        rmd=S3.remote(SROOT + '/evaluation-report-dev.Rmd'),
        graph=S3.remote(expand('{sroot}/{exp}.vgstats.txt', sroot=SROOT, exp=EXPS)),
        degree=S3.remote(expand('{sroot}/{exp}.vgstats-degree.tsv', sroot=SROOT, exp=EXPS)),
        diststats=S3.remote(expand('{sroot}/{exp}.dist-stats.tsv.gz', sroot=SROOT, exp=EXPS)),
        map=S3.remote(expand('{sroot}/{exp}.stats.txt', sroot=SROOT, exp=EXPS_MAP)),
        decon=S3.remote(expand('{sroot}/{exp}.decon.vcf.gz', sroot=SROOT, exp=EXPS))
    output: S3.remote(SROOT + '/evaluation-report-dev.html')
    params: exps=config['exp'].split()
    shell:
        """
        echo {params.exps} > files.info
        echo graph {input.graph} >> files.info
        echo degree {input.degree} >> files.info
        echo diststats {input.diststats} >> files.info
        echo map {input.map} >> files.info
        echo decon {input.decon} >> files.info
        cp {input.rmd} evaluation-report.Rmd
        Rscript -e 'rmarkdown::render("evaluation-report.Rmd")'
        cp evaluation-report.html {output}
        """

##
## Misc
##

rule gzip:
    input: '{file}'
    output: S3.remote(SROOT + '/{file}.gz')
    threads: 8
    shell:
        'pigz -c -p {threads} {input} > {output}'

rule bgzip_fa:
    input: '{file}.fa'
    output: S3.remote(SROOT + '/{file}.fa.gz')
    shell:
        'bgzip -c {input} > {output}'

# rule bgzip:
#     input: '{file}'
#     output: S3.remote(SROOT + '/{file}.bgz')
#     shell:
#         'bgzip -c {input} > {output}'

# rule bgzip_vcf:
#     input: '{vcf}.vcf'
#     output:
#         bgz=S3.remote(SROOT + '/{vcf}.vcf.bgz'),
#         tbi=S3.remote(SROOT + '/{vcf}.vcf.bgz.tbi')
#     shell:
#         """
#         bcftools sort {input} | bgzip > {output.bgz}
#         tabix {output.bgz}
#         """

rule bgzip_vcf:
    input: '{vcf}.vcf'
    output:
        bgz=S3.remote(SROOT + '/{vcf}.vcf.bgz'),
        tbi=S3.remote(SROOT + '/{vcf}.vcf.bgz.tbi')
    params:
        temp_vcf='{vcf}.temp.vcf'
    run:
        ## check if empty VCF
        empty_vcf = True
        with open(input[0], 'r') as invcf:
            for line in invcf:
                if line[0] != '#':
                    empty_vcf = False
                    break
        ## sort, compress and index
        shell('head -10000 {input} | grep "^#" >> {params.temp_vcf}')
        if not empty_vcf:
            shell('grep -v "^#" {input} | sort -k1,1d -k2,2n >> {params.temp_vcf}')
        shell('bgzip -c {params.temp_vcf} > {output.bgz}')
        shell('rm {params.temp_vcf}')
        shell('tabix -f {output.bgz}')

rule rename_contigs_ref:
    input: S3.remote(SROOT + '/hg38_' + CHR + '.fa.gz')
    output: 'hg38_' + CHR + '.renamed.fa'
    shell:
        """
        zcat {input} | sed 's/^>/>hg38_/' > {output}
        """

rule rename_contigs:
    input: S3.remote(SROOT + '/{dataset}_fasta/{samp}-' + CHR + '.fa.gz')
    output: '{dataset}_fasta/{samp}-' + CHR + '.renamed.fa'
    shell:
        """
        zcat {input} | sed 's/^>/>{wildcards.samp}_/' > {output}
        """

##
## Pangenome construction using seqwish
##

rule merge_fasta:
    input: ref=S3.remote(SROOT + '/hg38_' + CHR + '.renamed.fa.gz'),
           amb=S3.remote(expand(SROOT + '/{{dataset}}_fasta/{samp}-' + CHR + '.renamed.fa.gz', samp=SAMPS))
    output: '{dataset}-hg38.fa'
    shell:
        """
        zcat {input.ref} > {output}
        zcat {input.amb} >> {output}
        """

rule minimap2:
    input: S3.remote(SROOT + '/{dataset}.fa.gz')
    output: 'paf/{dataset}.asm{asmX,\d+}.paf'
    threads: 16
    benchmark: S3.remote(SROOT + '/benchmarks/{dataset}.minimap2.asm{asmX}.minimap2.benchmark.tsv')
    log: S3.remote(SROOT + '/logs/{dataset}.minimap2.asm{asmX}.minimap2.log')
    shell:
        "minimap2 -t {threads} -K 100M -r 10000 -c --cs -X -x asm{wildcards.asmX} {input} {input} > {output} 2> {log}"

rule filterpaf_fpa:
    input: S3.remote(SROOT + '/paf/{dataset}.asm{asmX}.paf.gz')
    output: 'paf/{dataset}.asm{asmX}.dropl{l}.paf'
    wildcard_constraints:
        l="\d+",
        asmX="\d+"
    benchmark: S3.remote(SROOT + '/benchmarks/{dataset}.fpa.asm{asmX}-dropl{l}.fpa.benchmark.tsv')
    log: S3.remote(SROOT + '/logs/{dataset}.fpa.asm{asmX}-dropl{l}.fpa.log')
    shell:
        'fpa -i {input} -o {output} -z no drop -l {wildcards.l} 2> {log}'

rule filterpaf_fpa_edyeet:
    input: S3.remote(SROOT + '/paf/{dataset}.{pafal}.paf.gz')
    output: 'paf/{dataset}.{pafal}-dropl{l}.paf'
    wildcard_constraints:
        l="\d+",
        pafal="edyeet_s\d+"
    benchmark: S3.remote(SROOT + '/benchmarks/{dataset}.fpa.{pafal}-dropl{l}.fpa.benchmark.tsv')
    log: S3.remote(SROOT + '/logs/{dataset}.fpa.{pafal}-dropl{l}.fpa.log')
    shell:
        'fpa -i {input} -o {output} -z no drop -l {wildcards.l} 2> {log}'

rule filterpaf_primary:
    input: 'paf/{dataset}.asm{asmX}.dropl{l}.paf'
    output: 'paf/{dataset}.asm{asmX}.dropl{l}-tpP.paf'
    wildcard_constraints:
        l="\d+",
        asmX="\d+"
    benchmark: S3.remote(SROOT + '/benchmarks/{dataset}.filter.asm{asmX}-dropl{l}.tpP.benchmark.tsv')
    shell:
        'grep "tp:A:P" {input} > {output}'

rule filterpaf_noselfalign:
    input: 'paf/{dataset}.{methpar}.paf'
    output: 'paf/{dataset}.{methpar}-noself.paf'
    benchmark: S3.remote(SROOT + '/benchmarks/{dataset}.filter.{methpar}.noself.benchmark.tsv')
    shell:
        "awk '{{if($1!=$6){{print $0}}}}' {input} > {output}"

rule filterpaf_ol:
    input:
        paf=S3.remote(SROOT + '/paf/{dataset}.asm{asmX}.sorted.paf.gz'),
        rscript=S3.remote(SROOT + '/filterPafByOverlaps.R'),
    output: S3.remote(SROOT + '/paf/{dataset}.asm{asmX}.ol{olmeth}.paf.gz')
    benchmark: S3.remote(SROOT + '/benchmarks/{dataset}.filter.asm{asmX}.ol{olmeth}.benchmark.tsv')
    params:
        temp_paf='{dataset}.asm{asmX}.ol{olmeth}.paf'
    threads: 32
    shell:
        """
        Rscript {input.rscript} {input.paf} {params.temp_paf} {threads} ol{wildcards.olmeth}
        pigz -c  -p {threads} {params.temp_paf} > {output}
        rm {params.temp_paf}
        """

rule seqwish_kl:
    input:
        fa=S3.remote(SROOT + '/{dataset}.fa.gz'),
        paf='paf/{dataset}.asm{asmX}.{filterpaf}.paf'
    output: 'seqwish/asm{asmX}-{filterpaf}-k{k}-l{l}/{dataset}.seqwish.asm{asmX}-{filterpaf}-k{k}-l{l}.gfa'
    wildcard_constraints:
        l="\d+",
        asmX="\d+"
    benchmark: S3.remote(SROOT + '/benchmarks/{dataset}.seqwish.asm{asmX}-{filterpaf}-k{k}-l{l}.seqwish.benchmark.tsv')
    log: S3.remote(SROOT + '/logs/{dataset}.seqwish.asm{asmX}-{filterpaf}-k{k}-l{l}.seqwish.log')
    threads: 16
    singularity:
        "docker://jmonlong/seqwish:9bbfa70"
    shell:
        'seqwish -P -t {threads} -k {wildcards.k} -l {wildcards.l} -s {input.fa} -p {input.paf} -g {output} 2> {log}'

rule seqwish_kl_gz:
    input:
        fa=S3.remote(SROOT + '/{dataset}.fa.gz'),
        paf=S3.remote(SROOT + '/paf/{dataset}.asm{asmX}.{filterpaf}.paf.gz')
    output: 'seqwish/asm{asmX}-{filterpaf}-k{k}-l{l}/{dataset}.seqwish.asm{asmX}-{filterpaf}-k{k}-l{l}.gfa'
    wildcard_constraints:
        l="\d+",
        asmX="\d+"
    benchmark: S3.remote(SROOT + '/benchmarks/{dataset}.seqwish.asm{asmX}-{filterpaf}-k{k}-l{l}.seqwish.benchmark.tsv')
    log: S3.remote(SROOT + '/logs/{dataset}.seqwish.asm{asmX}-{filterpaf}-k{k}-l{l}.seqwish.log')
    threads: 16
    singularity:
        "docker://jmonlong/seqwish:9bbfa70"
    shell:
        'seqwish -P -t {threads} -k {wildcards.k} -l {wildcards.l} -s {input.fa} -p {input.paf} -g {output} 2> {log}'

rule smooth_gfa:
    input: S3.remote(SROOT + '/seqwish/asm{asmX}-dropl{dropl}-k{k}-l{l}/{dataset}.seqwish.asm{asmX}-dropl{dropl}-k{k}-l{l}.gfa.gz')
    output: 'seqwish/asm{asmX}-dropl{dropl}-k{k}-l{l}-w{w}/{dataset}.seqwish.asm{asmX}-dropl{dropl}-k{k}-l{l}-w{w}.gfa'
    benchmark: S3.remote(SROOT + '/benchmarks/{dataset}.seqwish.asm{asmX}-dropl{dropl}-k{k}-l{l}-w{w}.smoothxg.benchmark.tsv')
    threads: 16
    params:
        unz_in_gfa='temp.{dataset}.seqwish.asm{asmX}-dropl{dropl}-k{k}-l{l}-w{w}.smoothxg.gfa'
    singularity:
        "docker://jmonlong/smoothxg:5308db3"
    shell:
        """
        gunzip -c {input} > {params.unz_in_gfa}
        echo smoothxg $SMOOTHXG_COMMIT
        smoothxg -g {params.unz_in_gfa} -w {wildcards.w} -t {threads} > {output}
        rm {params.unz_in_gfa}
        """

## edyeet version
rule edyeet:
    input: S3.remote(SROOT + '/{dataset}.fa.gz')
    output: 'paf/{dataset}.edyeet_s100000.paf'
    threads: 16
    benchmark: S3.remote(SROOT + '/benchmarks/{dataset}.edyeet.s100000.edyeet.benchmark.tsv')
    log: S3.remote(SROOT + '/logs/{dataset}.edyeet.s100000.edyeet.log')
    shell:
        "edyeet -t {threads} -X -p 85 -a 85 -n 10 -s 100000 {input} {input} > {output} 2> {log}"

rule seqwish_k_gz:
    input:
        fa=S3.remote(SROOT + '/{dataset}.fa.gz'),
        paf=S3.remote(SROOT + '/paf/{dataset}.{almeth}.paf.gz')
    output: 'seqwish/{almeth}-k{k}/{dataset}.seqwish.{almeth}-k{k}.gfa'
    benchmark: S3.remote(SROOT + '/benchmarks/{dataset}.seqwish.{almeth}-k{k}.seqwish.benchmark.tsv')
    log: S3.remote(SROOT + '/logs/{dataset}.seqwish.{almeth}-k{k}.seqwish.log')
    wildcard_constraints:
        k="\d+"
    threads: 16
    params:
        com=os.environ["SEQWISH_COMMIT"]
    singularity:
        "docker://jmonlong/seqwish:9bbfa70"
    shell:
        """
        echo seqwish {params.com}
        seqwish -s {input.fa} -p {input.paf} -t {threads} -k {wildcards.k} -g {output} -P 2> {log}
        """

rule smooth_gfa_ed:
    input: S3.remote(SROOT + '/seqwish/{almeth}-k{k}/{dataset}.seqwish.{almeth}-k{k}.gfa.gz')
    output: 'seqwish/{almeth}-k{k}-w{w}/{dataset}.seqwish.{almeth}-k{k}-w{w}.gfa'
    benchmark: S3.remote(SROOT + '/benchmarks/{dataset}.seqwish.{almeth}-k{k}-w{w}.smoothxg.benchmark.tsv')
    wildcard_constraints:
        k="\d+"
    threads: 16
    params:
        unz_in_gfa='temp.{dataset}.seqwish.{almeth}-k{k}-w{w}.smoothxg.gfa'
    singularity:
        "docker://jmonlong/smoothxg:5308db3"
    shell:
        """
        gunzip -c {input} > {params.unz_in_gfa}
        echo smoothxg $SMOOTHXG_COMMIT
        smoothxg -t {threads} -g {params.unz_in_gfa} -X 32 -U 0.1 -k 8 -w {wildcards.w} > {output}
        rm {params.unz_in_gfa}
        """

## pggb version
rule pggb:
    input: S3.remote(SROOT + '/{dataset}.fa.gz')
    output: 'seqwish/pggb-s{s}-k{k}-n{n}-pa{pa}/{dataset}.seqwish.pggb-s{s}-k{k}-n{n}-pa{pa}.gfa'
    threads: 16
    benchmark: S3.remote(SROOT + '/benchmarks/{dataset}.seqwish.pggb-s{s}-k{k}-n{n}-pa{pa}.benchmark.tsv')
    log: S3.remote(SROOT + '/logs/{dataset}.seqwish.pggb-s{s}-k{k}-n{n}-pa{pa}.log')
    singularity:
        "docker://jmonlong/pggb:970d2e5"
    shell:
        """
        pggb -i {input} -s {wildcards.s} -p {wildcards.pa} -a {wildcards.pa} -n {wildcards.n} -k {wildcards.k} -t {threads} 2> {log}
        mv {input}.*.smooth.gfa {output}
        """

# rule gfa_to_vg_seqwish:
#     input: S3.remote(SROOT + '/seqwish/{params}/{dataset}.seqwish.{params}.gfa.gz')
#     output: S3.remote(SROOT + '/seqwish/{params}/{dataset}.seqwish.{params}.vg')
#     threads: 4
#     benchmark: S3.remote(SROOT + '/benchmarks/{dataset}.seqwish.{params}.gfa-to-vg.benchmark.tsv')
#     singularity:
#         "docker://" + config['vg_docker']
#     shell:
#         'zcat {input} | vg convert -ga - | vg mod --unchop - | vg mod --chop 32 - > {output}'

##
## pangenome construction using minigraph
##

rule minigraph:
    input: S3.remote(SROOT + '/hg38_' + CHR + '.renamed.fa.gz'),
           S3.remote(expand(SROOT + '/{{dataset}}_fasta/{samp}-' + CHR + '.renamed.fa.gz', samp=SAMPS))
    output: 'minigraph/L50-l50k/{dataset}-hg38.minigraph.L50-l50k.gfa'
    threads: 16
    benchmark: S3.remote(SROOT + '/benchmarks/{dataset}-hg38.minigraph.L50-l50k.minigraph.benchmark.tsv')
    log: S3.remote(SROOT + '/logs/{dataset}-hg38.minigraph.L50-l50k.minigraph.log')
    shell:
        """
        minigraph -t {threads} -x ggs -L 50 -l 50k {input} -o {output} 2> {log}
        """

## for minigraph we also need to convert the segment names to numeric (e.g. s1 -> 1)
rule gfa_to_vg_minigraph:
    input: S3.remote(SROOT + '/minigraph/{params}/{dataset}.minigraph.{params}.gfa.gz')
    output: S3.remote(SROOT + '/minigraph/{params}/{dataset}.minigraph.{params}.vg')
    threads: 4
    benchmark: S3.remote(SROOT + '/benchmarks/{dataset}.minigraph.{params}.gfa-to-vg.benchmark.tsv')
    singularity:
        "docker://" + config['vg_docker']
    shell:
        "zcat {input} | python3 updateMinigraphGfa.py -r hg38_{CHR} | vg view -v -F - | vg mod --chop 32 - > {output}"

##
## Pangenome construction through VCF files using paftools
##

rule minimap2_ref_samp:
    input:
        ref=S3.remote(SROOT + '/hg38_' + CHR + '.renamed.fa.gz'),
        amb=S3.remote(SROOT + '/{dataset}_fasta/{samp}-' + CHR + '.renamed.fa.gz')
    output: 'paf/{dataset}-hg38.{samp}.asm{asmX,\d+}.paf'
    threads: 2
    benchmark: S3.remote(SROOT + '/benchmarks/{dataset}-hg38.minimap2.{samp}_asm{asmX}.minimap2.benchmark.tsv')
    log: S3.remote(SROOT + '/logs/{dataset}-hg38.minimap2.{samp}_asm{asmX}.minimap2.log')
    shell:
        "minimap2 -t {threads} -c --cs -X -x asm{wildcards.asmX} {input.ref} {input.amb} > {output} 2> {log}"

rule paftools_vcf:
    input:
        paf=S3.remote(SROOT + '/paf/{dataset}.{samp}.asm{asmX}.paf.gz'),
        ref=S3.remote(SROOT + '/hg38_' + CHR + '.renamed.fa.gz')
    output: 'paftools/asm{asmX}-L10k-l1k/{dataset}.paftools.asm{asmX}-L10k-l1k.{samp}.vcf'
    wildcard_constraints:
        asmX="\d+"
    threads: 2
    benchmark: S3.remote(SROOT + '/benchmarks/{dataset}.paftools.asm{asmX}-L10k-l1k-{samp}.paftools_call.benchmark.tsv')
    log: S3.remote(SROOT + '/logs/{dataset}.paftools.asm{asmX}-L10k-l1k-{samp}.paftools_call.log')
    shell:
        "zcat  {input.paf} | sort -k6,6 -k8,8n | paftools.js call -L 10000 -l 1000 -f {input.ref} -s {wildcards.samp} - > {output} 2> {log}"

rule merge_vcfs:
    input:
        bgz=S3.remote(expand(SROOT + '/paftools/{{params}}/{{dataset}}.paftools.{{params}}.{samp}.vcf.bgz', samp=SAMPS)),
        tbi=S3.remote(expand(SROOT + '/paftools/{{params}}/{{dataset}}.paftools.{{params}}.{samp}.vcf.bgz.tbi', samp=SAMPS))
    output: 'paftools/{params}/{dataset}.paftools.{params}.vcf'
    benchmark: S3.remote(SROOT + '/benchmarks/{dataset}.paftools.{params}.merge_vcfs.benchmark.tsv')
    log: S3.remote(SROOT + '/logs/{dataset}.paftools.{params}.merge_vcfs.log')
    shell:
        "bcftools merge -O v {input.bgz} > {output} 2> {log}"

rule add_ids_vcf:
    input: S3.remote(SROOT + '/paftools/{params}/{dataset}.paftools.{params}.{samp}.vcf.bgz')
    output: 'paftools/{params}/{dataset}.paftools.{params}.{samp}.ids.vcf'
    shell:
        "bcftools annotate -I '%CHROM\_%POS\_%END\_%QNAME\_%QSTART' {input} > {output}"

rule merge_vcf_svanalyzer:
    input:
        vcf=expand('paftools/{{params}}/{{dataset}}.paftools.{{params}}.{samp}.ids.vcf', samp=SAMPS),
        ref=S3.remote(SROOT + '/hg38_' + CHR + '.renamed.fa.gz')
    output: 'paftools/{params}-svanalyzer/{dataset}.paftools.{params}-svanalyzer.vcf'
    benchmark: S3.remote(SROOT + '/benchmarks/{dataset}.paftools.{params}-svanalyzer.merge_vcfs.benchmark.tsv')
    log: S3.remote(SROOT + '/logs/{dataset}.paftools.{params}-svanalyzer.merge_vcfs.log')
    params:
        temp_fof='temp_paftools.{dataset}.{params}_svanalyzer.fof',
        prefix='temp_out_paftools.{dataset}.{params}_svanalyzer'
    singularity:
        "docker://quay.io/biocontainers/svanalyzer:0.36--pl526_0"
    shell:
        """
        echo {input.vcf} | sed 's/ /\\n/g' > {params.temp_fof}
        svanalyzer merge --ref {input.ref} --fof {params.temp_fof} --prefix {params.prefix}
        mv {params.prefix}.clustered.vcf {output}
        """

rule merge_vcf_rol:
    input:
        rscript=S3.remote(SROOT + '/filterDupsVcf-rol.R'),
        vcf='paftools/{params}/{dataset}.paftools.{params}.vcf'
    output: 'paftools/{params}-rol50/{dataset}.paftools.{params}-rol50.vcf'
    benchmark: S3.remote(SROOT + '/benchmarks/{dataset}.paftools.{params}-rol50.merge_vcfs.benchmark.tsv')
    log: S3.remote(SROOT + '/logs/{dataset}.paftools.{params}-rol50.merge_vcfs.log')
    params:
        temp_vcf='temp_paftools.{dataset}.{params}_rol50.vcf'
    shell:
        """
        bcftools norm -m - {input.vcf} > {params.temp_vcf}
        Rscript {input.rscript} {params.temp_vcf} {output}
        rm {params.temp_vcf}
        """    

rule vg_construct_vcf:
    input:
        ref=S3.remote(SROOT + '/hg38_' + CHR + '.renamed.fa.gz'),
        vcf=S3.remote(SROOT + '/paftools/{params}/{dataset}.paftools.{params}.vcf.bgz'),
        vcf_tbi=S3.remote(SROOT + '/paftools/{params}/{dataset}.paftools.{params}.vcf.bgz.tbi')
    output: S3.remote(SROOT + '/paftools/{params}/{dataset}.paftools.{params}.vg')
    threads: 4
    benchmark: S3.remote(SROOT + '/benchmarks/{dataset}.paftools.{params}.vg_construct.benchmark.tsv')
    log: S3.remote(SROOT + '/logs/{dataset}.paftools.{params}.vg_construct.log')
    singularity:
        "docker://" + config['vg_docker']
    params:
        ref='hg38_' + CHR + '.renamed.fa'
    shell:
        """
        gunzip -c {input.ref} > {params.ref}
        vg construct -r {params.ref} -v {input.vcf} -a -S -t {threads} -p > {output} 2> {log}
        rm {params.ref}
        """

##
## Linear genome
##

rule linear:
    input: S3.remote(SROOT + '/hg38_' + CHR + '.renamed.fa.gz')
    output: S3.remote(SROOT + '/vg/linear/{dataset}.vg.linear.vg')
    params:
        ref='hg38_' + CHR + '.renamed.fa'
    singularity:
        "docker://" + config['vg_docker']
    shell:
        """
        gunzip -c {input} > {params.ref}
        vg construct -r {params.ref} > {output}
        rm {params.ref}
        """

##
## pangenome stats
##
# rule gfa_to_vg:
#     input: S3.remote(SROOT + '/{dir}/{file}.gfa.gz')
#     output: S3.remote(SROOT + '/{dir}/{file}.vg')
#     threads: 4
#     benchmark: S3.remote(SROOT + '/benchmarks/gfa_to_vg.{dir}.{file}.benchmark.tsv')
#     shell:
#         'zcat {input} | vg convert -g - | vg mod --chop 32 - > {output}'

def gfa_input(wildcards):
    in_path = '{method}/{params}/{dataset}.{method}.{params}.gfa.gz'.format(method=wildcards.method,
                                                                            params=wildcards.params,
                                                                            dataset=wildcards.dataset)
    if in_path in in_ln:
        return S3.remote(SROOT + '/' + in_ln[in_path])
    else:
        return S3.remote(SROOT + '/' + in_path)
rule gfa_to_vg:
    input: gfa_input
    output: S3.remote(SROOT + '/{method}/{params}/{dataset}.{method}.{params}.vg')
    threads: 4
    benchmark: S3.remote(SROOT + '/benchmarks/{dataset}.{method}.{params}.gfa-to-vg.benchmark.tsv')
    singularity:
        "docker://" + config['vg_docker']
    shell:
        'zcat {input} | vg convert -ga - | vg mod --unchop - | vg mod --chop 32 - > {output}'

def pg_input(wildcards):
    in_path = '{method}/{params}/{dataset}.{method}.{params}.pg'.format(method=wildcards.method,
                                                                        params=wildcards.params,
                                                                        dataset=wildcards.dataset)
    if in_path in in_ln:
        return S3.remote(SROOT + '/' + in_ln[in_path])
    else:
        return S3.remote(SROOT + '/' + in_path)
rule pg_to_vg:
    input: pg_input
    output: S3.remote(SROOT + '/{method}/{params}/{dataset}.{method}.{params}.vg')
    threads: 4
    benchmark: S3.remote(SROOT + '/benchmarks/{dataset}.{method}.{params}.pg-to-vg.benchmark.tsv')
    singularity:
        "docker://" + config['vg_docker']
    shell:
        'vg convert {input} | vg mod --chop 32 - > {output}'

rule vg_stats:
    input: S3.remote(SROOT + '/{method}/{params}/{dataset}.{method}.{params}.vg')
    output: S3.remote(SROOT + '/{method}/{params}/{dataset}.{method}.{params}.vgstats.txt')
    benchmark: S3.remote(SROOT + '/benchmarks/{dataset}.{method}.{params}.vg-stats.benchmark.tsv')
    singularity:
        "docker://" + config['vg_docker']
    shell:
        'vg stats -zl {input} > {output}'

rule vg_stats_degree:
    input: S3.remote(SROOT + '/{method}/{params}/{dataset}.{method}.{params}.vg')
    output: S3.remote(SROOT + '/{method}/{params}/{dataset}.{method}.{params}.vgstats-degree.tsv')
    benchmark: S3.remote(SROOT + '/benchmarks/{dataset}.{method}.{params}.vg-stats-degree.benchmark.tsv')
    singularity:
        "docker://" + config['vg_docker']
    shell:
        'vg stats -D {input} > {output}'

rule map_stats:
    input: S3.remote(SROOT + '/{method}/{params}/map/{dataset}.{method}.{params}.{read}.{mapper}.gaf.gz')
    output: S3.remote(SROOT + '/{method}/{params}/map/{dataset}.{method}.{params}.{read}.{mapper}.stats.txt')
    benchmark: S3.remote(SROOT + '/benchmarks/{dataset}.{method}.{params}.mapstats-{read}-{mapper}.benchmark.tsv')
    threads: 1
    shell:
        """
        gunzip -c {input} | awk '{{if($10=="*"){{$12=-1}};if($10!="*" && $10==$11){{perf="true"}}else{{perf="false"}};if($10!="*" && $10/$11>.99){{nn="true"}}else{{nn="false"}};print $12,perf,nn}}' | sort | uniq -c > {output}
        """

rule deconstruct_vcf_minigraph:
    input:
        py='gamToGFApath.py',
        vg=S3.remote(SROOT + '/minigraph/{params}/{dataset}.minigraph.{params}.vg'),
        xg=S3.remote(SROOT + '/minigraph/{params}/map/{dataset}.minigraph.{params}.xg'),
        snarls=S3.remote(SROOT + '/minigraph/{params}/map/{dataset}.minigraph.{params}.snarls'),
        gbwt=S3.remote(expand('{sroot}/minigraph/{{params}}/map/{{dataset}}.minigraph.{{params}}.N{n}.gbwt', sroot=SROOT, n=config['cover_n']))
    output: 'minigraph/{params}/{dataset}.minigraph.{params}.decon.vcf'
    benchmark: S3.remote(SROOT + '/benchmarks/{dataset}.minigraph.{params}.deconstruct.benchmark.tsv')
    log: S3.remote(SROOT + '/logs/{dataset}.minigraph.{params}.deconstruct.log')
    threads: 16
    params:
        gfa="temp_decon_{dataset}.minigraph.{params}.gfa",
        vgaug="temp_decon_{dataset}.minigraph.{params}.vg"
    singularity:
        "docker://" + config['vg_docker']
    shell:
        """
        vg view -g {input.vg} > {params.gfa}
        vg paths -X -x {input.xg} -g {input.gbwt} | vg view -a - | python3 {input.py} >> {params.gfa}
        vg convert -ga {params.gfa} > {params.vgaug}
        vg deconstruct -t {threads} -e -r {input.snarls} -P hg38 -P {CHR} {params.vgaug} > {output} 2> {log}
        """

rule deconstruct_vcf_paftools:
    input: S3.remote(SROOT + '/paftools/{params}/{dataset}.paftools.{params}.vcf.bgz')
    output: S3.remote(SROOT + '/paftools/{params}/{dataset}.paftools.{params}.decon.vcf.bgz')
    threads: 1
    shell:
        "cp {input} {output}"

rule deconstruct_vcf:
    input:
        xg=S3.remote(SROOT + '/{method}/{params}/map/{dataset}.{method}.{params}.xg'),
        snarls=S3.remote(SROOT + '/{method}/{params}/map/{dataset}.{method}.{params}.snarls')
    output: '{method}/{params}/{dataset}.{method}.{params}.decon.vcf'
    benchmark: S3.remote(SROOT + '/benchmarks/{dataset}.{method}.{params}.deconstruct.benchmark.tsv')
    log: S3.remote(SROOT + '/logs/{dataset}.{method}.{params}.deconstruct.log')
    threads: 16
    singularity:
        "docker://" + config['vg_docker']
    shell:
        "vg deconstruct -t {threads} -e -r {input.snarls} -P hg38 -P {CHR} {input.xg} > {output} 2> {log}"

## Deconstruct variants relative to the reference path
rule splitalts_deconstructed_vcf:
    input: S3.remote(SROOT + '/{exp}.decon.vcf.bgz')
    output: '{exp}.decon.altsplit.vcf'
    threads: 1
    shell:
        "bcftools norm -m - {input} > {output}"
    
rule dist_stats:
    input: S3.remote(SROOT + '/{method}/{params}/map/{dataset}.{method}.{params}.dist'),
    output: '{method}/{params}/{dataset}.{method}.{params}.dist-stats.tsv'
    singularity:
        "docker://" + config['vg_docker']
    shell:
        """
        echo -e "node_count\tdepth\tmin_length\tmax_length" > {output}
        vg view -B {input} | jq -r 'select(.type=="snarl") | select(.node_count>2) | [.node_count, .depth, .minimum_length, .maximum_length] | @tsv' >> {output}
        """

rule extract_region_dot:
    input: S3.remote(SROOT + '/{method}/{params}/map/{dataset}.{method}.{params}.xg')
    output: S3.remote(SROOT + '/{method}/{params}/subgraphs/{dataset}.{method}.{params}.{path}-{start}-{end}_c{context}.dot')
    singularity:
        "docker://" + config['vg_docker']
    shell:
        """
        CHRP=`vg paths -L -x {input} | grep -e hg38 -e "^{CHR}" | head -1`
        vg find -x {input} -c {wildcards.context} -p $CHRP:{wildcards.start}-{wildcards.end} | vg mod -Ou - | vg paths -r -Q $CHRP -v - | vg view -dup - > {output}
        """

rule extract_region_gfa:
    input: S3.remote(SROOT + '/{method}/{params}/map/{dataset}.{method}.{params}.xg')
    output: S3.remote(SROOT + '/{method}/{params}/subgraphs/{dataset}.{method}.{params}.{path}-{start}-{end}_c{context}.gfa')
    singularity:
        "docker://" + config['vg_docker']
    shell:
        """
        CHRP=`vg paths -L -x {input} | grep -e hg38 -e "^{CHR}" | head -1`
        vg find -x {input} -c {wildcards.context} -p $CHRP:{wildcards.start}-{wildcards.end} | vg mod -Ou - | vg view -g - > {output}')
        """

rule dot_to_png:
    input: S3.remote(SROOT + '/{method}/{params}/subgraphs/{dataset}.{method}.{params}.{path}-{start}-{end}_c{context}.dot')
    output: S3.remote(SROOT + '/{method}/{params}/subgraphs/{dataset}.{method}.{params}.{path}-{start}-{end}_c{context}.dot.png')
    shell:
        """
        dot -Tpng -o {output} {input}
        """

rule dot_to_svg:
    input: S3.remote(SROOT + '/{method}/{params}/subgraphs/{dataset}.{method}.{params}.{path}-{start}-{end}_c{context}.dot')
    output: S3.remote(SROOT + '/{method}/{params}/subgraphs/{dataset}.{method}.{params}.{path}-{start}-{end}_c{context}.dot.svg')
    shell:
        """
        dot -Tsvg -o {output} {input}
        """

rule odgi_viz:
    input: S3.remote(SROOT + '/{method}/{params}/subgraphs/{dataset}.{method}.{params}.{path}-{start}-{end}_c{context}.gfa')
    output: S3.remote(SROOT + '/{method}/{params}/subgraphs/{dataset}.{method}.{params}.{path}-{start}-{end}_c{context}.odgi.png')
    shell:
        """
        odgi build -g {input} -o - | odgi sort -i - -o - | odgi viz -i - -o {output} -x 1000 -y 500
        """

## Make test data
# rule gamsort:
#     input: S3.remote('s3://vg-k8s/users/jmonlong/giab6-v1.22/HG002/HG002-hs37d5-giab6.map.gam')
#     output:
#         gam=S3.remote('s3://vg-k8s/users/jmonlong/giab6-v1.22/HG002/HG002-hs37d5-giab6.map.sorted.gam'),
#         gai=S3.remote('s3://vg-k8s/users/jmonlong/giab6-v1.22/HG002/HG002-hs37d5-giab6.map.sorted.gam.gai')
#     threads: 8
#     shell:
#         "vg gamsort -t {threads} -p -i {output.gai} {input} > {output.gam}"

# rule real_reads_hg002:
#     input:
#         xg=S3.remote('s3://vg-k8s/users/jmonlong/giab6-v1.22/hs37d5-giab6.xg'),
#         gam=S3.remote('s3://vg-k8s/users/jmonlong/giab6-v1.22/HG002/HG002-hs37d5-giab6.map.sorted.gam'),
#         gai=S3.remote('s3://vg-k8s/users/jmonlong/giab6-v1.22/HG002/HG002-hs37d5-giab6.map.sorted.gam.gai')
#     output: S3.remote(SROOT + '/reads/HG002-sv_giab6.' + CHR + '.fastq.gz')
#     shell:
#         """
#         vg chunk -C -p ' + CHR + ' -x {input.xg} -a {input.gam} -b chunk
#         mv chunk_' + CHR + '.gam {output}
#         """

rule real_reads_hg002:
    input: S3.remote('s3://glennhickey/outstore/GIAB-FEB26/map-HG002/HG002_{chr}.gam')
    output: S3.remote(SROOT + '/reads/HG002-glenn_giabfeb26.chr{chr}.interleaved.fastq.gz')
    shell:
        """
        vg view -aX {input} | awk 'BEGIN{{delete ar; OFS="\\t"}}{{if(length(ar)==4){{print ar[1],ar[2],ar[3],ar[4]; delete ar}}; ar[length(ar)+1]=$0}}' | sort | awk '{{gsub("\\t","\\n", $0); print $0}}' | seqtk dropse | gzip > {output}
        """

rule dwl_hg002_ccs_bam:
    output:
        bam = 'reads/HG002.SequelII.merged_15kb_20kb.pbmm2.GRCh38.haplotag.10x.bam',
        bai = 'reads/HG002.SequelII.merged_15kb_20kb.pbmm2.GRCh38.haplotag.10x.bam.bai',
    shell:
        """
        wget -nv -O {output.bam} ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb_20kb_chemistry2/GRCh38/HG002.SequelII.merged_15kb_20kb.pbmm2.GRCh38.haplotag.10x.bam
        wget -nv -O {output.bai} ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb_20kb_chemistry2/GRCh38/HG002.SequelII.merged_15kb_20kb.pbmm2.GRCh38.haplotag.10x.bam.bai
        """

rule real_reads_hg002_ccs:
    input: 
        bam = 'reads/HG002.SequelII.merged_15kb_20kb.pbmm2.GRCh38.haplotag.10x.bam',
        bai = 'reads/HG002.SequelII.merged_15kb_20kb.pbmm2.GRCh38.haplotag.10x.bam.bai',
    output: S3.remote(SROOT + '/reads/HG002-SequelII-merged_15kb_20kb-pbmm2.{chr}.fastq.gz')
    shell:
        """
        samtools view -h {input.bam} {wildcards.chr} | samtools fastq - | gzip > {output}
        """


##
## Map reads using vg
##

# make xg index containing the alts paths
rule index_xg:
    input: S3.remote(SROOT + '/{method}/{params}/{dataset}.{method}.{params}.vg')
    output: S3.remote(SROOT + '/{method}/{params}/map/{dataset}.{method}.{params}.xg')
    threads: 16
    benchmark: S3.remote(SROOT + '/benchmarks/{dataset}.{method}.{params}.mapvg-xg.benchmark.tsv')
    singularity:
        "docker://" + config['vg_docker']
    shell:
        'vg index -L -t {threads} -x {output} {input}'

## GBWT with greedy path cover
rule index_gbwt_greedy:
    input: S3.remote(SROOT + '/{method}/{params}/map/{dataset}.{method}.{params}.xg')
    output:
        gg=S3.remote(SROOT + '/{method}/{params}/map/{dataset}.{method}.{params}.N{n}.gg'),
        gbwt=S3.remote(SROOT + '/{method}/{params}/map/{dataset}.{method}.{params}.N{n}.gbwt')
    threads: 1
    benchmark: S3.remote(SROOT + '/benchmarks/{dataset}.{method}.{params}.mapvg-gbwt-N{n}.benchmark.tsv')
    singularity:
        "docker://" + config['vg_docker']
    shell:
        'vg gbwt -n {wildcards.n} -g {output.gg} -o {output.gbwt} -x {input} -P'

## GBWT from the embedded paths in the pangenome
rule index_gbwt_paths:
    input: S3.remote(SROOT + '/{method}/{params}/{dataset}.{method}.{params}.vg')
    output: S3.remote(SROOT + '/{method}/{params}/map/{dataset}.{method}.{params}.paths.gbwt')
    threads: 1
    benchmark: S3.remote(SROOT + '/benchmarks/{dataset}.{method}.{params}.mapvg-gbwt-paths.benchmark.tsv')
    singularity:
        "docker://" + config['vg_docker']
    shell:
        'vg index -T -G {output} {input}'

rule index_minimizer:
    input:
        xg=S3.remote(SROOT + '/{method}/{params}/map/{dataset}.{method}.{params}.xg'),
        gbwt=S3.remote(SROOT + '/{method}/{params}/map/{dataset}.{method}.{params}.{n}.gbwt')
    output: S3.remote(SROOT + '/{method}/{params}/map/{dataset}.{method}.{params}.k{k}.w{w}.{n}.min')
    threads: 16
    benchmark: S3.remote(SROOT + '/benchmarks/{dataset}.{method}.{params}.mapvg-min-k{k}-w{w}-{n}.benchmark.tsv')
    singularity:
        "docker://" + config['vg_docker']
    shell:
        'vg minimizer -k {wildcards.k} -w {wildcards.w} -t {threads} -i {output} -g {input.gbwt} {input.xg}'

rule index_trivial_snarls:
    input: S3.remote(SROOT + '/{method}/{params}/map/{dataset}.{method}.{params}.xg')
    output: S3.remote(SROOT + '/{method}/{params}/map/{dataset}.{method}.{params}.trivial.snarls')
    threads: 16
    benchmark: S3.remote(SROOT + '/benchmarks/{dataset}.{method}.{params}.mapvg-trivialsnarls.benchmark.tsv')
    singularity:
        "docker://" + config['vg_docker']
    shell:
        'vg snarls -t {threads} --include-trivial -A integrated {input} > {output}'

rule index_snarls:
    input: S3.remote(SROOT + '/{method}/{params}/map/{dataset}.{method}.{params}.xg')
    output: S3.remote(SROOT + '/{method}/{params}/map/{dataset}.{method}.{params}.snarls')
    threads: 16
    benchmark: S3.remote(SROOT + '/benchmarks/{dataset}.{method}.{params}.mapvg-snarls.benchmark.tsv')
    singularity:
        "docker://" + config['vg_docker']
    shell:
        'vg snarls -t {threads} -A integrated -m 1000 {input} > {output}'

rule index_distance:
    input:
        vg=S3.remote(SROOT + '/{method}/{params}/{dataset}.{method}.{params}.vg'),
        snarls=S3.remote(SROOT + '/{method}/{params}/map/{dataset}.{method}.{params}.trivial.snarls')
    output: S3.remote(SROOT + '/{method}/{params}/map/{dataset}.{method}.{params}.dist')
    threads: 16
    benchmark: S3.remote(SROOT + '/benchmarks/{dataset}.{method}.{params}.mapvg-dist.benchmark.tsv')
    singularity:
        "docker://" + config['vg_docker']
    shell:
        'vg index -t {threads} -j {output} -s {input.snarls} {input.vg}'

rule prune_vg:
    input: S3.remote(SROOT + '/{method}/{params}/{dataset}.{method}.{params}.vg')
    output: S3.remote(SROOT + '/{method}/{params}/map/{dataset}.{method}.{params}.pruned.vg')
    threads: 1
    benchmark: S3.remote(SROOT + '/benchmarks/{dataset}.{method}.{params}.mapvg-prune.benchmark.tsv')
    log: S3.remote(SROOT + '/logs/{dataset}.{method}.{params}.mapvg-prune.log')
    singularity:
        "docker://" + config['vg_docker']
    shell:
        'vg prune -t {threads} -M 32 --restore-paths {input} > {output} 2> {log}'

rule index_gcsa:
    input: S3.remote(SROOT + '/{method}/{params}/map/{dataset}.{method}.{params}.pruned.vg')
    output:
        gcsa=S3.remote(SROOT + '/{method}/{params}/map/{dataset}.{method}.{params}.gcsa'),
        gcsalcp=S3.remote(SROOT + '/{method}/{params}/map/{dataset}.{method}.{params}.gcsa.lcp')
    threads: 16
    benchmark: S3.remote(SROOT + '/benchmarks/{dataset}.{method}.{params}.mapvg-gcsa.benchmark.tsv')
    log: S3.remote(SROOT + '/logs/{dataset}.{method}.{params}.mapvg-gcsa.log')
    params:
        tmp_dir="temp_gsca_{method}_{params}_{dataset}"
    singularity:
        "docker://" + config['vg_docker']
    shell:
        """
        mkdir -p {params.tmp_dir}
        vg index --temp-dir {params.tmp_dir} -p -t {threads} -g {output.gcsa} {input} 2> {log}
        rm -r {params.tmp_dir}
        """

# map reads to the graph using giraffe
rule map_giraffe:
    input:
        fastq=S3.remote(SROOT + '/reads/{read}.' + CHR + '.interleaved.fastq.gz'),
        xg=S3.remote(SROOT + '/{method}/{params}/map/{dataset}.{method}.{params}.xg'),
        min=S3.remote(SROOT + '/{method}/{params}/map/{dataset}.{method}.{params}.k{k}.w{w}.N{n}.min'),
        dist=S3.remote(SROOT + '/{method}/{params}/map/{dataset}.{method}.{params}.dist'),
        gbwt=S3.remote(SROOT + '/{method}/{params}/map/{dataset}.{method}.{params}.N{n}.gbwt')
    output: '{method}/{params}/map/{dataset}.{method}.{params}.{read}.giraffe{k}k{w}w{n}N.gaf'
    threads: config['map_cores']
    benchmark: S3.remote(SROOT + '/benchmarks/{dataset}.{method}.{params}.mapvg-giraffe{k}k{w}w{n}N-{read}.benchmark.tsv')
    log: S3.remote(SROOT + '/logs/{dataset}.{method}.{params}.mapvg-giraffe{k}k{w}w{n}N-{read}.log')
    singularity:
        "docker://" + config['vg_docker']
    shell:
        'vg giraffe -o gaf -p -t {threads} -b fast -m {input.min} -d {input.dist} --gbwt-name {input.gbwt} -x {input.xg} -N {wildcards.read} -i -f {input.fastq} > {output} 2> {log}'

rule map_giraffe_paths:
    input:
        fastq=S3.remote(SROOT + '/reads/{read}.' + CHR + '.interleaved.fastq.gz'),
        xg=S3.remote(SROOT + '/{method}/{params}/map/{dataset}.{method}.{params}.xg'),
        min=S3.remote(SROOT + '/{method}/{params}/map/{dataset}.{method}.{params}.k{k}.w{w}.paths.min'),
        dist=S3.remote(SROOT + '/{method}/{params}/map/{dataset}.{method}.{params}.dist'),
        gbwt=S3.remote(SROOT + '/{method}/{params}/map/{dataset}.{method}.{params}.paths.gbwt')
    output: '{method}/{params}/map/{dataset}.{method}.{params}.{read}.giraffe{k}k{w}wPaths.gaf'
    threads: config['map_cores']
    benchmark: S3.remote(SROOT + '/benchmarks/{dataset}.{method}.{params}.mapvg-giraffe{k}k{w}wPaths-{read}.benchmark.tsv')
    log: S3.remote(SROOT + '/logs/{dataset}.{method}.{params}.mapvg-giraffe{k}k{w}wPaths-{read}.log')
    singularity:
        "docker://" + config['vg_docker']
    shell:
        'vg giraffe -o gaf -p -t {threads} -b fast -m {input.min} -d {input.dist} --gbwt-name {input.gbwt} -x {input.xg} -N {wildcards.read} -i -f {input.fastq} > {output} 2> {log}'

# map reads to the graph using map
rule map_vgmap:
    input:
        fastq=S3.remote(SROOT + '/reads/{read}.' + CHR + '.interleaved.fastq.gz'),
        xg=S3.remote(SROOT + '/{method}/{params}/map/{dataset}.{method}.{params}.xg'),
        gcsa=S3.remote(SROOT + '/{method}/{params}/map/{dataset}.{method}.{params}.gcsa'),
        gcsa_lcp=S3.remote(SROOT + '/{method}/{params}/map/{dataset}.{method}.{params}.gcsa.lcp')
    output: '{method}/{params}/map/{dataset}.{method}.{params}.{read}.map.gaf'
    threads: config['map_cores']
    benchmark: S3.remote(SROOT + '/benchmarks/{dataset}.{method}.{params}.mapvg-map-{read}.benchmark.tsv')
    log: S3.remote(SROOT + '/logs/{dataset}.{method}.{params}.mapvg-map-{read}.log')
    singularity:
        "docker://" + config['vg_docker']
    shell:
        'vg map --gaf -t {threads} -x {input.xg} -g {input.gcsa} -i -f {input.fastq} -N {wildcards.read} > {output} 2> {log}'

# pack coverage from the reads
rule pack_ga:
    input:
        gaf=S3.remote(SROOT + '/{method}/{params}/map/{dataset}.{method}.{params}.{read}.ga.gaf.gz'),
        xg=S3.remote(SROOT + '/{method}/{params}/map/{dataset}.{method}.{params}.xg')
    output: S3.remote(SROOT + '/{method}/{params}/map/{dataset}.{method}.{params}.{read}.ga.pack')
    threads: 4
    benchmark: S3.remote(SROOT + '/benchmarks/{dataset}.{method}.{params}.mapvg-ga-pack-{read}.benchmark.tsv')
    log: S3.remote(SROOT + '/logs/{dataset}.{method}.{params}.mapvg-ga-pack-{read}.log')
    singularity:
        "docker://" + config['vg_docker']
    shell:
        "vg pack -x {input.xg} -a {input.gaf} -Q 0 -t {threads} -o {output} 2> {log}"

# pack coverage from the reads
rule pack:
    input:
        gaf=S3.remote(SROOT + '/{method}/{params}/map/{dataset}.{method}.{params}.{read}.{mapper}.gaf.gz'),
        xg=S3.remote(SROOT + '/{method}/{params}/map/{dataset}.{method}.{params}.xg')
    output: S3.remote(SROOT + '/{method}/{params}/map/{dataset}.{method}.{params}.{read}.{mapper}.pack')
    threads: 4
    benchmark: S3.remote(SROOT + '/benchmarks/{dataset}.{method}.{params}.mapvg-{mapper}-pack-{read}.benchmark.tsv')
    log: S3.remote(SROOT + '/logs/{dataset}.{method}.{params}.mapvg-{mapper}-pack-{read}.log')
    singularity:
        "docker://" + config['vg_docker']
    shell:
        "vg pack -x {input.xg} -a {input.gaf} -Q 5 -t {threads} -o {output} 2> {log}"

# call variants from the packed read coverage
rule call_novcf:
    input:
        pack=S3.remote(SROOT + '/{method}/{params}/map/{dataset}.{method}.{params}.{read}.{mapper}.pack'),
        xg=S3.remote(SROOT + '/{method}/{params}/map/{dataset}.{method}.{params}.xg'),
        snarls=S3.remote(SROOT + '/{method}/{params}/map/{dataset}.{method}.{params}.snarls')
    output: '{method}/{params}/map/{dataset}.{method}.{params}.{read}.{mapper}.call.vcf'
    threads: 8
    benchmark: S3.remote(SROOT + '/benchmarks/{dataset}.{method}.{params}.mapvg-{mapper}-call-{read}.benchmark.tsv')
    log: S3.remote(SROOT + '/logs/{dataset}.{method}.{params}.mapvg-{mapper}-call-{read}.log')
    singularity:
        "docker://" + config['vg_docker']
    shell:
        """
        CHRP=`vg paths -L -x {input.xg} | grep -e hg38 -e "^{CHR}" | head -1`
        vg call -k {input.pack} -t {threads} -s HG002 --snarls {input.snarls} -p $CHRP {input.xg} > {output} 2> {log}
        """

##
## Map reads using GraphAligner
##

rule map_graphaligner:
    input:
        fastq=S3.remote(SROOT + '/reads/{read}.' + CHR + '.fastq.gz'),
        vg=S3.remote(SROOT + '/{method}/{params}/{dataset}.{method}.{params}.vg')
    output: '{method}/{params}/map/{dataset}.{method}.{params}.{read}.ga.gaf'
    threads: config['map_cores']
    benchmark: S3.remote(SROOT + '/benchmarks/{dataset}.{method}.{params}.graphaligner-{read}.benchmark.tsv')
    log: S3.remote(SROOT + '/logs/{dataset}.{method}.{params}.graphaligner-{read}.log')
    params:
        gfa = 'temp.{dataset}.{method}.{params}.gfa'
    shell:
        """
        vg view -g {input.vg} > {params.gfa}
        GraphAligner -g {params.gfa} -f {input.fastq} -a {output} -x vg -t {threads} 2> {log}
        rm -f {params.gfa}
        """

# rule map_graphaligner:
#     input:
#         fastq=S3.remote(SROOT + '/reads/{read}.' + CHR + '.fastq.gz'),
#         vg=S3.remote(SROOT + '/{method}/{params}/{dataset}.{method}.{params}.vg')
#     output: '{method}/{params}/map/{dataset}.{method}.{params}.{read}.ga.gaf'
#     threads: config['map_cores']
#     benchmark: S3.remote(SROOT + '/benchmarks/{dataset}.{method}.{params}.graphaligner-{read}.benchmark.tsv')
#     log: S3.remote(SROOT + '/logs/{dataset}.{method}.{params}.graphaligner-{read}.log')
#     params:
#         gfa = 'temp.{dataset}.{method}.{params}.gfa',
#         pg = 'temp.{dataset}.{method}.{params}.pg'
#     shell:
#         """
#         vg convert -p {input.vg} > {params.pg}
#         vg view -g {params.pg} > {params.gfa}
#         GraphAligner -g {params.gfa} -f {input.fastq} -a {output} -x vg -t {threads} 2> {log}
#         rm {params.gfa} {params.pg}
#         """

# Preference to specific rules if multiple could be used
ruleorder: bgzip_fa > gzip
ruleorder: gfa_to_vg_minigraph > gfa_to_vg
ruleorder: seqwish_kl_gz > seqwish_kl
ruleorder: deconstruct_vcf_minigraph > deconstruct_vcf
ruleorder: deconstruct_vcf_paftools > bgzip_vcf
ruleorder: pack_ga > pack
