#!/usr/bin/env python3

###########
# GLOBALS #
###########

#containers
plink_container = 'docker://biocontainers/plink1.9:v1.90b6.6-181012-1-deb_cv1'
vcftools_container = 'docker://biocontainers/vcftools:v0.1.16-1-deb_cv1'
bioconductor_container = 'library://sinwood/bioconductor/bioconductor_3.12:0.0.1'

#########
# RULES #
#########

rule target:
    input:
        expand('output/02_plink/parasitism_{vcf}/parasitism_{vcf}.assoc', vcf=['all_samples', 'all_samples_pruned', 'dunedin', 'ruakura']),
        expand('output/02_plink/attack_status_{vcf}/attack_status_{vcf}.assoc', vcf=['all_samples', 'all_samples_pruned', 'dunedin', 'ruakura']),
        expand('output/02_plink/location_{vcf}/location_{vcf}.assoc', vcf=['all_samples', 'all_samples_pruned']),
        #have to specify list without dunedin/ruakura_location as plotting crashes on them
        expand('output/02_plink/{analysis}/{analysis}_manhattan.pdf', analysis=['parasitism_all_samples', 'parasitism_all_samples_pruned', 'parasitism_dunedin', 'parasitism_ruakura', 'attack_status_all_samples', 'attack_status_all_samples_pruned', 'attack_status_dunedin', 'attack_status_ruakura', 'location_all_samples', 'location_all_samples_pruned']),
        'output/02_plink/sig_associations.out'

##helpful resources
#https://evomics.org/learning/population-and-speciation-genomics/2016-population-and-speciation-genomics/fileformats-vcftools-plink/#ex2.3
#https://angus.readthedocs.io/en/2017/GWAS.html

#--all-pheno    For basic association tests, loop through all phenotypes in --pheno file.
#--make-pheno   Generate a binary phenotype from a text table.
#--pheno    Load phenotype data from the given file instead of the main input fileset.

#############################
## 01 format vcf for plink ##
#############################

##edit .fam files to have second column of FID with sample name again
#then can load phenotype from alternate file using --pheno filename
##covariates file with FID and IID then columns for covariates, call file with --covar


#If a covariate file is also specified, then all covariates in that file will be included in the regression model, labelled COV1, COV2, etc.
#This is different to other commands which take only a single covariate (possibly working in conjunction with the --mcovar option).

## .fam ##
#--no-fid
#--no-parents
#--no-sex
#--no-pheno
#These allow you to use .fam or .ped files which lack family ID, parental ID, sex, and/or phenotype columns. - which I don't have

#plink gwas limited to 95 chromosomes - so don't know if it will work?
#don't use .fam files at all

rule grep_results:
    input:
        expand('output/logs/assoc_analysis/{analysis}.log', analysis=['parasitism_all_samples', 'parasitism_all_samples_pruned', 'parasitism_dunedin', 'parasitism_ruakura', 'attack_status_all_samples', 'attack_status_all_samples_pruned', 'attack_status_dunedin', 'attack_status_ruakura', 'location_all_samples', 'location_all_samples_pruned'])
    output:
        'output/02_plink/sig_associations.out'
    shell:
        'grep "significant assoc.s after Bonferroni correction" {input} > {output}'

# this isn't taking into account sample relatedness
    # can input PCAs as COVARs but would muck up location analysis?
rule assoc_analysis:
    input:
        assoc = 'output/02_plink/{analysis}/{analysis}.assoc'
    output:
        adjusted = 'output/02_plink/{analysis}/{analysis}.assoc.adjusted',
        qq_plot = 'output/02_plink/{analysis}/{analysis}_qq.pdf',
        manhattan_plot = 'output/02_plink/{analysis}/{analysis}_manhattan.pdf',
        log = 'output/logs/assoc_analysis/{analysis}.log'
    log:
        'output/logs/assoc_analysis/{analysis}.log'
    singularity:
        bioconductor_container
    script:
        'src/assoc_analysis.R'

##assoc test doesn't use covars
        #'--covar {input.covar} '
rule plink_assoc_location:
    input:
        ped = 'output/01_format_vcfs/{vcf}/{vcf}.ped',
        map = 'output/01_format_vcfs/{vcf}/{vcf}.map',
        pheno = 'data/phenotypes_numeric.txt'
    output:
        'output/02_plink/location_{vcf}/location_{vcf}.assoc'
    params:
        outdir = 'output/02_plink/location_{vcf}/location_{vcf}'
    log:
        'output/logs/plink_assoc_location/plink_{vcf}.log'
    singularity:
        plink_container
    shell:
        '/usr/bin/plink1.9 '
        '--ped {input.ped} '
        '--map {input.map} '
        '--allow-no-sex '
        '--pheno {input.pheno} '
        '--pheno-name location '
        '--allow-extra-chr '
        '--assoc '
        '--out {params.outdir} '
        '2> {log}'

rule plink_assoc_parasitism:
    input:
        ped = 'output/01_format_vcfs/{vcf}/{vcf}.ped',
        map = 'output/01_format_vcfs/{vcf}/{vcf}.map',
        pheno = 'data/phenotypes_numeric.txt'
    output:
        'output/02_plink/parasitism_{vcf}/parasitism_{vcf}.assoc'
    params:
        outdir = 'output/02_plink/parasitism_{vcf}/parasitism_{vcf}'
    log:
        'output/logs/plink_assoc_parasitism/plink_{vcf}.log'
    singularity:
        plink_container
    shell:
        '/usr/bin/plink1.9 '
        '--ped {input.ped} '
        '--map {input.map} '
        '--allow-no-sex '
        '--pheno {input.pheno} '
        '--pheno-name parasitism '
        '--allow-extra-chr '
        '--assoc '
        '--out {params.outdir} '
        '2> {log}'

rule plink_assoc_attack:
    input:
        ped = 'output/01_format_vcfs/{vcf}/{vcf}.ped',
        map = 'output/01_format_vcfs/{vcf}/{vcf}.map',
        pheno = 'data/phenotypes_numeric.txt'
    output:
        'output/02_plink/attack_status_{vcf}/attack_status_{vcf}.assoc'
    params:
        outdir = 'output/02_plink/attack_status_{vcf}/attack_status_{vcf}'
    log:
        'output/logs/plink_assoc_attack/plink_{vcf}.log'
    singularity:
        plink_container
    shell:
        '/usr/bin/plink1.9 '
        '--ped {input.ped} '
        '--map {input.map} '
        '--allow-no-sex '
        '--pheno {input.pheno} '
        '--pheno-name attack_status '
        '--allow-extra-chr '
        '--assoc '
        '--out {params.outdir} '
        '2> {log}'

##makes .ped and .map
rule vcftools_plink:
    input:
        vcf = 'data/asw_vcfs/{vcf}.vcf.gz'
    output:
        map = 'output/01_format_vcfs/{vcf}/{vcf}.map',
        ped = 'output/01_format_vcfs/{vcf}/{vcf}.ped'
    params:
        outdir = 'output/01_format_vcfs/{vcf}/{vcf}'
    log:
        'output/logs/vcftools_plink_{vcf}.log'
    singularity:
        vcftools_container
    shell:
        'vcftools '
        '--gzvcf {input.vcf} '
        '--plink '
        '--out {params.outdir} '
        '2> {log}'

# Hardy Weinburg Equil. - plink --hwe  -what effect does this have?
# relatedness - IBD for all sample pairs (using pruned SNPs) plink --genome
# multidimensional scaling for sample relatedness - plink --genome --min




# to make .ped and .map - but it didn't want .fam??
    # plink --recode 




#  plink wants .fam, .ped, and .map


#vcftools can maybe make map and ped using vcftools --vcf input --plink --out out_prefix
