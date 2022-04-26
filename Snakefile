#!/usr/bin/env python3

###########
# GLOBALS #
###########

#containers
vcftools_container = 'docker://biocontainers/vcftools:v0.1.16-1-deb_cv1'
bioconductor_container = 'library://sinwood/bioconductor/bioconductor_3.12:0.0.1'
plink_v2_container = 'library://sinwood/plink/plink_2.0:0.0.1'

all_gwas = ["covars/all_samples_pruned/all_samples_pruned.attack_status",
            "covars/all_samples_pruned/all_samples_pruned.location",
            "covars/all_samples_pruned/all_samples_pruned.parasitism",
            "covars/all_samples/all_samples.attack_status",
            "covars/all_samples/all_samples.location",
            "covars/all_samples/all_samples.parasitism",
            "covars/dunedin/dunedin.attack_status",
            "covars/dunedin/dunedin.parasitism",
            "covars/ruakura/ruakura.attack_status",
            "no_covars/all_samples_pruned/all_samples_pruned.attack_status",
            "no_covars/all_samples_pruned/all_samples_pruned.location",
            "no_covars/all_samples_pruned/all_samples_pruned.parasitism",
            "no_covars/all_samples/all_samples.attack_status",
            "no_covars/all_samples/all_samples.location",
            "no_covars/all_samples/all_samples.parasitism",
            "no_covars/dunedin/dunedin.attack_status",
            "no_covars/dunedin/dunedin.parasitism",
            "no_covars/ruakura/ruakura.attack_status", 
            "no_covars/ruakura/ruakura.parasitism"]

#########
# RULES #
#########

rule target:
    input:
        # plink2 - attack & para
        expand('output/plink2/glm_analysis/covars/{vcf}/{vcf}.attack_status.glm.logistic.adjusted',  vcf=['all_samples', 'all_samples_pruned', 'dunedin', 'ruakura']),
        expand('output/plink2/glm_analysis/no_covars/{vcf}/{vcf}.attack_status.glm.logistic.adjusted',  vcf=['all_samples', 'all_samples_pruned', 'dunedin', 'ruakura']),
        # plink2 - location
        expand('output/plink2/glm_analysis/covars/{vcf}/{vcf}.location.glm.logistic.adjusted', vcf=['all_samples', 'all_samples_pruned']),
        expand('output/plink2/glm_analysis/no_covars/{vcf}/{vcf}.location.glm.logistic.adjusted', vcf=['all_samples', 'all_samples_pruned']),
        'output/plink2/glm_analysis/sig_associations.out'

## helpful resources
    # https://evomics.org/learning/population-and-speciation-genomics/2016-population-and-speciation-genomics/fileformats-vcftools-plink/#ex2.3
    # https://angus.readthedocs.io/en/2017/GWAS.html
    # https://plink.readthedocs.io/en/latest/GWAS/
    # https://yosuketanigawa.com/posts/2020/09/PLINK2
    # https://www.cog-genomics.org/plink/2.0/strat#pca

###############
## Plink 2.0 ##
###############

rule grep_plink2_results:
    input:
        expand('output/plink2/glm_analysis/{gwas}_assoc_analysis.log', gwas=all_gwas)
    output:
        'output/plink2/glm_analysis/sig_associations.out'
    shell:
        'grep "significant assoc.s after Bonferroni correction" {input} > {output}'

rule plink2_results_analysis:
    input:
        results = 'output/plink2/glm_analysis/{gwas}.glm.logistic.adjusted'
    output:
        qq_plot = 'output/plink2/glm_analysis/{gwas}_qq.pdf',
        manhattan_plot = 'output/plink2/glm_analysis/{gwas}_manhattan.pdf'
    log:
        'output/plink2/glm_analysis/{gwas}_assoc_analysis.log'
    singularity:
        bioconductor_container
    script:
        'src/assoc_analysis.R'

##### have to split location analysis from attack/para as can't do location analysis on ruakura/invermay only vcfs #####

rule plink_glm_logistic_covars_location:
    input:
        pgen = 'output/plink2/plink_input/{vcf}/{vcf}.pgen',
        pheno = 'data/plink2/location_numeric.txt',
        covars = 'output/plink2/plink_input/{vcf}/all_covars.tsv'
    output:
        'output/plink2/glm_analysis/covars/{vcf}/{vcf}.location.glm.logistic.adjusted'
    params:
        outdir = 'output/plink2/glm_analysis/covars/{vcf}/{vcf}',
        input_pref = 'output/plink2/plink_input/{vcf}/{vcf}'
    log:
        'output/logs/plink2/plink_glm_logistic_covars_location/plink_{vcf}.log'
    singularity:
        plink_v2_container
    threads:
        10
    shell:
        '/plink/plink2 '
        '--pfile {params.input_pref} '
        '--glm no-x-sex pheno-ids '
        '--allow-extra-chr '
        '--covar {input.covars} '
        '--pheno {input.pheno} '
        '--adjust '
        '--ci 0.95 '
        '--out {params.outdir} '
        '--threads {threads} '
        '2> {log}'

rule plink_glm_logistic_location:
    input:
        pgen = 'output/plink2/plink_input/{vcf}/{vcf}.pgen',
        pheno = 'data/plink2/location_numeric.txt',
    output:
        'output/plink2/glm_analysis/no_covars/{vcf}/{vcf}.location.glm.logistic.adjusted'
    params:
        outdir = 'output/plink2/glm_analysis/no_covars/{vcf}/{vcf}',
        input_pref = 'output/plink2/plink_input/{vcf}/{vcf}'
    log:
        'output/logs/plink2/{vcf}/plink_glm_logistic_location.log'
    singularity:
        plink_v2_container
    threads:
        10
    shell:
        '/plink/plink2 '
        '--pfile {params.input_pref} '
        '--glm no-x-sex pheno-ids '
        '--allow-extra-chr '
        '--pheno {input.pheno} '
        '--adjust '
        '--ci 0.95 '
        '--out {params.outdir} '
        '--threads {threads} '
        '2> {log}'

#### Attack & Para analyses ####

rule plink_glm_logistic_covars:
    input:
        pgen = 'output/plink2/plink_input/{vcf}/{vcf}.pgen',
        pheno = 'data/plink2/phenotypes_numeric.txt',
        covars = 'output/plink2/plink_input/{vcf}/all_covars.tsv'
    output:
        attack = 'output/plink2/glm_analysis/covars/{vcf}/{vcf}.attack_status.glm.logistic.adjusted' # skips para for ruakura-only as more covar columns than para samples
    params:
        outdir = 'output/plink2/glm_analysis/covars/{vcf}/{vcf}',
        input_pref = 'output/plink2/plink_input/{vcf}/{vcf}'
    log:
        'output/logs/plink2/plink_glm_logistic_covars/plink_{vcf}.log'
    singularity:
        plink_v2_container
    threads:
        10
    shell:
        '/plink/plink2 '
        '--pfile {params.input_pref} '
        '--glm no-x-sex pheno-ids '
        '--allow-extra-chr '
        '--covar {input.covars} '
        '--pheno {input.pheno} '
        '--adjust '
        '--ci 0.95 '
        '--out {params.outdir} '
        '--threads {threads} '
        '2> {log}'

rule plink_glm_logistic:
    input:
        pgen = 'output/plink2/plink_input/{vcf}/{vcf}.pgen',
        pheno = 'data/plink2/phenotypes_numeric.txt', # plink2 can have categorical covars but NOT categorical phenotypes
    output:
        attack = 'output/plink2/glm_analysis/no_covars/{vcf}/{vcf}.attack_status.glm.logistic.adjusted'
    params:
        outdir = 'output/plink2/glm_analysis/no_covars/{vcf}/{vcf}',
        input_pref = 'output/plink2/plink_input/{vcf}/{vcf}'
    log:
        'output/logs/plink2/{vcf}/plink_glm_logistic.log'
    singularity:
        plink_v2_container
    threads:
        10
    shell:
        '/plink/plink2 '
        '--pfile {params.input_pref} '
        '--glm no-x-sex pheno-ids ' # don't use sex, & write list of samples used in each regression
        '--allow-extra-chr '
        '--pheno {input.pheno} '
        '--adjust '
        '--ci 0.95 '
        '--out {params.outdir} '
        '--threads {threads} '
        '2> {log}'

#### plink2 set up ####

rule merge_covars:
    input:
        pca_covars_file = 'output/plink2/plink_input/{vcf}/{vcf}.eigenvec',
        expt_covars_file = 'data/plink2/covariates.tsv'
    output:
        merged_covars = 'output/plink2/plink_input/{vcf}/all_covars.tsv'
    singularity:
        bioconductor_container
    log:
        'output/plink2/{vcf}/merge_covars.log'
    script:
        'src/merge_covars.R'

rule plink2_pca:
    input:
        'output/plink2/plink_input/{vcf}/{vcf}.rel.id'
    output:
        'output/plink2/plink_input/{vcf}/{vcf}.eigenvec'
    params:
        prefix = 'output/plink2/plink_input/{vcf}/{vcf}',
    log:
        'output/logs/plink2/{vcf}/plink2_pca_for_gwas.log'
    singularity:
        plink_v2_container
    shell:
        '/plink/plink2 '
        '--pfile {params.prefix} '
        '--pca '
        '--allow-extra-chr '
        '-out {params.prefix} '
        '2> {log}'

rule plink2_make_rel:
    input:
        'output/plink2/plink_input/{vcf}/{vcf}.pgen'
    output:
        'output/plink2/plink_input/{vcf}/{vcf}.rel.id'
    params:
        prefix = 'output/plink2/plink_input/{vcf}/{vcf}'
    log:
        'output/logs/plink2/{vcf}/plink2_make_rel.log'
    singularity:
        plink_v2_container
    shell:
        '/plink/plink2 '
        '--pfile {params.prefix} '
        '--make-rel '
        '--allow-extra-chr '
        '-out {params.prefix} '
        '2> {log}'

rule convert_vcf_to_plink2:
    input:
        vcf = 'data/asw_vcfs/{vcf}.vcf.gz'
    output:
        pgen = 'output/plink2/plink_input/{vcf}/{vcf}.pgen'
    params:
        prefix = 'output/plink2/plink_input/{vcf}/{vcf}'
    log:
        'output/logs/plink2/{vcf}/convertvcf_to_plink2.log'
    singularity:
        plink_v2_container
    shell:
        '/plink/plink2 '
        '--vcf {input.vcf} '
        '--allow-extra-chr '
        '--set-missing-var-ids @:# '
        '--make-pgen '
        '--out {params.prefix} '
        '2> {log}'

