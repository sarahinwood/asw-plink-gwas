#!/usr/bin/env python3

###########
# GLOBALS #
###########

#containers
plink_container = 'docker://biocontainers/plink1.9:v1.90b6.6-181012-1-deb_cv1'
vcftools_container = 'docker://biocontainers/vcftools:v0.1.16-1-deb_cv1'
bioconductor_container = 'library://sinwood/bioconductor/bioconductor_3.12:0.0.1'
plink_v2_container = 'library://sinwood/plink/plink_2.0:0.0.1'

old_all_gwas = ["para_attack/all_samples/all_samples.parasitism",
            "para_attack/all_samples/all_samples.attack_status",
            "para_attack/all_samples_pruned/all_samples_pruned.parasitism",
            "para_attack/all_samples_pruned/all_samples_pruned.attack_status",
            "para_attack/dunedin/dunedin.parasitism",
            "para_attack/dunedin/dunedin.attack_status",
            "para_attack/ruakura/ruakura.parasitism",
            "para_attack/ruakura/ruakura.attack_status",
            "location/all_samples_pruned/all_samples_pruned.location",
            "location/all_samples/all_samples.location"]

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
        # assoc analysis
        #expand('output/02_plink_assoc/para_attack/{vcf}/{vcf}.parasitism.assoc.adjusted', vcf=['all_samples', 'all_samples_pruned', 'dunedin', 'ruakura']),
        #expand('output/02_plink_assoc/location/{vcf}/{vcf}.location.assoc.adjusted', vcf=['all_samples', 'all_samples_pruned']),
        'output/02_plink_assoc/sig_associations.out',
        # logistic analysis
        #expand('output/03_plink_logistic/para_attack/{vcf}/{vcf}.parasitism.assoc.logistic.adjusted',  vcf=['all_samples', 'all_samples_pruned', 'dunedin', 'ruakura']),
        #expand('output/03_plink_logistic/location/{vcf}/{vcf}.location.assoc.logistic.adjusted',  vcf=['all_samples', 'all_samples_pruned']),
        'output/03_plink_logistic/sig_associations.out',
        ##testing stuff
        'test_output/no_covars/all_samples_pruned.assoc.logistic',
        # plink2 - attack & para
        expand('output/plink2/glm_analysis/covars/{vcf}/{vcf}.attack_status.glm.logistic.adjusted',  vcf=['all_samples', 'all_samples_pruned', 'dunedin', 'ruakura']),
        expand('output/plink2/glm_analysis/no_covars/{vcf}/{vcf}.attack_status.glm.logistic.adjusted',  vcf=['all_samples', 'all_samples_pruned', 'dunedin', 'ruakura']),
        # plink2 - location
        expand('output/plink2/glm_analysis/covars/{vcf}/{vcf}.location.glm.logistic.adjusted', vcf=['all_samples', 'all_samples_pruned']),
        expand('output/plink2/glm_analysis/no_covars/{vcf}/{vcf}.location.glm.logistic.adjusted', vcf=['all_samples', 'all_samples_pruned']),
        'output/plink2/glm_analysis/sig_associations.out'

        #expand('output/02_plink_assoc/attack_status_{vcf}/attack_status_{vcf}.assoc', vcf=['all_samples', 'all_samples_pruned', 'dunedin', 'ruakura']),
        #expand('output/02_plink_assoc/location_{vcf}/location_{vcf}.assoc', vcf=['all_samples', 'all_samples_pruned']),
        #have to specify list without dunedin/ruakura_location as plotting crashes on them
        #expand('output/02_plink_assoc/{analysis}/{analysis}_manhattan.pdf', analysis=['parasitism_all_samples', 'parasitism_all_samples_pruned', 'parasitism_dunedin', 'parasitism_ruakura', 'attack_status_all_samples', 'attack_status_all_samples_pruned', 'attack_status_dunedin', 'attack_status_ruakura', 'location_all_samples', 'location_all_samples_pruned']),
        #'output/02_plink_assoc/sig_associations.out',
        #expand('output/03_plink_logistic/parasitism_{vcf}/parasitism_{vcf}.assoc.logistic', vcf=['all_samples'])

##helpful resources
#https://evomics.org/learning/population-and-speciation-genomics/2016-population-and-speciation-genomics/fileformats-vcftools-plink/#ex2.3
#https://angus.readthedocs.io/en/2017/GWAS.html

###############
## Plink 2.0 ##
###############

# https://plink.readthedocs.io/en/latest/GWAS/
# https://yosuketanigawa.com/posts/2020/09/PLINK2
# https://www.cog-genomics.org/plink/2.0/strat#pca

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

# i thought while playing around this had found the same para SNP again without covars
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

###############
## Plink 1.9 ##
###############

rule plink_testing_why_para_notsig:
    input:
        ped = 'output/01_format_vcfs/all_samples_pruned/all_samples_pruned.ped',
        map = 'output/01_format_vcfs/all_samples_pruned/all_samples_pruned.map',
        pheno = 'data/phenotypes_numeric.txt',
        #covars = 'data/covariates_numeric.tsv'
    output:
        parasitism = 'test_output/no_covars/all_samples_pruned.assoc.logistic'
    params:
        outdir = 'test_output/no_covars/all_samples_pruned'
    log:
        'test_output/no_covars/plink_all_samples_pruned.log'
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
        '--logistic --adjust '
        '--ci 0.95 '
        #'--covar {input.covars} '
        #'--covar-name COV1 '
        '--out {params.outdir} '
        '2> {log}'


# --assoc isn't taking into account sample relatedness, logistic can with covariates
    # can input PCAs as COVARs but would muck up location analysis?


#plink association tests byond basic are linear (quantitative pheno) and logistic (qualitative pheno)

#########################
## 03 - plink logistic ## can we include location numeric covariate from PCA/alike analysis?
#########################

rule grep_logistic_results:
    input:
        expand('output/03_plink_logistic/{gwas}_assoc_analysis.log', gwas=all_gwas)
    output:
        'output/03_plink_logistic/sig_associations.out'
    shell:
        'grep "significant assoc.s after Bonferroni correction" {input} > {output}'

rule logistic_analysis:
    input:
        results = 'output/03_plink_logistic/{gwas}.assoc.logistic.adjusted'
    output:
        qq_plot = 'output/03_plink_logistic/{gwas}_qq.pdf',
        manhattan_plot = 'output/03_plink_logistic/{gwas}_manhattan.pdf'
    log:
        'output/03_plink_logistic/{gwas}_assoc_analysis.log'
    singularity:
        bioconductor_container
    script:
        'src/assoc_analysis.R'

##assoc analysis with location
rule plink_logistic_location:
    input:
        ped = 'output/01_format_vcfs/{vcf}/{vcf}.ped',
        map = 'output/01_format_vcfs/{vcf}/{vcf}.map',
        pheno = 'data/location_numeric.txt',
        covars = 'data/covariates_numeric.tsv'
    output:
        location = 'output/03_plink_logistic/location/{vcf}/{vcf}.location.assoc.logistic.adjusted'
    params:
        outdir = 'output/03_plink_logistic/location/{vcf}/{vcf}'
    log:
        'output/logs/plink_assoc/plink_{vcf}.log'
    singularity:
        plink_container
    shell:
        '/usr/bin/plink1.9 '
        '--ped {input.ped} '
        '--map {input.map} '
        '--allow-no-sex '
        '--pheno {input.pheno} '
        '--allow-extra-chr '
        '--logistic --adjust '
        '--all-pheno '
        '--ci 0.95 '
        '--covar {input.covars} '
        '--out {params.outdir} '
        '2> {log}'

##phenotype and covariates must be numeric for plink 1.9
##logistic analysis with para/attack_status
rule plink_logistic_attack_para:
    input:
        ped = 'output/01_format_vcfs/{vcf}/{vcf}.ped',
        map = 'output/01_format_vcfs/{vcf}/{vcf}.map',
        pheno = 'data/phenotypes_numeric.txt',
        covars = 'data/covariates_numeric.tsv'
    output:
        parasitism = 'output/03_plink_logistic/para_attack/{vcf}/{vcf}.parasitism.assoc.logistic.adjusted',
        attac_status = 'output/03_plink_logistic/para_attack/{vcf}/{vcf}.attack_status.assoc.logistic.adjusted'
    params:
        outdir = 'output/03_plink_logistic/para_attack/{vcf}/{vcf}'
    log:
        'output/logs/plink_logistic/plink_{vcf}.log'
    singularity:
        plink_container
    shell:
        '/usr/bin/plink1.9 '
        '--ped {input.ped} '
        '--map {input.map} '
        '--allow-no-sex '
        '--pheno {input.pheno} '
        '--allow-extra-chr '
        '--logistic --adjust '
        '--all-pheno '
        '--ci 0.95 '
        '--covar {input.covars} '
        '--out {params.outdir} '
        '2> {log}'

######################
## 02 - plink assoc ## doesn't use covariates
######################

rule grep_assoc_results:
    input:
        expand('output/02_plink_assoc/{gwas}_assoc_analysis.log', gwas=all_gwas)
    output:
        'output/02_plink_assoc/sig_associations.out'
    shell:
        'grep "significant assoc.s after Bonferroni correction" {input} > {output}'

rule assoc_analysis:
    input:
        results = 'output/02_plink_assoc/{gwas}.assoc.adjusted'
    output:
        qq_plot = 'output/02_plink_assoc/{gwas}_qq_para.pdf',
        manhattan_plot = 'output/02_plink_assoc/{gwas}_manhattan_para.pdf'
    log:
        'output/02_plink_assoc/{gwas}_assoc_analysis.log'
    singularity:
        bioconductor_container
    script:
        'src/assoc_analysis.R'

##assoc analysis with location
rule plink_assoc_location:
    input:
        ped = 'output/01_format_vcfs/{vcf}/{vcf}.ped',
        map = 'output/01_format_vcfs/{vcf}/{vcf}.map',
        pheno = 'data/location_numeric.txt'
    output:
        location = 'output/02_plink_assoc/location/{vcf}/{vcf}.location.assoc.adjusted'
    params:
        outdir = 'output/02_plink_assoc/location/{vcf}/{vcf}'
    log:
        'output/logs/plink_assoc/plink_{vcf}.log'
    singularity:
        plink_container
    shell:
        '/usr/bin/plink1.9 '
        '--ped {input.ped} '
        '--map {input.map} '
        '--allow-no-sex '
        '--pheno {input.pheno} '
        '--all-pheno '
        '--allow-extra-chr '
        '--assoc --adjust '
        '--out {params.outdir} '
        '2> {log}'

##assoc analysis with para & attack_status
rule plink_assoc_attack_para:
    input:
        ped = 'output/01_format_vcfs/{vcf}/{vcf}.ped',
        map = 'output/01_format_vcfs/{vcf}/{vcf}.map',
        pheno = 'data/phenotypes_numeric.txt'
    output:
        parasitism = 'output/02_plink_assoc/para_attack/{vcf}/{vcf}.parasitism.assoc.adjusted',
        attack_status = 'output/02_plink_assoc/para_attack/{vcf}/{vcf}.attack_status.assoc.adjusted'
    params:
        outdir = 'output/02_plink_assoc/para_attack/{vcf}/{vcf}'
    log:
        'output/logs/plink_assoc/plink_{vcf}.log'
    singularity:
        plink_container
    shell:
        '/usr/bin/plink1.9 '
        '--ped {input.ped} '
        '--map {input.map} '
        '--allow-no-sex ' # need to tell plink to ignore no sex info/chromosomes
        '--pheno {input.pheno} '
        '--all-pheno ' # test all phenotypes
        '--allow-extra-chr ' # expects human number
        '--assoc --adjust ' # adjust = correct p values for multiple testing
        '--out {params.outdir} '
        '2> {log}'

#######################
## 01 - format files ##
#######################

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


##to include sample relatedness in gwas analysis can perform plink mds analysis to check for outliers and use main components as covariates
##‐‐ mperm can be used with assoc to do permutation to help with correcting for multiple testing

# Hardy Weinburg Equil. - plink --hwe  -what effect does this have?
# relatedness - IBD for all sample pairs (using pruned SNPs) plink --genome
# multidimensional scaling for sample relatedness - plink --genome --min
