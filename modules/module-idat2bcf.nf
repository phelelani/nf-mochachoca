#!/usr/bin/env nextflow

ref            = file(params.ref, type: 'file')
bpm_manifest   = file(params.bpm_manifest, type: 'file')
csv_manifest   = file(params.csv_manifest, type: 'file')
egt_cluster    = file(params.egt_cluster, type: 'file')
passed_snps    = file(params.passed_snps, type: 'file')
passed_samples = file(params.passed_samples, type: 'file')
awigen_ids     = file(params.awigen_ids, type: 'file')
dup            = file(params.dup, type: 'file')
outdir         = file(params.outdir, type: 'dir')

// CONVERT iDAT FILE PAIRS TO GTC FILES
process run_iDAT2GTC  {
    tag { sample }
    errorStrategy 'ignore'

    input:
    tuple val(sample), path(red), path(green)

    output:
    path("*.gtc"), emit: gtcs
    
    """
    /bin/hostname
    LANG="en_US.UTF-8" /opt/iaap-cli/iaap-cli gencall \
        ${bpm_manifest} ${egt_cluster} . \
        --idat-folder . \
        --output-gtc \
        --num-threads 6 \
        --gender-estimate-call-rate-threshold -0.1
    """
}

// CONVERT GTC FILES TO VCF FILES
process run_GTC2BCF  {
    tag { 'gtc2vcf' }
    cpus 12
    publishDir "${outdir}", mode: 'copy', overwrite: false

    input:
    path(gtcs)

    output:
    path('baf_lrr.tsv'), emit: tsv
    tuple val('bcf_pair'), path("baf_lrr.bcf"), path("baf_lrr.bcf.csi"), emit: bcf
    
    """
    bcftools +gtc2vcf --no-version -Ou \
        --threads 12 \
        --bpm ${bpm_manifest} \
        --csv ${csv_manifest} \
        --egt ${egt_cluster} \
        --gtcs . \
        --fasta-ref ${ref} \
        --extra baf_lrr.tsv | \
        bcftools sort -Ou -T ./bcftools-sort.XXXXXX | \
        bcftools norm --no-version --threads 12 -Ob -c x -f ${ref} | \
        bcftools annotate --no-version --threads 12 -Ob -x ^INFO/ALLELE_A,^INFO/ALLELE_B,^INFO/GC,^FMT/GT,^FMT/BAF,^FMT/LRR | \
        tee baf_lrr.bcf | bcftools index --force --threads 12 --output baf_lrr.bcf.csi
    """
}

process run_FilterBCF {
    tag { "Filter BCF"}
    cpus 12
    publishDir "${outdir}", mode: 'copy', overwrite: false

    input:
    tuple val(bcf_pair), path(bcf), path(index)

    output:
    tuple val("passed_bcf_pair"), path("baf_lrr_passed.bcf"), path("baf_lrr_passed.bcf.csi"), emit: passed_bcf
    
    """
    bcftools view -R ${passed_snps} ${bcf} --threads 12 -Ou | \
        bcftools view -S ${passed_samples} --threads 12 --force-samples -Ob | \
        tee baf_lrr_passed.bcf | bcftools index --force --threads 12 --output baf_lrr_passed.bcf.csi
    """
}

process run_RenameIDs {
    tag { "rename_ids"}
    cpus 12
    publishDir "${out_dir}", mode: 'copy', overwrite: false

    input:
    tuple val(bcf_pair), path(bcf), path(index)

    output:
    tuple val(bcf_pair), path("baf_lrr_passed_unphased.bcf"), path("baf_lrr_passed_unphased.bcf.csi"), emit: passed_unphased_bcf
    
    """
    cut -f 1,2 ${awigen_ids} | sed 's/:/	/; 1d' > annotation.tab
    bgzip annotation.tab
    tabix -s1 -b2 -e2 annotation.tab.gz
    bcftools annotate -a annotation.tab.gz -c CHROM,POS,ID --threads ${task.cpus} -Ob ${bcf}| \
        tee baf_lrr_passed_unphased.bcf | bcftools index --force --threads ${task.cpus} --output baf_lrr_passed_unphased.bcf.csi 
    """
}

// GENERATE A LIST OF VARIANTS THAT WILL BE EXCLUDED FROM MODELING BY BOTH EAGLE AND MOCHA 
process run_ExcludeVariants {
    tag { "exclude_varianrs"}
    cpus 12
    publishDir "${out_dir}", mode: 'copy', overwrite: false
    
    input:
    path(tsv)
    tuple val(bcf_pair), path(bcf), path(idx)

    output:
    path('call_rate.txt'), emit: call_rate
    path('samples_xcl_list.txt'), emit: sample_xcl
    tuple val(bcf_pair), path('baf_lrr_passed_unphased.xcl.bcf'), path('baf_lrr_passed_unphased.xcl.bcf.csi'), emit: excluded_variants_bcf

    """
    awk -F\"\\t\" '{ print \$1"\\t"\$21}' ${tsv} | sed 's/.gtc//' > call_rate.txt

    awk -F\"\\t\" '\$2<.97 {print \$1}' call_rate.txt > samples_xcl_list.txt; echo '##INFO=<ID=JK,Number=1,Type=Float,Description=\"Jukes Cantor\">' | \
        bcftools annotate --no-version --threads ${task.cpus} -Ou -a ${dup} -c CHROM,FROM,TO,JK -h /dev/stdin ${bcf} | \
        bcftools view --no-version --threads ${task.cpus} -Ou -S ^samples_xcl_list.txt --force-samples | \
        bcftools +fill-tags --no-version -Ou -t ^Y,MT,chrY,chrM -- -t ExcHet,F_MISSING | \
        bcftools view --no-version --threads ${task.cpus} -Ou -G | \
        bcftools annotate --no-version --threads ${task.cpus} -Ob \
            -i 'FILTER!=\".\" && FILTER!=\"PASS\" || INFO/JK<.02 || INFO/ExcHet<1e-6 || INFO/F_MISSING>1-.97' \
            -x ^INFO/JK,^INFO/ExcHet,^INFO/F_MISSING | \
        tee baf_lrr_passed_unphased.xcl.bcf | \
        bcftools index --force --threads ${task.cpus} --output baf_lrr_passed_unphased.xcl.bcf.csi
    """
}
