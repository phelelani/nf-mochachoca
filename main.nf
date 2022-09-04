#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// iDAT FILES
Channel.fromFilePairs([
        '/dataE/AWIGenGWAS/Batch1-3/Batch1/Batch1_iDATS/**/*{Red,Grn}.idat',
        '/dataE/AWIGenGWAS/Batch1-3/Batch2/Batch2_iDATS/**/*{Red,Grn}.idat',
        '/dataE/AWIGenGWAS/Batch1-3/Batch3/Batch3_iDATS/**/*{Red,Grn}.idat',
        '/dataE/AWIGenGWAS/Batch4/Batch4_iDATS/**/*{Red,Grn}.idat'      
    ])
    .map { it -> [ it[0], it[1][0], it[1][1] ] }
    .set { idats }

Channel.fromFilePairs('/home/phelelani/projects/jonathan/RUN_2022_05_15/baf_lrr_passed_unphased{.bcf,.bcf.csi}', size: 2)
    .map { it -> [ it[0], it[1][0], it[1][1] ] }.set { unphased }
Channel.fromFilePairs('/home/phelelani/projects/jonathan/RUN_2022_05_15/baf_lrr_passed_unphased.xcl{.bcf,.bcf.csi}', size: 2)
    .map { it -> [ it[0], it[1][0], it[1][1] ] }.set { exclude }

// REQUIRED FILES
params.mode           = '' 
params.outdir         = '/home/phelelani/projects/jonathan/RUN_MOCHA_11_09_2022'
params.ref            = '/home/phelelani/nf-workflows/nf-mochachoca/data/GRCh37/human_g1k_v37.fasta'
params.bpm_manifest   = '/dataE/AWIGenGWAS/aux/H3Africa_2017_20021485_A3.bpm'
params.csv_manifest   = '/dataE/AWIGenGWAS/aux/H3Africa_2017_20021485_A3.csv'
params.egt_cluster    = '/home/phelelani/nf-workflows/nf-mochachoca/data/GenomeStudio-H3Africa-array-clusters-HapMap2-186-samples.egt'
params.passed_snps    = '/home/phelelani/nf-workflows/nf-mochachoca/data/awigen/passed_snps.list'
params.passed_samples = '/home/phelelani/nf-workflows/nf-mochachoca/data/awigen/passed_samples.list'
params.awigen_ids     = '/home/phelelani/nf-workflows/nf-mochachoca/data/awigen/h3achip_dbsnp150.tsv'
params.dup            = '/home/phelelani/nf-workflows/nf-mochachoca/data/GRCh37/segdups.bed.gz'
params.map            = '/home/phelelani/nf-workflows/nf-mochachoca/data/GRCh37/genetic_map_hg19_withX.txt.gz'
params.paneldir       = '/home/phelelani/nf-workflows/nf-mochachoca/data/GRCh37'
params.cnps           = '/home/phelelani/nf-workflows/nf-mochachoca/data/GRCh37/cnps.bed'
params.input_stats    = '/home/phelelani/nf-workflows/nf-mochachoca/data/GRCh37/input_stats.tsv'

// OUTPUT DIR
outdir = file(params.outdir, type: 'dir')
outdir.mkdir()

include { run_iDAT2GTC; run_GTC2BCF; run_FilterBCF; run_RenameIDs; run_ExcludeVariants } from './modules/module-idat2bcf.nf'
include { run_ExtractGenotypes; run_PhaseBCF; run_MergeBCFs; run_ImportPhased; run_ChromAlt } from './modules/module-phasing.nf'

// CONVERT iDAT FILE PAIRS TO GTC FILES 
workflow RUN_IDAT2BCF {
    take:
        idats
    main:
        run_iDAT2GTC(idats)
        run_GTC2BCF(run_iDAT2GTC.out.gtcs.collect())
        run_FilterBCF(run_GTC2BCF.out.bcf)
        run_RenameIDs(run_FilterBCF.out.passed_bcf)
        run_ExcludeVariants(run_RenameIDs.out.passed_unphased_bcf, run_GTC2BCF.out.tsv)
}

// CONVERT GTC FILES TO BCF
workflow RUN_PHASING {
    take:
        unphased
        exclude
    main:
        run_ExtractGenotypes(unphased, exclude)
        run_PhaseBCF(run_ExtractGenotypes.out.chr_genotypes.flatten().map { it -> [ it.baseName[6..-1], it ] })
        run_MergeBCFs(run_PhaseBCF.out.pgt.collect())
        run_ImportPhased(run_MergeBCFs.out.comb_pgt, unphased)
    run_ChromAlt(exclude, run_ImportPhased.out.phased_bcf).view()
}

// PICK AND CHOOSE 
workflow {
    mode = params.mode
    switch (mode) {
        case['idat2bcf']:
            RUN_IDAT2BCF(idats)
            break
            // =====
        case['phasing']:
            RUN_PHASING(unphased, exclude)
            break
            // =====
        default:
            exit 1, """
OOOPS!! SEEMS LIE WE HAVE A WORFLOW ERROR!
No workflow \'mode\' give! Please use one of the following options for workflows:
    --mode idat2bcf          // To run the 
    --mode phasing           // To run the 
"""
            break
            // =====}
    }
}
