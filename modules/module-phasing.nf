#!/usr/bin/env nextflow

paneldir       = file(params.paneldir, type: 'dir')
map            = file(params.map, type: 'file')
cnps           = file(params.cnps, type: 'file')
input_stats    = file(params.input_stats, type: 'file')
outdir         = file(params.outdir, type: 'dir')


// Extract genotypes and split by autosomes and chromosome X
process run_ExtractGenotypes {
    tag { "extract_genotypes" }
    publishDir "${outdir}/extracted_genotypes", mode: 'copy', overwrite: false

    input:
    tuple val(unphased), path(bcf_unphased), path(index_unphased)
    tuple val(exclude), path(bcf_exclude), path(index_exclude)
       
    output:
    path("mocha_chr*.bcf"), emit: chr_genotypes
    
    """
    bcftools isec --no-version --threads ${task.cpus} -Ou --complement --exclude "N_ALT>1" --write 1 ${bcf_unphased} ${bcf_exclude} | \
        bcftools annotate --no-version --threads ${task.cpus}  -Ou --remove ID,QUAL,INFO,^FMT/GT  | \
        bcftools +scatter --no-version -Ob --output . --scatter \$(echo {1..22} | tr ' ' ',') --prefix mocha_chr
    """
}

// Phase VCF file by chromosome with SHAPEIT4
process run_PhaseBCF {
    tag { "phasing: ${chr_num}" }
    maxForks 10
    memory '100 GB'
    clusterOptions = '--constraint=avx2'
    publishDir "${outdir}/phased_genotypes", mode: 'copy', overwrite: false

    input:
    tuple val(chr_num), path(chr)

    output:
    path("mocha_${chr_num}.pgt{.bcf,.log}"), emit: pgt
    path("genetic_map.${chr_num}.txt"), emit: gen_map
    
    """
    bcftools index --force ${chr}
    zcat ${map} | sed 's/^23/X/' | awk -v chr=${chr_num[3..-1]} '\$1==chr {print \$2,\$3,\$4}' > genetic_map.${chr_num}.txt
    /opt/shapeit4-4.2.1/bin/shapeit4.2 --thread ${task.cpus} \
        --input ${chr} \
        --reference ${paneldir}/ALL.${chr_num}.phase3_integrated.20130502.genotypes.bcf \
        --map genetic_map.${chr_num}.txt \
        --region ${chr_num[3..-1]} \
        --output mocha_${chr_num}.pgt.bcf \
        --log mocha_${chr_num}.pgt.log
   """
}

// Concatenate phased output into a single VCF file
process run_MergeBCFs {
    tag { "merge pgt" }
    publishDir "${outdir}/phased_genotypes", mode: 'copy', overwrite: false

    input:
    path(pgt)

    output:
    path("mocha_pgt{.bcf,.bcf.csi}"), emit: comb_pgt

    """
    bcftools concat --no-version --threads ${task.cpus} -Ob mocha_chr{1..22}.pgt.bcf | \
        tee mocha_pgt.bcf | \
        bcftools index --threads ${task.cpus} --force --output mocha_pgt.bcf.csi
    """
}

// Import phased genotypes in the original VCF without changing missing genotypes
process run_ImportPhased {
    tag { "phasing: ${chr_num}" }
    publishDir "${outdir}/phased_genotypes", mode: 'copy', overwrite: false
    
    input:
    path(comb_pgt)
    tuple val(unphased), path(bcf_unphased), path(index_unphased)
    
    output:
    path("mocha_pgt.pg{.bcf,.bcf.csi}"), emit: phased_bcf
    
    """
    bcftools annotate --no-version --threads ${task.cpus} -Ob --annotations mocha_pgt.bcf --columns -FMT/GT ${bcf_unphased} | \
        tee mocha_pgt.pg.bcf | \
        bcftools index --force --threads ${task.cpus} --output mocha_pgt.pg.bcf.csi
    """
}

// Chromosomal alterations
process run_ChromAlt {
    tag { "chrom_alt" }
    memory '100 GB'
    publishDir "${outdir}/chromosomal_alterations", mode: 'copy', overwrite: false

    input:
    path(phased_bcf)
    tuple val(exclude), path(bcf_exclude), path(index_exclude)
    
    output:
    path("mocha_pgt.pg.as*"), emit: chrom_alt
    
    """
    bcftools +mocha \
        --genome GRCh37 \
        --input-stats ${input_stats} \
        --no-version \
        --threads ${task.cpus}
        --output - \
        --output-type b \
        --variants ^${bcf_exclude} \
        --calls mocha_pgt.pg.as.calls.tsv \
        --stats mocha_pgt.pg.as.stats.tsv \
        --ucsc-bed mocha_pgt.pg.as.ucsc.bed \
        --cnp ${cnp} \
        --mhc 6:27486711-33448264 --kir 19:54574747-55504099 \
        mocha_pgt.pg.bcf | \
        tee mocha_pgt.pg.as.bcf | \
        bcftools index --force --threads ${task.cpus} --output mocha_pgt.pg.as.bcf.csi
    """
}


