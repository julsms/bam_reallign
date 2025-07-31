process EXPAND_REGIONS {
    conda "bioconda::bedtools=${params.bedtools_version}"
    
    tag "${chromosome}"

    publishDir "${params.outdir}/expanded_regions", mode: 'copy'
    
    input:
    path genome_file
    tuple val(chromosome), path(intervals), path(bed_file)
    
    output:
    tuple val(chromosome), path("${chromosome}.expanded.bed"), emit: expanded_regions
    
    script:
    """
    # Expand regions by specified number of base pairs
    if [ -s ${bed_file} ]; then
        bedtools slop -i ${bed_file} -g ${genome_file} -b ${params.expand_bp} > ${chromosome}.expanded.bed
    else
        # Create empty expanded BED file if input is empty
        touch ${chromosome}.expanded.bed
    fi
    """
}
