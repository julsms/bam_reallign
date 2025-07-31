process SPLIT_BAM {
    conda "bioconda::samtools=${params.samtools_version}"
    
    publishDir "${params.outdir}/per_chr", mode: 'copy'
    
    input:
    path input_bam
    path input_bai
    
    output:
    path "*.bam", emit: chr_bams  // Simplified output - just the BAM files
    
    script:
    """
    # Get chromosome list from BAM index
    samtools idxstats ${input_bam} | cut -f1 | grep -v '*' > chromosomes.txt
    
    # Split BAM by chromosome
    while read CHR; do
        echo "Extracting chromosome: \$CHR"
        samtools view -@ 80 -b ${input_bam} "\$CHR" > "\${CHR}.bam"
        samtools index "\${CHR}.bam"
    done < chromosomes.txt
    """
}
