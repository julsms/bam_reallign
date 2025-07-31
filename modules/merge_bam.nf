process MERGE_BAM {
    conda "bioconda::samtools=${params.samtools_version}"
    
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path realigned_bams
    
    output:
    path "final.realigned.bam", emit: final_bam
    path "final.realigned.bam.bai", emit: final_bai
    
    script:
    """
    # Merge realigned per-chromosome BAMs into one file
    samtools merge -@ ${task.cpus} final.realigned.bam ${realigned_bams}
    
    # Index the final merged BAM
    samtools index final.realigned.bam
    """
}
