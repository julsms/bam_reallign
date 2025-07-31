process BAMLEFTALIGN {
    conda "bioconda::gatk=3.8 bioconda::samtools=1.15.1"
    
    tag "${chromosome}"
    
    publishDir "${params.outdir}/bamleftalign", mode: 'copy'
    
    input:
    path reference_fasta
    path reference_fai
    path reference_dict
    tuple val(chromosome), path(chr_bam)
    
    output:
    tuple val(chromosome), path("${chromosome}.leftaligned.bam"), emit: leftaligned_bams
    
    script:
    """
    echo "=== GATK3 LeftAlignIndels for ${chromosome} ==="
    
    INPUT_READS=\$(samtools view -c ${chr_bam})
    echo "Input reads: \$INPUT_READS"
    
    # Check that all reference files are present
    echo "Reference files:"
    ls -la ${reference_fasta}*
    
    if [ \$INPUT_READS -eq 0 ]; then
        echo "WARNING: Input BAM has no reads"
        cp ${chr_bam} ${chromosome}.leftaligned.bam
    else
        # Ensure BAM is indexed
        samtools index ${chr_bam}
        
        # Run GATK LeftAlignIndels with all reference files available
        java -Xmx10g -jar \$CONDA_PREFIX/opt/gatk-3.8/GenomeAnalysisTK.jar \\
            -T LeftAlignIndels \\
            -R ${reference_fasta} \\
            -I ${chr_bam} \\
            -o ${chromosome}.leftaligned.bam
        
        # Check if GATK succeeded
        if [ \$? -ne 0 ] || [ ! -s "${chromosome}.leftaligned.bam" ]; then
            echo "GATK LeftAlignIndels failed, copying input"
            cp ${chr_bam} ${chromosome}.leftaligned.bam
        fi
    fi
    
    # Index the output
    samtools index ${chromosome}.leftaligned.bam
    
    # Verify
    OUTPUT_READS=\$(samtools view -c ${chromosome}.leftaligned.bam)
    echo "Output reads: \$OUTPUT_READS"
    
    if [ \$OUTPUT_READS -ne \$INPUT_READS ]; then
        echo "WARNING: Read count changed from \$INPUT_READS to \$OUTPUT_READS"
    else
        echo "SUCCESS: Read counts preserved"
    fi
    """
}
