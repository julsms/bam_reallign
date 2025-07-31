process ADD_READ_GROUPS {
    conda "bioconda::picard=2.23.0"
    
    tag "${chromosome}"
    
    publishDir "${params.outdir}/add_read_groups", mode: 'copy'
    
    input:
    tuple val(chromosome), path(chr_bam)
    path reference_fasta
    path reference_fai
    path reference_dict
    
    output:
    tuple val(chromosome), path("${chromosome}.rg.bam"), path("${chromosome}.rg.bam.bai"), emit: rg_bams
    
    script:
    """
    echo "=== ADD_READ_GROUPS for ${chromosome} ==="
    
    INPUT_READS=\$(samtools view -c ${chr_bam})
    echo "Input reads: \$INPUT_READS"
    
    if [ \$INPUT_READS -eq 0 ]; then
        echo "WARNING: Input BAM has no reads"
        cp ${chr_bam} ${chromosome}.rg.bam
    else
        # Add read groups directly (let Picard handle sequence dictionary)
        java -Xmx10g -jar \$CONDA_PREFIX/share/picard*/picard.jar AddOrReplaceReadGroups \\
            INPUT=${chr_bam} \\
            OUTPUT=${chromosome}.rg.bam \\
            RGID=${chromosome} \\
            RGLB=lib1 \\
            RGPL=illumina \\
            RGPU=unit1 \\
            RGSM=sample1 \\
            VALIDATION_STRINGENCY=LENIENT \\
            CREATE_INDEX=false
    fi
    
    # Create index with samtools
    samtools index ${chromosome}.rg.bam
    
    # Verify output
    OUTPUT_READS=\$(samtools view -c ${chromosome}.rg.bam)
    echo "Output reads: \$OUTPUT_READS"
    """
}
