process ABRA_REALIGNMENT {
    conda "bioconda::abra2=2.23"
    
    tag "${chromosome}"
    
    publishDir "${params.outdir}/abra", mode: 'copy'
    
    input:
    path reference_fasta
    tuple val(chromosome), path(leftaligned_bam), path(expanded_bed)
    
    output:
    tuple val(chromosome), path("${chromosome}.realigned.bam"), emit: realigned_bams
    
    script:
    """
    echo "=== ABRA_REALIGNMENT for ${chromosome} ==="
    
    # Check if we have target regions
    if [ -s ${expanded_bed} ]; then
        echo "Running ABRA2 for ${chromosome}..."
        
        # Create a chromosome-specific reference (extract only this chromosome)
        samtools faidx ${reference_fasta} ${chromosome} > ${chromosome}.fa
        samtools faidx ${chromosome}.fa
        
        # Filter BED to only include current chromosome
        grep "^${chromosome}\t" ${expanded_bed} > ${chromosome}_targets.bed || touch ${chromosome}_targets.bed
        
        if [ -s ${chromosome}_targets.bed ]; then
            # Run ABRA with chromosome-specific reference
            abra2 --in ${leftaligned_bam} \\
                  --out ${chromosome}.realigned.bam \\
                  --ref ${chromosome}.fa \\
                  --targets ${chromosome}_targets.bed \\
                  --threads ${task.cpus} \\
                  --tmpdir ./abra_tmp
            
            # Check if ABRA succeeded
            if [ \$? -ne 0 ] || [ ! -f "${chromosome}.realigned.bam" ]; then
                echo "ABRA failed for ${chromosome}, copying input BAM"
                cp ${leftaligned_bam} ${chromosome}.realigned.bam
            fi
        else
            echo "No targets for ${chromosome}, copying input BAM"
            cp ${leftaligned_bam} ${chromosome}.realigned.bam
        fi
        
        # Clean up
        rm -rf ./abra_tmp
        rm -f ${chromosome}.fa ${chromosome}.fa.fai ${chromosome}_targets.bed
    else
        echo "No target regions for ${chromosome}, copying input BAM"
        cp ${leftaligned_bam} ${chromosome}.realigned.bam
    fi
    
    # Index the realigned BAM
    samtools index ${chromosome}.realigned.bam
    """
}
