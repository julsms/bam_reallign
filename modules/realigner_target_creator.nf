process REALIGNER_TARGET_CREATOR {
    conda "bioconda::gatk=3.8"
    
    tag "${chromosome}"
    
    publishDir "${params.outdir}/realigner_targets", mode: 'copy'
    
    input:
    path reference_fasta
    path reference_fai
    path reference_dict
    tuple val(chromosome), path(rg_bam)
    
    output:
    tuple val(chromosome), path("${chromosome}.intervals"), path("${chromosome}.bed"), emit: target_intervals
    
    script:
    """
    echo "=== REALIGNER_TARGET_CREATOR for ${chromosome} ==="
    
    # Find realignment targets with GATK
    java -Xmx10g -jar \$CONDA_PREFIX/opt/gatk-3.8/GenomeAnalysisTK.jar \\
        -T RealignerTargetCreator \\
        -R ${reference_fasta} \\
        -I ${rg_bam} \\
        -o ${chromosome}.intervals
    
    # Convert intervals to BED format with proper error handling
    if [ -s ${chromosome}.intervals ]; then
        echo "Converting intervals to BED format..."
        
        # Convert GATK intervals to BED format, handling single positions properly
        while read line; do
            if [[ \$line == *":"* ]]; then
                # Extract chromosome
                CHR=\$(echo \$line | cut -d':' -f1)
                
                # Extract coordinates
                COORDS=\$(echo \$line | cut -d':' -f2)
                
                if [[ \$COORDS == *"-"* ]]; then
                    # Range format: start-end
                    START=\$(echo \$COORDS | cut -d'-' -f1)
                    END=\$(echo \$COORDS | cut -d'-' -f2)
                    # Convert to 0-based for BED
                    START=\$((START - 1))
                    
                    # Only output if we have valid coordinates
                    if [[ \$START =~ ^[0-9]+\$ ]] && [[ \$END =~ ^[0-9]+\$ ]]; then
                        echo -e "\$CHR\\t\$START\\t\$END"
                    fi
                else
                    # Single position - make it a 1bp interval
                    if [[ \$COORDS =~ ^[0-9]+\$ ]]; then
                        START=\$((COORDS - 1))
                        END=\$COORDS
                        echo -e "\$CHR\\t\$START\\t\$END"
                    fi
                fi
            fi
        done < ${chromosome}.intervals > ${chromosome}.bed.tmp
        
        # Filter out any lines that don't have exactly 3 columns
        awk 'NF==3 {print}' ${chromosome}.bed.tmp > ${chromosome}.bed
        rm ${chromosome}.bed.tmp
        
        echo "BED file created with \$(wc -l < ${chromosome}.bed) valid intervals"
        
        # Show first few lines for verification
        if [ -s ${chromosome}.bed ]; then
            echo "First few intervals:"
            head -5 ${chromosome}.bed
        fi
    else
        echo "No intervals found, creating empty BED file"
        touch ${chromosome}.bed
    fi
    
    # Final verification
    if [ -s ${chromosome}.bed ]; then
        # Check that all lines have 3 columns
        BAD_LINES=\$(awk 'NF!=3 {print NR}' ${chromosome}.bed | wc -l)
        if [ \$BAD_LINES -gt 0 ]; then
            echo "WARNING: Found \$BAD_LINES malformed lines, filtering them out"
            awk 'NF==3 {print}' ${chromosome}.bed > ${chromosome}.bed.clean
            mv ${chromosome}.bed.clean ${chromosome}.bed
        fi
        echo "Final BED file has \$(wc -l < ${chromosome}.bed) valid intervals"
    fi
    """
}
