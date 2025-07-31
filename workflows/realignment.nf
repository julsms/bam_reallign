include { SPLIT_BAM } from '../modules/split_bam'
include { BAMLEFTALIGN } from '../modules/bamleftalign'
include { ADD_READ_GROUPS } from '../modules/add_read_groups'
include { REALIGNER_TARGET_CREATOR } from '../modules/realigner_target_creator'
include { EXPAND_REGIONS } from '../modules/expand_regions'
include { ABRA_REALIGNMENT } from '../modules/abra_realignment'
include { MERGE_BAM } from '../modules/merge_bam'

workflow REALIGNMENT {
    
    // Validate required parameters
    if (!params.input_bam) {
        error "Please provide --input_bam parameter"
    }
    if (!params.reference_fasta) {
        error "Please provide --reference_fasta parameter"
    }
    if (!params.genome_file) {
        error "Please provide --genome_file parameter (chr\tlength format)"
    }
    
    // Create input channels
    input_bam_ch = Channel.fromPath(params.input_bam, checkIfExists: true)
    input_bai_ch = Channel.fromPath("${params.input_bam}.bai", checkIfExists: true)
    reference_fasta_ch = Channel.fromPath(params.reference_fasta, checkIfExists: true)
    reference_fai_ch = Channel.fromPath("${params.reference_fasta}.fai", checkIfExists: true)
    reference_dict_ch = Channel.fromPath("${params.reference_fasta.replaceAll(/\.fa(sta)?$/, '.dict')}", checkIfExists: true)
    genome_file_ch = Channel.fromPath(params.genome_file, checkIfExists: true)
    
    // Step 1: Split BAM by chromosome
    SPLIT_BAM(input_bam_ch, input_bai_ch)
    
    // Step 2: Create chromosome-specific channels for parallel processing
    chr_bams_ch = SPLIT_BAM.out.chr_bams
        .flatten()
        .filter { it.name.endsWith('.bam') }
        .map { bam_file -> 
            def chr_name = bam_file.baseName
            [chr_name, bam_file]
        }
    
    // Step 3: Add read groups FIRST (PARALLEL by chromosome)
    ADD_READ_GROUPS(
        chr_bams_ch,
        reference_fasta_ch.first(),
        reference_fai_ch.first(),
        reference_dict_ch.first()
    )
    
    // Step 4: Left-align indels AFTER adding read groups (PARALLEL by chromosome)
    BAMLEFTALIGN(
        reference_fasta_ch.first(),
        reference_fai_ch.first(),
        reference_dict_ch.first(),
        ADD_READ_GROUPS.out.rg_bams.map { chr, bam, bai -> [chr, bam] }
    )
    
    // Step 5: Identify realignment targets (PARALLEL by chromosome)
    REALIGNER_TARGET_CREATOR(
        reference_fasta_ch.first(),
        reference_fai_ch.first(),
        reference_dict_ch.first(),
        BAMLEFTALIGN.out.leftaligned_bams
    )
    
    // Step 6: Expand regions (PARALLEL by chromosome)
    EXPAND_REGIONS(
        genome_file_ch.first(),
        REALIGNER_TARGET_CREATOR.out.target_intervals
    )
    
    // Step 7: Combine left-aligned BAMs and expanded regions for ABRA
    abra_input_ch = BAMLEFTALIGN.out.leftaligned_bams
        .join(EXPAND_REGIONS.out.expanded_regions, by: 0)
        .map { chr_name, leftaligned_bam, expanded_bed -> 
            [chr_name, leftaligned_bam, expanded_bed]
        }
    
    // Step 8: ABRA realignment (PARALLEL by chromosome)
    ABRA_REALIGNMENT(
        reference_fasta_ch.first(),
        abra_input_ch
    )
    
    // Step 9: Collect and merge
    realigned_bams_ch = ABRA_REALIGNMENT.out.realigned_bams
        .map { chr_name, realigned_bam -> realigned_bam }
        .collect()
    
    MERGE_BAM(realigned_bams_ch)
    
    emit:
    final_bam = MERGE_BAM.out.final_bam
    final_bai = MERGE_BAM.out.final_bai
}
