
/*
 * Mark duplicates
 * Tool: MarkDuplicatesSpark...? / Picard MarkDuplicates + SortSam
 * Input: Channel aligned_bam_ch
 * Output: Channel bam_dups_marked_ch


process mark_dups {
    publishDir "${outdir}data/bam"
    tag "sample"
    module 'gatk'
    cpus 4
    memory "8G"
    
    input: 
    tuple val(sample), val(lane), path("*.bam") from aligned_reads
    
    output: 
    tuple val(sample), file("${sample}_sorted_dedup.bam") into bam_dedup_merge_ch

    """

    gatk MarkDuplicatesSpark \
    -I "${sample}_*_aligned.bam" \
    -M ${sample}_dedup_metrics.txt \
    -O ${sample}_sorted_dedup.bam    \
    --conf 'spark.executor.cores=4'
    
    """

}
*/

/*
 * Check the quality of the BAMs
 * Tool: qualimap, picard
 * Input: Channel bam_for_qc
 * Outout: Channel initial_bam_qc


process bam_qc{
    publishDir "${outdir}qc/qualimap"
    tag "$sample_$lane"
    module 'qualimap'
    module 'picard'
    
    input:
    tuple val(sample), val(lane), path(aligned_bam) from bam_for_qc
    output:
    tuple val(sample), val(lane), file("${sample}_${lane}_bamqc.html") into initial_bam_qc
    
    """
    qualimap bamqc \
    -bam ${aligned_bam} \
    -outdir ${outdir}qc/qualimap
    
    picard CollectMultipleMetrics \
    -I ${aligned_bam} \
    -O ${sample}_picard_metrics
    """
}
*/


/*
 * Merge files from both lanes
 * Tool: bash
 * Input: Channel trimmed_reads_for_align_ch
 * Output: Channel merged_reads_ch
*/

/*
process merge_reads {
//where exactly do I want this? probably safest to merge bams...

    publishDir "${outdir}data/merged"
    tag "$sample_$lane"
    
    input:
    tuple val(sample), val(lane), path(trimmed_reads) from trimmed_reads_for_align_ch
    
    output:
    tuple val(sample), file("${sample}_trimmed_merged_{1,2}P.fq.gz") into merged_reads_ch, get_read_groups_ch
    
    script:
    
    """
    zcat ${sample}_*_trimmed_1P.fq.gz > ${sample}_trimmed_merged_1P.fq.gz
    zcat ${sample}_*_trimmed_2P.fq.gz > ${sample}_trimmed_merged_2P.fq.gz
    """

}
*/

