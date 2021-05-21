/* Contamination pipeline
 * Usage: nextflow run main.nf --bam --sample --basedir --threads --subsample_frac --n_partitions
 *
 * Author: RRLove < rrlove@email.unc.edu >
 * University of North Carolina Chapel Hill, 2021
*/

bam = params.bam
sample = params.sample
basedir = params.basedir
threads = params.threads
subsample_frac = params.subsample_frac
n_partitions = params.n_partitions

log.info """
--------
bam  : $bam
sample  : $sample
basedir : $basedir
threads : $threads
subsample_frac : $subsample_frac
current directory:  $PWD
"""

Channel
    .fromPath( bam )
    .ifEmpty{ error "Cannot find any file matching: ${bam}" }
    .set { bam_ch }
    
/*
 * Extract unmapped reads
 * Tool: samtools
 * Input: Channel bam_ch
 * Output: Channel unmapped_fasta_ch
 */

process extract_unmapped_reads {
    publishDir "${basedir}", mode: "copy"
    tag "$sample"
    module 'samtools'
    
    input: 
    path(bam) from bam_ch
    
    output: 
    tuple val(sample), path("${sample}_unmapped.fasta") into unmapped_fasta_ch

    script: 

    """
    samtools view -b -f 4 $bam | samtools sort -n | samtools fasta > ${sample}_unmapped.fasta
    """

}

/*
 * Subsample reads
 * Tool: seqtk
 * Input: Channel unmapped_fasta_ch
 * Output: Channel subsampled_ch
*/

process subsample_reads {
    publishDir "${basedir}", mode: "copy"
    tag "$sample"
    module 'seqtk'

    input: 
    tuple val(sample), path("${sample}_unmapped.fasta") from unmapped_fasta_ch
    
    output: 
    tuple val(sample), path("${sample}_subsampled.fasta") into subsampled_ch
    
    script:
    
    """
    seqtk sample ${sample}_unmapped.fasta ${subsample_frac} > ${sample}_subsampled.fasta
    
    """
}
    
Channel
    .of(0.."$n_partitions".toInteger())
    //.view()
    .map{ it: it.toString().padLeft( 2, '0' ) }
    .collect()
    .set{ blast_partitions_suffix_ch }

/*
 * BLAST unmapped reads
 * Tool: BLAST
 * Input: Channel subsampled_ch, blast_partitions_suffix_ch
 * Output: Channel blast_ch
*/

process blast_unmapped_reads {
    publishDir "${basedir}", mode: "copy"
    tag "$sample"
    module 'blast'
    //memory '8G'
    //clusterOptions "-n $threads"

    input: 
    tuple val(sample), path("${sample}_subsampled.fasta") from subsampled_ch
    each suffix from blast_partitions_suffix_ch

    output: 
    tuple val(sample), path("${sample}_BLAST_hits.${suffix}.txt") into blast_ch

    script: 
    
    """

    blastn \
    -db nt.${suffix} \
    -query ${sample}_subsampled.fasta \
    -outfmt "6 std staxid ssciname scomname"\
    -num_threads ${threads} \
    -max_target_seqs 1 \
    -perc_identity 90 \
    -out ${sample}_BLAST_hits.${suffix}.txt
    """
}

