/* GATK Variant Calling Pipeline
 * Usage: nextflow run main.nf --reference --reads --etc
 *
 * Author: RRLove < rrlove@email.unc.edu >
 * University of North Carolina Chapel Hill, 2021

*/


/*
 * running list of packages for the conda environment:
 * fastqc
 * trimmomatic
 * bwa-mem2
 * gatk
 
 */

params.reads = null 
if( !params.sample ) error "Missing sample parameter"
//println "sample: $params.sample"
params.outdir = "${baseDir}/${prefix}"
params.ref = "${baseDir}/path/to/fasta.fa"

ref = path(params.ref)
outdir = val(params.outdir)

log.info """
VARAEDES
--------
genome  : $params.ref
reads   : $params.reads
sample  : $params.sample

"""

Channel
    .fromFilePairs( param.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .into{ reads_for_qc_ch; reads_for_trimming_ch }

/*
 * First quality check
 * Tool: fastqc
 * Input: Channel reads_for_qc_ch
 * Output: Channel fastqc_pre_trim_ch
 */

process initialQC {
    tag "$pair_id"
    module 'fastQC'
    publishDir "$baseDir/data/fastqc"
    
    input: 
    tuple val(pair_id), path(reads) from reads_for_qc_ch
    
    output: 
    tuple val(pair_id), path("${pair_id}.fastqc.html") into fastqc_pre_trim_ch

    script: 

    """
    fastqc ${reads[0]} ${reads[1]} -o .
    """

}

/*
Trim reads
Tool: trimmomatic (?)
Input: Channel reads_for_trimming_ch
Output: Channel trimmed_reads_for_qc_ch ; trimmed_reads_for_align_ch 
*/

process trim {
    tag "$pair_id"
    module 'trimmomatic'
    publishDir

    input: 
    tuple val(pair_id), path(reads) from reads_for_trimming_ch
    output: 
    tuple val(pair_id), path(trimmed_reads) into trimmed_reads_for_qc_ch, trimmed_reads_for_align_ch
    
    script:
    
    """
    trimmomatic -basein ${reads[0]} ${reads[1]} -baseout ${pair_id}_trimmed
    """
}

/*
 * Post-trim quality check
 * Tool: fastqc
 * Input: Channel trimmed_reads_ch
 * Output: Channel fastqc_post_trim_ch
*/

process secondQC {
    tag "$pair_id"
    module 'fastQC'
    publishDir "$baseDir/data/fastqc"
    
    input: 
    tuple(pair_id), path(trimmed_reads) from trimmed_reads_for_qc_ch
    
    output: 
    tuple val(pair_id), path("${pair_id}_trimmed.fastqc.html") into fastqc_post_trim_ch

    script: 
    
    """
    fastqc ${trimmed_reads[0]} ${trimmed_reads[1]} -o .
    """

}

/*
 * Merge files from both lanes
 * Tool: bash
 * Input: Channel trimmed_reads_for_align_ch
 * Output: Channel merged_reads_ch
*/

process mergereads {


}


/*
 * Align reads
 * Tool: bwa-mem
 * Input: Channel merged_reads_ch
 * Output: Channel aligned_bam_ch
*/

// need to check for the genome index file, and if it doesn't exist, index it

process align {
    tag "$pair_id"
    module 'bwa'
    publishDir
    cpus 4
    
    input: 
    tuple(pair_id), path(trimmed_reads) from trimmed_reads_ch
    
    output: 
    tuple(pair_id), path("${pair_id}.aligned.bam") into aligned_bam_ch
    
    script:
    
    """
    bwa-mem2 \
    -R \// includes complete read groups
    -t 4 \
    ${ref} ${trimmed_reads} \ | samtools view -1bh > ${pair_id}.aligned.bam
    """

}

/*
 * Mark duplicates
 * Tool: MarkDuplicatesSpark...? / Picard MarkDuplicates + SortSam
 * Input: Channel aligned_bam_ch
 * Output: Channel bam_dups_marked_ch
*/

process mark_dups{
    tag "$pair_id"
    module 'gatk'
    cpus 4
    memory "8G"
    
    input: 
    tuple(pair_id), path(aligned_bam) from aligned_bam_ch
    
    output: 
    tuple(pair_id), path(deduped_bam) into bam_dups_marked_ch
    
    script:
    
    """
    gatk MarkDuplicatesSpark \
    -I ${aligned_bam} \
    -M ${pair_id}_dedup_metrics.txt \
    -O ${pair_id}_sorted_dedup.bam \
    --conf 'spark.executor.cores=4'
    
    """

}

/*
 * Call variants
 * Tool: GATK HaplotypeCaller
 * Input: Channel bam_dups_marked_ch
 * Output:
*/

process call_variants{
    
    publishDir
    
    input:
    tuple(pair_id), path(sorted_deduped_bam) from bam_dups_marked_ch
    
    output:
    tuple(pair_id), path(genotyped_variants) into genotyped_variants_ch

    script:
    """
    gatk HaplotypeCaller -R $ref -I input.bam -O output.vcf
    """
}

/*
 * Genotype variants jointly
 * Tool: GATK HaplotypeCaller
 * Input: files..., reference
 * Output: files...
*/













