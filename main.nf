/* GATK Variant Calling Pipeline
 * Usage: nextflow run main.nf --ref --reads --sample --basedir
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
 * picard
 * qualimap
 */
 
//set default values for parameters
params.sequenced = false
//params.sequenced indicates whether or not it is raw read data from UNC, meaning it has the read group information intact, or was downloaded, meaning it does not
params.singleEnd = null
params.trim = 100
params.radseq = null
params.confident_ref = null

if( !params.reads ) error "Missing reads parameter"
if( !params.sample ) error "Missing sample parameter"
if( !params.sequenced && (!params.platform || params.platform instanceof Boolean )) error \
"If read group information not present in headers, you must specify the sequencing platform (e.g. Illumina)"

ref = file(params.ref)
sample = params.sample
trim = params.trim
intervals = params.intervals
outdir = "${params.basedir}/${params.sample}/"
adapter_dir = "/nas/longleaf/home/rrlove/.conda/envs/varaedes/opt/bbmap-38.90-0/resources/"

log.info """
VARAEDES
--------
genome  : $params.ref
reads   : $params.reads
sample  : $params.sample
current directory:  "$PWD"
basedir : "$outdir"
single-ended  :   $params.singleEnd
sequenced   :   $params.sequenced
confident_ref   :   $params.confident_ref
"""

reads = file(params.reads)

if ( reads.size() > 2 ) {

    Channel
    .fromPath( reads )
    .ifEmpty{ error "Cannot find any reads matching: ${params.reads}" }
    .view()
    .map { path -> subtags = (path.getBaseName() =~ /[\w\d\-\_]+_L0*(\d+)/)[0]; [sample, subtags[1], path] }
    .view()
    .groupTuple(by : [0, 1], sort : true)
    .view()
    .into{ reads_for_qc_ch; reads_for_trimming_ch }
}
else {

    Channel
    .fromPath ( reads )
    .ifEmpty{ error "Cannot find any reads matching: ${params.reads}" }
    .view()
    .map { path -> [sample, "1", path] }
    .view()
    .groupTuple(by : [0, 1], sort : true)
    .view()
    .into{ reads_for_qc_ch; reads_for_trimming_ch }
}

//reads_for_trimming_ch.view()

/*
 * First quality check
 * Tool: fastqc
 * Input: Channel reads_for_qc_ch
 * Output: Channel fastqc_pre_trim_ch
 */

process initialQC {
    publishDir "${outdir}qc/fastqc", mode: "copy"
    tag "$sample_$lane"
    //module 'fastqc'
    
    input: 
    tuple val(sample), val(lane), path(reads) from reads_for_qc_ch
    
    output: 
    tuple val(sample), val(lane), file("${sample}*fastqc.html"), file("${sample}*fastqc.zip")  into fastqc_pre_trim_ch

    script: 
    
    if ( params.singleEnd ){

    """
    fastqc ${reads}
    """
    }
    
    else {
    """
    fastqc ${reads[0]} ${reads[1]}
    """
    }
    
}

/*
Trim reads
Tool: trimmomatic
Input: Channel reads_for_trimming_ch
Output: Channel trimmed_reads_for_qc_ch ; trimmed_reads_for_align_ch ; 
trimmed_reads_for_rg_ch ; trimmed_reads_for_multiqc_ch
*/

process trim {
    publishDir "${outdir}data/trimmed_reads"
    tag "$sample_$lane"
    //module 'trimmomatic'

    input: 
    tuple val(sample), val(lane), path(reads) from reads_for_trimming_ch
    output: 
    tuple val(sample), val(lane), path("${sample}_${lane}_trimmed*P.fq.gz") \
    into trimmed_reads_for_qc_ch, trimmed_reads_for_align_ch, \
    trimmed_reads_for_rg_ch, trimmed_reads_for_multiqc_ch
    
    script:
    
    if ( params.singleEnd ) {
    """
    trimmomatic SE ${reads} ${sample}_${lane}_trimmed_P.fq.gz \
    ILLUMINACLIP:${adapter_dir}adapters.fa:2:30:10 LEADING:10 TRAILING:10 MINLEN:${trim}
    """
    } 
    
    else {
    """
    trimmomatic PE ${reads[0]} ${reads[1]} \
    -baseout ${sample}_${lane}_trimmed.fq.gz \
    ILLUMINACLIP:${adapter_dir}adapters.fa:2:30:10 LEADING:10 TRAILING:10 MINLEN:${trim}
    """
    }
}

/*
 * Post-trim quality check
 * Tool: fastqc
 * Input: Channel trimmed_reads_for_qc_ch
 * Output: Channel fastqc_post_trim_ch
*/

process secondQC {
    publishDir "${outdir}qc/fastqc", mode: "copy"
    tag "$sample_$lane"
    //module 'fastqc'

    input: 
    tuple val(sample), val(lane), path(trimmed_reads) from trimmed_reads_for_qc_ch

    output: 
    tuple val(sample), val(lane), file("${sample}*trimmed*fastqc.html"), file("${sample}*trimmed*fastqc.zip") into fastqc_post_trim_ch

    script: 
    
    if ( params.singleEnd ){

    """
    fastqc ${trimmed_reads}
    """
    }
    
    else {
    """
    fastqc ${trimmed_reads[0]} ${trimmed_reads[1]}
    """
    }
}

/*
 * Extract read group information to be passed to bwa
 * Tool: bash
 * Input: Channel trimmed_reads_for_rg_ch
 * Output: Channel rg_for_alignment_ch
*/

process extract_read_groups {
    tag "$sample_$lane"
    
    input:
    tuple val(sample), val(lane), path(merged_reads) from trimmed_reads_for_rg_ch
    
    output:
    tuple val(sample), val(lane), env(RG) into rg_for_alignment_ch
    
    shell:
    
    if ( params.sequenced ) {
    
    '''
    header=$(zcat !{sample}_!{lane}_trimmed*.fq.gz | head -n 1)
    IFS=':'; array=($header); unset IFS
    sample=!{sample}
    lb="lib-"${sample}
    pl="ILLUMINA"
    id=${array[2]}.${array[3]}.${array[4]}
    pu=${array[3]}.${array[2]}.${array[-1]}
    RG="@RG\\tID:"${id}"\\tLB:"${lb}"\\tPL:"${pl}"\\tSM:"${sample}"\\tPU:"${pu}
    
    '''
    }
    
    else {
    '''
    RG="@RG\\tID:"!{sample}"\\tLB:"!{sample}"\\tPL:"!{params.platform}"\\tSM:"!{sample}"\\tPU:"!{sample}
    '''
    }

//IFS=':' read -r -a ARRAY <<< (zcat "!{sample}_!{lane}_trimmed_*P.fq.gz" | head -n 1)


//@RG\tID:CL1000\tLB:Library_ID\tPL:HELICOS\tSM:SAMPLE_ID\tPU:Index_seq
//"@RG\tID:\tSM:\tPL:/tLB:\tPI:"

//ID: run + flowcell barcode + lane number
//SM: sample
//LB: "lib-"sample
//PL: ILLUMINA
//PU: flowcell barcode + lane number + sample barcode

//so, need to extract the run, the flowcell barcode, 
// the lane number, and the sample barcode

//in these data, these have the following fields (1-indexed):
//run: 2
//flowcell barcode: 3
//lane number: 4
//sample barcode: -1

//ID: 2 + 3 + 4
//SM: sample
//LB: "lib-"sample
//PL: ILLUMINA
//PU: 3 + 2 + -1
}

/*
 * Align reads
 * Tool: bwa-mem2
 * Input: Channel rg_for_alignment_ch, Channel trimmed_reads_for_align_ch
 * Output: Channel aligned_reads_ch, Channel bam_for_qc_ch
*/

process align_reads {
    tag "$sample_$lane"
    cpus 4
    memory "32G"
    //module 'samtools'
    
    input:
    tuple val(sample), val(lane), path(trimmed_reads) from trimmed_reads_for_align_ch
    tuple val(sample), val(lane), val(read_group) from rg_for_alignment_ch
    
    output:
    file("${sample}_${lane}_aligned.bam") into aligned_reads_ch
    
    """
    bwa-mem2 mem \
    -R "${read_group}" \
    -t 4 \
    ${ref} ${trimmed_reads} | samtools view -1bh > ${sample}_${lane}_aligned.bam
    """
}

aligned_reads_ch
    .collect()
    .set { aligned_reads_path_ch }
    //.view()
    
/*
 * Merge BAM files
 * Tool: samtools
 * Input: Channel aligned_reads
 * Output: Channel merged_bam
 * NB: this extra step is because of the difficulty in getting Picard to take multiple files from nextflow
*/
    
process merge_bams {
    tag "$sample"
    //module 'samtools'
    
    input:
    path(aligned_reads) from aligned_reads_path_ch
    
    output:
    tuple val(sample), file("${sample}_merged.bam") into merged_reads_ch
    
    """
    samtools merge ${sample}_merged.bam ${aligned_reads}
    """
}

/*
 * Mark duplicates
 * Tool: MarkDuplicatesSpark
 * Input: Channel merged_reads_ch
 * Output: Channel bam_dups_marked_ch
*/

/* MarkDuplicates shouldn't be used with radseq data: https://gatk.broadinstitute.org/hc/en-us/community/posts/4403892377371-Difference-in-RAD-seq-and-resequencing-data-analysis
*/

process mark_dups {
    publishDir "${outdir}data/bam"
    tag "sample"
    //module 'gatk'
    //module 'samtools'
    cpus 4
    memory "16G"
    
    input:
    tuple val(sample), path(merged_bam) from merged_reads_ch
    
    output: 
    tuple val(sample), file("${sample}_sorted_dedup.bam") into bam_dedup_merge_ch, bam_for_qualimap_ch, bam_for_picard_ch
    tuple val(sample), file("${sample}_sorted_dedup.bam.bai") into bam_index_ch

    script:
    
    if ( params.radseq ) {
    """
    picard SortSam \
    I=${merged_bam} \
    O=${sample}_sorted_dedup.bam \
    SORT_ORDER=coordinate
    
    samtools index ${sample}_sorted_dedup.bam
    """
    
    }

    else {
    """    
    gatk MarkDuplicatesSpark \
    -I ${merged_bam} \
    -M ${sample}_dedup_metrics.txt \
    -O ${sample}_sorted_dedup.bam    \
    --conf 'spark.executor.cores=4'
    
    samtools index ${sample}_sorted_dedup.bam
    """
    }

}

/*
 * Check the quality of the BAMs with bamqc
 * Tool: qualimap
 * Input: Channel bam_for_qualimap_ch
 * Outout: Channel qualimap_results_ch
*/

process bamqc_qc{
    publishDir "${outdir}qc/qualimap", mode: "copy"
    tag "$sample"
    memory { 4.GB * task.attempt }
    time '6d'
    errorStrategy 'retry'
    maxRetries 3
    //module 'qualimap'

    input:
    tuple val(sample), path(aligned_bam) from bam_for_qualimap_ch
    
    output:
    path("${sample}_sorted_dedup_stats*") into qualimap_results_ch

    """
    qualimap bamqc \
    -bam ${aligned_bam} \
    --java-mem-size=4G

    """
}

/*
 * Check the quality of the BAMs with picard
 * Tool: picard
 * Input: Channel bam_for_picard_ch
 * Outout: Channel picard_results_ch
*/

process picard_qc{
    publishDir "${outdir}qc/picard", mode: "copy"
    tag "$sample"
    memory "4G"
    time '6d'
    //module 'picard'
    
    input:
    tuple val(sample), path(aligned_bam) from bam_for_picard_ch
    
    output:
    path("${sample}_picard_metrics*") into picard_results_ch

    """
    picard CollectMultipleMetrics \
    -I ${aligned_bam} \
    -O ${sample}_picard_metrics
    """
}

/*
 * Genotype variants in the specimen, splitting across chromosomes/scaffolds
 * Tool: GATK
 * Input: Channel bam_dedup_merge_ch, intervals_ch
 * Output: Channel genotyped_variant_chunks_ch, realigned_bam_chunks_ch
*/

Channel
    .fromPath( intervals )
    //.view()
    .set{ intervals_ch }

process call_variants{
    //publishDir "${outdir}data/vcf/chunked", pattern: '*vcf*'
    publishDir "${outdir}data/vcf"
    //publishDir "${outdir}data/bam", pattern: '*bam*'
    tag "$sample"
    //module 'gatk'
    memory "8G"
    
    input:
    tuple val(sample), path(bam) from bam_dedup_merge_ch
    path( intervals ) from intervals_ch

    output:
    //this should be: tuple val(sample), file("${sample}.${intervals.simpleName}.g.vcf.gz") into genotyped_variant_chunks_ch
    tuple val(sample), file("${sample}.${intervals.simpleName}.g.vcf.gz") into genotyped_variant_chunks_ch

    script:
    
    if ( params.confident_ref ) {
    """
    gatk HaplotypeCaller \
    -R $ref \
    -I ${bam} \
    -O ${sample}.${intervals.simpleName}.g.vcf.gz \
    -ERC GVCF \
    --output-mode EMIT_ALL_CONFIDENT_SITES \
    -L ${intervals}
    """
    }
    
    else {
    """
    gatk HaplotypeCaller \
    -R $ref \
    -I ${bam} \
    -O ${sample}.${intervals.simpleName}.g.vcf.gz \
    -ERC GVCF \
    -L ${intervals}
    """
    }
}

/*
 * Merge gVCFs
 * Tool: bash, picard
 * Input: Channel genotyped_variant_chunks_ch
 * Output: Channel merged_gvcfs_ch

//modified on 02-01-2022 to take merge_gvcfs out of the pipeline, as genotyping now runs on the per-interval gVCFs

process merge_gvcfs{
    publishDir "${outdir}data/vcf"
    tag "$sample"
    memory "8G"
    
    input:
    file(vcfsList) from genotyped_variant_chunks_ch.toSortedList( { path -> path.getBaseName() } )
    
    output:
    tuple val(sample), path("${sample}.g.vcf.gz") into merged_gvcfs_ch

    shell:

    '''
    for file in !{vcfsList}
    do
        echo ${file} >> !{sample}_vcfs_to_merge.list
    
    done

    picard GatherVcfs \
    --INPUT !{sample}_vcfs_to_merge.list \
    --OUTPUT !{sample}.g.vcf.gz \
    '''
}
*/

/*
 * Index gVCF chunks
 * Tool: GATK
 * Input: Channel genotyped_variant_chunks_ch
 * Output: Channel indexed_gvcfs_ch
*/

process index_gvcfs{
    publishDir "${outdir}data/vcf"
    tag "$sample"
    memory "8G"
    
    input:
    tuple val(sample), path("${sample}.*.g.vcf.gz") from genotyped_variant_chunks_ch
    
    output:
    tuple val(sample), path("${sample}.*.g.vcf.gz.tbi") into indexed_gvcfs_ch
    
    script:
    """
    gatk IndexFeatureFile \
    -I ${sample}.*.g.vcf.gz
    """
}

/*
 * Compile the quality reports
 * Tool: multiqc
 * Input: Channel fastqc_pre_trim_ch, Channel trimmed_reads_for_multiqc_ch, Channel fastqc_post_trim_ch, Channel qualimap_results_ch, Channel picard_results_ch
 * Output: Channel multiqc_output_ch
*/

process multiqc{
    publishDir "${outdir}qc/multiqc", mode: "copy"
    tag "${sample}"
    
    input:
    tuple val(sample), val(lane), file("${sample}_L00${lane}_*fastqc.html"), file("${sample}_L00${lane}_*fastqc.zip") from fastqc_pre_trim_ch
    tuple val(sample), val(lane), path("${sample}_${lane}_trimmed_*P.fq.gz") from trimmed_reads_for_multiqc_ch
    tuple val(sample), val(lane), file("${sample}_${lane}_trimmed_*P*fastqc.html"), file("${sample}_${lane}_trimmed_*P*fastqc.zip") from fastqc_post_trim_ch
    path("${sample}_sorted_dedup_stats*/qualimapReport.html") from qualimap_results_ch
    path("${sample}_picard_metrics*") from picard_results_ch
    
    output:
    path("${sample}_multiqc_report.html") into multiqc_output_ch
    
    script:
    
    """
    multiqc \
    -n "${sample}"_multiqc_report.html \
    -c /overflow/dschridelab/users/rrlove/aedes/pipeline/multiqc_config.yaml \
    "${outdir}"qc/fastqc \
    "${outdir}"qc/picard \
    "${outdir}"qc/qualimap \
    "${outdir}"data/trimmed_reads \
    "${outdir}"data/bam
    """
    
}

//make output not modifiable
//validate choice of trimmomatic parameters
//find good heterozygosity parameters for HaplotypeCaller