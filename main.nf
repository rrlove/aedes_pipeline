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

if( !params.reads ) error "Missing reads parameter"
if( !params.sample ) error "Missing sample parameter"
println "sample: $params.sample"

ref = file(params.ref)
println ref
outdir = "${baseDir}/${params.sample}/"
scratchdir = "/pine/scr/r/r/rrlove/aedes_pipeline/"
println outdir
adapter_dir = "/nas/longleaf/apps/trimmomatic/0.36/Trimmomatic-0.36/adapters/"
sample = params.sample

log.info """
VARAEDES
--------
genome  : $params.ref
reads   : $params.reads
sample  : $params.sample
current directory:  "$PWD"
basedir : "$outdir"

"""


Channel
    .fromPath( params.reads )
    .ifEmpty{ error "Cannot find any reads matching: ${params.reads}" }
    //.view()
    .map { path -> subtags = (path.getBaseName() =~ /([\w\d\-]+_S\d+)_L0+(\d+)/)[0]; [sample, subtags[2], path] }
    //.view()
    .groupTuple(by : [0, 1])
    //.view()
    .into{ reads_for_qc_ch; reads_for_trimming_ch }
    
/*
 * First quality check
 * Tool: fastqc
 * Input: Channel reads_for_qc_ch
 * Output: Channel fastqc_pre_trim_ch
 */


process initialQC {
    publishDir "${outdir}qc/fastqc"
    tag "$sample_$lane"
    //module 'fastqc'
    
    input: 
    tuple val(sample), val(lane), path(reads) from reads_for_qc_ch
    
    output: 
    tuple val(sample), val(lane), file("${sample}_L00${lane}_*fastqc.html") into fastqc_pre_trim_ch

    script: 

    """
    fastqc ${reads[0]} ${reads[1]}
    """

}

/*
Trim reads
Tool: trimmomatic
Input: Channel reads_for_trimming_ch
Output: Channel trimmed_reads_for_qc_ch ; trimmed_reads_for_align_ch 
*/

process trim {
    publishDir "${outdir}data/trimmed_reads"
    tag "$sample_$lane"
    //module 'trimmomatic'

    input: 
    tuple val(sample), val(lane), path(reads) from reads_for_trimming_ch
    output: 
    tuple val(sample), val(lane), path("${sample}_${lane}_trimmed_*P.fq.gz") into trimmed_reads_for_qc_ch, trimmed_reads_for_align_ch, trimmed_reads_for_rg_ch
    
    script:
    
    """
    trimmomatic PE ${reads[0]} ${reads[1]} -baseout ${sample}_${lane}_trimmed.fq.gz \
    ILLUMINACLIP:${adapter_dir}TruSeq3-PE.fa:2:30:10 LEADING:10 TRAILING:10 MINLEN:70
    """
}

/*
 * Post-trim quality check
 * Tool: fastqc
 * Input: Channel trimmed_reads_for_qc_ch
 * Output: Channel fastqc_post_trim_ch
*/

process secondQC {
    publishDir "${outdir}qc/fastqc"
    tag "$sample_$lane"
    //module 'fastqc'

    input: 
    tuple val(sample), val(lane), path(trimmed_reads) from trimmed_reads_for_qc_ch

    output: 
    tuple val(sample), val(lane), file("${sample}_${lane}_trimmed_*fastqc.html") into fastqc_post_trim_ch

    script: 
    
    """
    fastqc ${trimmed_reads[0]} ${trimmed_reads[1]}
    """
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
    
    '''
    header=$(zcat !{sample}_!{lane}_trimmed_*P.fq.gz | head -n 1)
    IFS=':'; array=($header); unset IFS
    sample=!{sample}
    lb="lib-"${sample}
    pl="ILLUMINA"
    id=${array[2]}.${array[3]}.${array[4]}
    pu=${array[3]}.${array[2]}.${array[-1]}
    RG="@RG\\tID:"${id}"\\tLB:"${lb}"\\tPL:"${pl}"\\tSM:"${sample}"\\tPU:"${pu}
    '''


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
    publishDir "${outdir}data/bam"
    tag "$sample_$lane"
    cpus 4
    memory "16G"
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
    //.flatten()
    //.filter( Path )
    .view()
    .set { aligned_reads_path_ch }
    
/*
 * Merge BAM files
 * Tool: samtools
 * Input: Channel aligned_reads
 * Output: Channel merged_bam
 * NB: this extra step is because of the difficulty in getting Picard to take multiple files from nextflow
*/

process merge_reads {
    publishDir "${outdir}data/bam"
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
    tuple val(sample), file("${sample}_sorted_dedup.bam") into bam_dedup_merge_ch, bam_for_qc_ch
    tuple val(sample), file("${sample}_sorted_dedup.bam.bai") into bam_index_ch

    """
    gatk MarkDuplicatesSpark \
    -I ${merged_bam} \
    -M ${sample}_dedup_metrics.txt \
    -O ${sample}_sorted_dedup.bam    \
    --conf 'spark.executor.cores=4'
    
    samtools index ${sample}_sorted_dedup.bam
    """

}

/*
 * Check the quality of the BAMs
 * Tool: qualimap, picard
 * Input: Channel bam_for_qc_ch
 * Outout: Channel initial_bam_qc_ch
*/

process bam_qc{
    publishDir "${outdir}qc/qualimap", pattern: '*stats*'
    publishDir "${outdir}qc/picard", pattern: '*picard*'
    tag "$sample"
    memory "4G"
    //module 'qualimap'
    //module 'picard'
    
    input:
    tuple val(sample), path(aligned_bam) from bam_for_qc_ch
    
    output:
    tuple path("${sample}_sorted_dedup_stats*"), path("${sample}_picard_metrics*") into initial_bam_qc_ch

    """
    qualimap bamqc \
    -bam ${aligned_bam} \
    --java-mem-size=4G
    
    picard CollectMultipleMetrics \
    -I ${aligned_bam} \
    -O ${sample}_picard_metrics
    """
}

/*
 * Genotype variants in the specimen
 * Tool: GATK
 * Input: Channel bam_dedup_merge_ch
 * Output: Channel genotyped_variants_ch, realigned_bam_ch
*/

process call_variants{
    publishDir "${outdir}data/vcf", pattern: '*vcf*'
    publishDir "${outdir}data/bam", pattern: '*bam*'
    tag "$sample"
    //module 'gatk'
    memory "8G"
    
    input:
    tuple val(sample), path(bam) from bam_dedup_merge_ch

    output:
    tuple val(sample), path("${sample}.g.vcf.gz") into genotyped_variants_ch
    tuple val(sample), path("${sample}_realigned.bam") into realigned_bam_ch

    script:
    """
    gatk HaplotypeCaller \
    -R $ref \
    -I ${bam} \
    -O ${sample}.g.vcf.gz \
    -ERC GVCF \
    -bamout ${sample}_realigned.bam
    """
}

//set up conda environment and remove "module" calls
//make output not modifiable
//validate choice of trimmomatic parameters

