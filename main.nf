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
 
 */

if( !params.reads ) error "Missing reads parameter"
if( !params.sample ) error "Missing sample parameter"
println "sample: $params.sample"

ref = file(params.ref)
println ref
outdir = "${baseDir}/${params.sample}/"
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
    publishDir "${outdir}/qc/fastqc"
    tag "$sample_$lane"
    module 'fastqc'
    
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
    publishDir "${outdir}/data/trimmed_reads"
    tag "$sample_$lane"
    module 'trimmomatic'

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
    publishDir "${outdir}/qc/fastqc"
    tag "$sample_$lane"
    module 'fastqc'

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
 * Output: Channel aligned_reads
*/

process align_reads{
    publishDir "${outdir}/data/bam"
    tag "$sample_$lane"
    cpus 4
    module 'samtools'
    
    input:
    tuple val(sample), val(lane), path(trimmed_reads) from trimmed_reads_for_align_ch
    tuple val(sample), val(lane), val(read_group) from rg_for_alignment_ch
    
    output:
    tuple val(sample). val(lane), file("${sample}_${lane}_aligned.bam") into aligned_reads
    
    """
    bwa-mem2 mem \
    -R "${read_group}" \
    -t 4 \
    ${ref} ${trimmed_reads} | samtools view -1bh > ${sample}_${lane}_aligned.bam
    """

}



/*
 * Merge files from both lanes
 * Tool: bash
 * Input: Channel trimmed_reads_for_align_ch
 * Output: Channel merged_reads_ch
*/

/*
process merge_reads {
//where exactly do I want this? probably safest to merge bams...

    publishDir "${outdir}/data/merged"
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




//NEED TO EXTRACT READ GROUPS
//make output not modifiable
//validate choice of trimmomatic parameters

