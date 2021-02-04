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
    tuple val(sample), val(lane), path("${sample}_${lane}_trimmed_*P.fq.gz") into trimmed_reads_for_qc_ch, trimmed_reads_for_align_ch
    
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
 * Merge files from both lanes
 * Tool: bash
 * Input: Channel trimmed_reads_for_align_ch
 * Output: Channel merged_reads_ch
*/

process mergereads {
    publishDir "${outdir}/data/merged"
    tag "$sample_$lane"
    
    input:
    tuple val(sample), val(lane), path(trimmed_reads) from trimmed_reads_for_align_ch
    
    output:
    tuple val(sample), file("${sample}_trimmed_merged_{1,2}P.fq.gz") into merged_reads_ch
    
    script:
    
    """
    zcat ${sample}_*_trimmed_1P.fq.gz > ${sample}_trimmed_merged_1P.fq.gz
    zcat ${sample}_*_trimmed_2P.fq.gz > ${sample}_trimmed_merged_2P.fq.gz
    """


}


//NEED TO EXTRACT READ GROUPS
