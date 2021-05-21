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
    .view()
    .map { path -> subtags = (path.getBaseName() =~ /([\w\d\-]+_S\d+)_L0+(\d+)/)[0]; [sample, subtags[2], path] }
    .view()
    .groupTuple(by : [0, 1])
    .view()
    .into{ reads_for_qc_ch; reads_for_trimming_ch }
    
