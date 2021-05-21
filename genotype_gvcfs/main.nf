/* GATK Variant Genotyping Pipeline
 * Usage: nextflow run main.nf --ref --basedir --intervals --update --mapfile --workspace_path --lazy
 *
 * Author: RRLove < rrlove@email.unc.edu >
 * University of North Carolina Chapel Hill, 2021
*/

//general use:

//collate GVCFs from all specimens into a GenomicsDB
//genotype GVCF

params.update = false
params.lazy = false

ref = file(params.ref)
basedir = params.basedir
intervals = params.intervals
mapfile = params.mapfile
workspace_path = params.workspace_path

if( !params.ref || params.ref instanceof Boolean ) error "Missing reference genome"
if( !params.basedir || params.basedir instanceof Boolean ) error "Missing base directory path"
if( !params.intervals || params.intervals instanceof Boolean ) error "Missing path to genomic intervals"
if( !params.mapfile || params.mapfile instanceof Boolean ) error "Missing path to sample map file"
if( params.update && (!params.workspace_path || params.workspace_path instanceof Boolean )) error "Missing path to desired output location"
if( params.lazy && !params.update ) error "Lazy sample addition can only be used on updates"

outdir = "${basedir}/vcf"
scratchdir = "/pine/scr/r/r/rrlove/aedes_pipeline/"
date = new java.util.Date().format( 'MMddyy' )


log.info """
GENOTYPE_GVCFs
--------
genome  : $params.ref
current directory:  "$PWD"
basedir : $basedir
outdir  : $outdir
update? : $params.update
lazy?   :   $params.lazy

"""

//Channel
//    .fromPath( intervals )
    //.view()
//    .into{ intervals_builddb_ch; intervals_updatedb_ch }
    
if ( mapfile ) {
    Channel
        .fromPath( mapfile )
        .set{ mapfile_ch }
} 
    

/*
 * When run in lazy mode, collect all g.vcf.gz files that are newer than the existing map file
 * Tool: shell scripting
 * Input: none
 * Output: Channel map_ch
*/

process get_changes {
    output:
    path("lazy_samples.map") into map_ch
    
    when:
    params.lazy

    shell:
    
    '''
    since=$(stat -c %y !{mapfile})
    
    find . -name "*.g.vcf.gz" -newermt ${since} -not -path "*/work/*" -exec sh -c "for f do name=$(basename '$f' .g.vcf.gz); echo -e '$name''\t''$f' ; done" find-sh {} ; > lazy_samples.map

    '''
}


/*
 * Build the database if it does not yet exist
 * Tool: GATK GenomicsDBImport
 * Input: Channel intervals_ch
 * Output: Channel database_ch
*/

process build_database {
    publishDir "${outdir}"
    memory "16G"
    
    input:
    path(mapfile) from map_ch.mix(mapfile_ch)
    
    output:
    path("${workspace_path}") into database_built_ch
    
    when:
    !params.update

    script:
    
    """
    
    gatk --java-options "-Xmx13G" GenomicsDBImport \
    --genomicsdb-workspace-path ${workspace_path} \
    --sample-name-map ${mapfile} \
    --intervals ${intervals}
    """
}

/*
 * Otherwise, update the database if it already exists
 * Tool: GATK GenomicsDBImport
 * Input: Channel intervals_ch
 * Output: Channel database_ch
*/

process update_database {
    publishDir "${outdir}"
    memory "16G"
    
    input:
    path( mapfile ) from map_ch
    
    output:
    path("${workspace_path}") into database_updated_ch
    
    when:
    params.update

    script:
    
    """
    gatk --java-options "-Xmx13G" GenomicsDBImport \
    --genomicsdb-update-workspace-path "${workspace_path}" \
    --sample-name-map ${mapfile} \
    --intervals ${intervals}
    
    """
}


/*
 * Genotype the variants using the database
 * Tool: GATK GenotypeGVCFs
 * Input: Channel database_ch
 * Output: genotyped_vcf_ch
*/

process genotype_GVCFs {
    publishDir "${outdir}"
    memory "16G"
    
    input:
    path(database) from database_built_ch.mix(database_updated_ch)
    
    output:
    path("merged*.vcf.gz") into merged_vcf_ch
    
    script:
    
    """
    gatk --java-options "-Xmx13G" GenotypeGVCFs \
    -R ${ref} \
    -V gendb://${database} \
    -O "merged_${date}.vcf.gz"
    
    """
}


