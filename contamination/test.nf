/*
stop = 25

println "$stop"

Channel
    .of(1.."$stop".toInteger())
    .view()
*/
    
//stop = $(ls /nas/longleaf/data/blast/nt*md5 | wc -l)

//println "$stop"

//val { id_in.toString().padLeft( 2, '0' ) } into transformed

n_partitions = 25

Channel
    .of(0.."$n_partitions".toInteger())
    //.view()
    .map{ it: it.toString().padLeft( 2, '0' ) }
    .collect()
    .view()
    .set{ blast_partitions_suffix_ch }