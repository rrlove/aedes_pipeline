This repository hosts the nextflow scripts used to align reads and call and genotype variants for Aedes aegypti sequenced as part of Love et al. 2022.

The first script, main.nf, is used to align reads, call variants, and perform various QC steps. An example slurm script showing its use can be found in examples/align_and_call_variants.job.

The second script, genotype_gvcfs/main.nf, is used to batch genotype variants from the output of the first script. It is set up to parallelize genotyping in several chunks per chromosome. An example slurm script showing its use can be found in examples/genotype_variants.job.

The nextflow.config scripts are specific to our HPC setup, but may be useful as examples. varaedes.txt shows the contents of the associated conda package used with these nextflow workflows.
