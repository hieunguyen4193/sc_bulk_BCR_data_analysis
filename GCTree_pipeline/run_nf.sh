path_to_samplesheet=$1;
nextflow run GCtree_pipeline_input_SampleSheet.nf --samplesheet ${path_to_samplesheet} -resume  