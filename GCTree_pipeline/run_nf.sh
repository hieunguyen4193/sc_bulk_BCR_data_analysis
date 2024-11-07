path_to_samplesheet=$1;
output=$2;
deduplicate_src=$3;

nextflow run GCtree_pipeline_input_SampleSheet.nf \
--samplesheet ${path_to_samplesheet} \
--deduplicate_src ${deduplicate_src} \
--output ${output} -resume  