path_to_samplesheet=$1;
output=$2;
deduplicate_src=$3;
path_to_work=$4;

nextflow run GCtree_pipeline_input_SampleSheet.nf \
--samplesheet ${path_to_samplesheet} \
--deduplicate_src ${deduplicate_src} \
--output ${output} -resume  -w ${path_to_work}

# for mouseID in m1 m2 m3;do for inputCase in all without_colon_sample;do \
# bash run_nf.sh \
# /home/hieu/src/sc_bulk_BCR_data_analysis/SampleSheet_GCTree_nextflow_240826_BSimons/SampleSheet_GCTree_240826_BSimons_${mouseID}_${inputCase}.nextflow.csv \
# /home/hieu/outdir/ \
# /home/hieu/src/sc_bulk_BCR_data_analysis/GCTree_pipeline/deduplicated.py \
# /home/hieu/src/sc_bulk_BCR_data_analysis/GCTree_pipeline/work_${mouseID}_${inputCase};
# rm -rf /home/hieu/src/sc_bulk_BCR_data_analysis/GCTree_pipeline/work_${mouseID}_${inputCase};
# done;done