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
# /home/hieu/outdir/${mouseID}_${inputCase} \
# /home/hieu/src/sc_bulk_BCR_data_analysis/GCTree_pipeline/deduplicated.py \
# /home/hieu/src/sc_bulk_BCR_data_analysis/GCTree_pipeline/work_${mouseID}_${inputCase};
# rm -rf /home/hieu/src/sc_bulk_BCR_data_analysis/GCTree_pipeline/work_${mouseID}_${inputCase};
# done;done

# sudo bash run_docker_GCTree.sh 8181 rstudio1 genov4
# mouseID="m3";
# mouseID="m7";
# bash run_nf.sh \
# /home/hieu/src/241031_${mouseID}/sc_bulk_BCR_data_analysis/SampleSheet_GCTree_nextflow_241031_BSimons/SampleSheet_GCTree_241031_BSimons_${mouseID}.nextflow.csv \
# /home/hieu/outdir/241031_BSimons_${mouseID} \
# /home/hieu/src/241031_${mouseID}/sc_bulk_BCR_data_analysis/GCTree_pipeline/deduplicated.py \
# /home/hieu/src/241031_${mouseID}/sc_bulk_BCR_data_analysis/work_${mouseID};
# rm -rf /home/hieu/src/241031_${mouseID}/sc_bulk_BCR_data_analysis/work_${mouseID};

# inputCase="all";
# mouseID="m2";
# bash run_nf.sh \
# /home/hieu/src/${mouseID}_all/sc_bulk_BCR_data_analysis/SampleSheet_GCTree_nextflow_240826_BSimons/SampleSheet_GCTree_240826_BSimons_${mouseID}_${inputCase}.nextflow.csv \
# /home/hieu/outdir/${mouseID}_${inputCase} \
# /home/hieu/src/${mouseID}_all/sc_bulk_BCR_data_analysis/GCTree_pipeline/deduplicated.py \
# /home/hieu/src/${mouseID}_all/sc_bulk_BCR_data_analysis/GCTree_pipeline/work_${mouseID}_${inputCase};
# rm -rf /home/hieu/src/${mouseID}_all/sc_bulk_BCR_data_analysis/GCTree_pipeline/work_${mouseID}_${inputCase};
