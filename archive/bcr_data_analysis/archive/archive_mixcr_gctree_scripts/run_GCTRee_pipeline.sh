input="/home/hieu/outdir/mixcr_pipeline_output/data_analysis/02_output/CDR3_0.15";
# input="/home/hieu/src/BCRTree_release/gctree/example";
output="/home/hieu/outdir/gctree_output/nextflow_output";
deduplicate_src="/home/hieu/src/BCRTree_release/gctree/deduplicated.py";
modify_tree_colors="/home/hieu/src/BCRTree_release/gctree/modify_tree_colors.py";
color_path="/home/hieu/src/BCRTree_release/gctree/hex_color.csv";
work=$output/work;
nextflow run GCtree_pipeline.nf \
--input $input \
--output $output \
--deduplicate_src $deduplicate_src \
--modify_tree_colors $modify_tree_colors \
--color_path $color_path \
--file_pattern "m*/*.fasta" \
-resume -w $work
