######-------------------------------------------------------------------------------------######
##### Run the full pipeline using GCTREE to generate phylogenetics tree for BCR sequences
#####-------------------------------------------------------------------------------------######

# Draft; bash pipeline.sh /home/hieu/src/BCRTree_release/gctree/example/m11_IGHV1-26\*01_IGHJ2\*01_45.aln.fasta /home/hieu/src/BCRTree_release/gctree/example/output/
thres=0.15;
mixcrdir=/home/hieu/outdir/mixcr_pipeline_output/data_analysis/02_output/CDR3_${thres};
all_fasta=$(ls ${mixcrdir}/m*/*.fasta | xargs -n 1 basename);
pipeline_outputdir="/home/hieu/outdir/gctree_output";

for f in ${all_fasta};do \
mouseID=$(echo $f | cut -d'_' -f1);
OUTPUT_FILE=${pipeline_outputdir}/${mouseID}/${f%.aln.fasta*}/${f%.aln.fasta*}.color.withLegend.svg;
echo $OUTPUT_FILE
if [ -f "$OUTPUT_FILE" ]; then
echo $OUTPUT_FILE exists;
else
mkdir -p ${pipeline_outputdir}/${mouseID};
echo -e "Working on mouse " $mouseID "and sample " $f;
bash pipeline_single_fasta.sh ${mixcrdir}/${mouseID}/${f} ${pipeline_outputdir}/${mouseID};
fi;
done
