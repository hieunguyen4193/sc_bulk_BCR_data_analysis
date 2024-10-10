mouse_based_list="/home/hieu/src/BCRTree_release/mouse_based_MID_list";
outputdir="/home/hieu/outdir/mixcr_pipeline_output/mouse_based_output";
mid_based_output="/home/hieu/outdir/mixcr_pipeline_output/mid_based_output";
mkdir -p ${outputdir}
files=$(ls ${mouse_based_list});

for file in $files;do \
list=$(cat ${mouse_based_list}/$file | cut -f1) \
&& for sampleid in $list;do \
mkdir -p ${outputdir}/${file%.csv*} && \
rsync -ahv --progress ${mid_based_output}/${sampleid}/*.clns \
${outputdir}/${file%.csv*} --exclude=*.reassigned*;done;done
