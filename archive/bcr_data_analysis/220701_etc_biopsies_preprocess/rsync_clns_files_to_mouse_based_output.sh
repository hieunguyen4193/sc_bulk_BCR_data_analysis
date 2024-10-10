##### Example input args
# mouse_based_list="/home/hieu/src/BCRTree_release/mouse_based_MID_list";
# outputdir="/home/hieu/outdir/mixcr_pipeline_output/mouse_based_output";
# mid_based_output="/home/hieu/outdir/mixcr_pipeline_output/mid_based_output";

##### explanation
# mid_based_output: path to the input MID-based outputs obtained from mixcr_pipeline.
# outputdir: path to save output
# group
while getopts "i:o:t:r:" opt; do
  case ${opt} in
    i )
      mid_based_output=$OPTARG
      ;;
    o )
      outputdir=$OPTARG
      ;;
    l )
      mouse_based_list=$OPTARG
      ;;
    n )
      group_name=$OPTARG
      ;;
    \? )
      echo "Usage: cmd [-i] mid_based_output [-o] outputdir [-l] mouse_based_list [-n] group_name"
      exit 1
      ;;
  esac
done

outputdir=${outputdir}/${group_name}
mkdir -p ${outputdir}
files=$(ls ${mouse_based_list});

for file in $files;do \
list=$(cat ${mouse_based_list}/$file | cut -f1) \
&& for sampleid in $list;do \
mkdir -p ${outputdir}/${file%.csv*} && \
rsync -ahv --progress ${mid_based_output}/${sampleid}/*.clns \
${outputdir}/${file%.csv*} --exclude=*.reassigned*;done;done
