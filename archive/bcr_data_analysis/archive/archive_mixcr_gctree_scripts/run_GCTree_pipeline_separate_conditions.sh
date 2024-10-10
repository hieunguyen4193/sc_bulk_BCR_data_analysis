outdir="/home/hieunguyen/CRC1382/outdir/molmed_server/gctree";
output_version="20240723";
project_src="/home/hieunguyen/CRC1382/src_2023/bcrtree";

output=${outdir}/${output_version};
mkdir -p ${output};

# samplesheet="/home/hieunguyen/CRC1382/src_2023/bcrtree/SampleSheets/m14_all_YFP_SampleSheet.csv"

all_sample_sheets=$(ls ./SampleSheets/*.csv);

samplesheet="/home/hieunguyen/CRC1382/src_2023/bcrtree/all_fasta_SampleSheet.csv"
samplesheet_name=$(echo $samplesheet | xargs -n 1 basename);

deduplicate_src=${project_src}/deduplicated.py;
modify_tree_colors=${project_src}/modify_tree_colors.py;
color_path=${project_src}/hex_color.csv;

work=${outdir}/gctree_pipeline_work_${samplesheet_name%.csv*};

nextflow run GCtree_pipeline_input_SampleSheet.nf \
--samplesheet ${samplesheet} \
--output ${output}/${samplesheet_name%.csv*} \
--deduplicate_src ${deduplicate_src} \
--modify_tree_colors ${modify_tree_colors} \
--color_path ${color_path} -resume -w ${work};

rm -rf $work;