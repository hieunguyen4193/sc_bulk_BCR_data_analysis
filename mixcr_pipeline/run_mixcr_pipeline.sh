##### batch 1: 220701_etc_biopsies
# path_to_input="/mnt/storage/data/local/mol_med/bcr/220701_etc_biopsies/samples";
##### batch 2:
# path_to_input="/home/hieu/storage/240826_Simons_pabst_MolMedcine_AmpliSeq/Fastq";
#####: 240826_BSimons dataset
##### batch 3:
# path_to_input="/home/hieu/storage/241031_Simons_Pabst_MolMedicine_FASTQ";

while getopts "i:o:e:" opt; do
  case ${opt} in
    i )
      path_to_input=$OPTARG
      ;;
    o )
      outputdir=$OPTARG
      ;;
    e )
      ext=$OPTARG
      ;;
    \? )
      echo "Usage: cmd [-i] path_to_input [-o] outputdir [-e] ext"
      exit 1
      ;;
  esac
done


mkdir -p $outputdir;

files=$(ls ${path_to_input}/*R1*.fastq*);
for file in $files;do \
filename=$(echo $file | xargs -n 1 basename);
sampleid=${filename%_S*};
fastq1=${file}
fastq2=${file%R1*}R2_001.${ext}; # add .gz depending on the input fastq format ext
echo -e "------------------------------------------"
echo -e "WORKING ON SAMPLE " $filename "\n"
echo -e "------------------------------------------"
echo -e $fastq1;
echo -e $fastq2;

bash mixcr_pipeline.sh -i $sampleid -o $outputdir -f $fastq1 -r $fastq2;
done