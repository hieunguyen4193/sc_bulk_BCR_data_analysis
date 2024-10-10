path_to_input="/mnt/storage/data/local/mol_med/bcr/220701_etc_biopsies/samples";

files=$(ls ${path_to_input}/*R1*.fastq*);
for file in $files;do \
filename=$(echo $file | xargs -n 1 basename);
sampleid=${filename%_S*};
fastq1=${file}
fastq2=${file%R1*}R2_001.fastq;
echo -e "------------------------------------------"
echo -e "WORKING ON SAMPLE " $sampleid "\n"
echo -e "------------------------------------------"
echo $sampleid;
echo -e $fastq1;
echo -e $fastq2;

bash mixcr_pipeline.sh $sampleid $fastq1 $fastq2;
done

echo -e  "Sync and group MID data into mouse based data";
bash rsync_clns_files_to_mouse_based_output.sh;
echo -e "finished syncing data to moused based data";

echo -e "generating moused based trees ..."
bash generate_mouse_based_trees.sh;

