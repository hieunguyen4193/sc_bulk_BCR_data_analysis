# Trong-Hieu Nguyen, last modified on 24.08.2022.

current_dir=$(pwd)

path_to_cellranger="/home/hieunguyen/CRC1382/src/src_pipeline/CellRanger/cellranger-6.1.2"
export PATH=${path_to_cellranger}:$PATH

path_to_save_output="/media/hieunguyen/HD0/ext_HDD/OUTPUT/BSimons_datasets/220504_Simons_Pabst_MolMedizin_scVDJseq/rerun_CellRanger_with_YFP_no_linebreak_fa"
mkdir -p ${path_to_save_output}

##############################################################################
# Write logs for this run
logfile=${path_to_save_output}/logs.txt
echo -e "###################################################################" >> ${logfile}
echo -e "CellRanger version: " ${path_to_cellranger} >> ${logfile}
echo -e "Cell ranger configs: " >> ${logfile}
echo -e "----- ----- ----- ----- -----" >> ${logfile}
echo -e "Sample 132" >> ${logfile}
cat multi_config_132.csv >> ${logfile}
echo -e "Sample 133" >> ${logfile}
cat multi_config_133.csv >> ${logfile}
echo "----- ----- ----- ----- -----" >> ${logfile}

echo -e "Command to download INPUT FASTQ files: wget -r --user simons --password 760QRxUgr0bH https://genomics.rwth-aachen.de/data/220504_Simons_Pabst_MolMedizin_scVDJseq/1_Raw_data/fastq/"  >> ${logfile}
echo -e "Path to outputs: " ${path_to_save_output} >> ${logfile}

echo -e "\n"

echo -e "Notes: In this analysis, we run the CellRanger multi pipeline for the two samples from B. Simons's second dataset. A YFP gene was added to the MM10 reference genome, see multi_config.csv files for paths to these modified references." >> ${logfile}

echo -e "Outputs are temporarily stored at " ${current_dir} "and then move to " ${path_to_save_output} >> ${logfile}

echo -e "\n"
echo -e "###################################################################" >> ${logfile}

##############################################################################

mkdir -p ${path_to_save_output}/BSimons_sample_132_added_YFP
mkdir -p ${path_to_save_output}/BSimons_sample_133_added_YFP

# M A I N - C E L L R A N G E R - C O M M A N D S
cellranger multi --id=BSimons_sample_132_added_YFP --csv=multi_config_132.csv --localcores=25

rsync -avh --progress ${current_dir}/BSimons_sample_132_added_YFP ${path_to_save_output}/BSimons_sample_132_added_YFP 

rm -rf ${current_dir}/BSimons_sample_132_added_YFP

cellranger multi --id=BSimons_sample_133_added_YFP --csv=multi_config_133.csv --localcores=25
rsync -avh --progress ${current_dir}/BSimons_sample_133_added_YFP ${path_to_save_output}/BSimons_sample_133_added_YFP 

rm -rf ${current_dir}/BSimons_sample_133_added_YFP

# EOF

