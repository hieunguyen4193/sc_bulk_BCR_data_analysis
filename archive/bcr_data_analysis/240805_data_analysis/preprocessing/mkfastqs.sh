path_to_cellranger="/home/hieunguyen/CRC1382/src_2023/CellRanger/cellranger-7.1.0"
export PATH=${path_to_cellranger}:$PATH
sample_sheet="/home/hieunguyen/CRC1382/src_2023/240805_BSimons/SampleSheet.csv";
inputdir="/home/hieunguyen/CRC1382/outdir/240805_BSimons/tar"
cellranger mkfastq --run ${inputdir} --sample-sheet ${sample_sheet};


