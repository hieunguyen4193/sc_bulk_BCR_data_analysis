path_to_cellranger="/home/hieunguyen/CRC1382/src_2023/CellRanger/cellranger-7.1.0"
export PATH=${path_to_cellranger}:$PATH
sample_sheet="/home/hieunguyen/CRC1382/outdir/compressed_tar/SampleSheet.csv";
inputdir="/home/hieunguyen/CRC1382/outdir/compressed_tar/241002_A01742_0300_BHTL7MDRX3"
cellranger mkfastq --run ${inputdir} --sample-sheet ${sample_sheet};


