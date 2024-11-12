path_to_cellranger="/home/hieunguyen/CRC1382/src_2023/CellRanger/cellranger-7.1.0"
export PATH=${path_to_cellranger}:$PATH
for sampleid in PP7;do cellranger multi --id=${sampleid} --csv=config.csv --localcores=14;done

