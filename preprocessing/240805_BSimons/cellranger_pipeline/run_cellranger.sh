path_to_cellranger="/home/hieunguyen/CRC1382/src_2023/CellRanger/cellranger-7.1.0"
export PATH=${path_to_cellranger}:$PATH
for sampleid in M2 M3 P1 P2 P3;do cellranger multi --id=${sampleid} --csv=multi_config_${sampleid}.csv --localcores=14;
rsync -ahv --progress ${sampleid} /media/hieunguyen/HD01/240805_BSimons/cellranger_output && rm -rf ${sampleid};done
