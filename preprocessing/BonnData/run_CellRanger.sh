# Trong-Hieu Nguyen, last modified on 24.08.2022.

current_dir=$(pwd)

path_to_cellranger="/work/uk104163/CRC1382/cellranger/cellranger-7.1.0"
export PATH=${path_to_cellranger}:$PATH



##############################################################################
# M A I N - C E L L R A N G E R - C O M M A N D S
cellranger multi --id=GF_7w_14w --csv=multi_config.csv --localcores=40
# EOF


