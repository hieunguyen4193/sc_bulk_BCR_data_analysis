#> Download the germline FASTA files for generating input to GCTREE.
#> link: https://www.imgt.org/IMGTindex/Fasta.php
#> also see immrep/code/gctree.single.sh (Fabio's code)
#> or should we use https://www.imgt.org/vquest/refseqh.html IMGT/GENE-DB reference directory sets

library("Biostrings")
library(stringr)
# path.to.fasta <- "/home/hieu/src/BCRTree_release/gctree/data_analysis"
path.to.fasta <- "/media/hieunguyen/HNSD01/src/bcr_data_analysis/1st_2nd_datasets"

s.V.genes <- readDNAStringSet(file.path(path.to.fasta, "IGHV.fasta")) 
names(s.V.genes) <- lapply(names(s.V.genes), function(x){
  str_split(x, "[|]")[[1]][[2]]
})

s.J.genes <- readDNAStringSet(file.path(path.to.fasta, "IGHJ.fasta")) 
names(s.J.genes) <- lapply(names(s.J.genes), function(x){
  str_split(x, "[|]")[[1]][[2]]
})

cellranger.ref <- readDNAStringSet("/media/hieunguyen/HNSD01/src/bcr_data_analysis/refdata-cellranger-vdj-GRCm38-alts-ensembl-7.0.0/fasta/regions.fa")

