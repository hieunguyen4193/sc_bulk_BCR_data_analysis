#> Download the germline FASTA files for generating input to GCTREE.
#> link: https://www.imgt.org/IMGTindex/Fasta.php
#> also see immrep/code/gctree.single.sh (Fabio's code)
#> or should we use https://www.imgt.org/vquest/refseqh.html IMGT/GENE-DB reference directory sets

  
library("Biostrings")
library(stringr)
# path.to.fasta <- "/home/hieu/src/BCRTree_release/gctree/data_analysis"
path.to.fasta <- "/home/hieunguyen/CRC1382/src_2023/bcrtree/data_analysis"

s.V.genes <- readDNAStringSet(file.path(path.to.fasta, "IGHV.fasta")) 
names(s.V.genes) <- lapply(names(s.V.genes), function(x){
  str_split(x, "[|]")[[1]][[2]]
})

# s.D.genes <- readDNAStringSet(file.path(path.to.fasta, "IGHD.fasta")) 
# names(s.D.genes) <- lapply(names(s.D.genes), function(x){
#   str_split(x, "[|]")[[1]][[2]]
# })

s.J.genes <- readDNAStringSet(file.path(path.to.fasta, "IGHJ.fasta")) 
names(s.J.genes) <- lapply(names(s.J.genes), function(x){
  str_split(x, "[|]")[[1]][[2]]
})
