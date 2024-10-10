# General code-repo for the BCR-bulk and single-cell data analysis.

## Description

- `first_two_scRNA_datasets`: analysis scripts for the first two batches of scRNA + BCR datasets, including `220504_Simons_Pabst_MolMedizin_scVDJseq`, `211216_BSimons`, `211221_BSimons`, `220104_BSimons`. The sub-folder `run_CellRanger_Bonn_data` contains preprocessing scripts and `CellRanger` pipeline for the external datasets (Bonn) `BSimons_Bonn_data`. Raw data and `CellRanger` outputs are stored in NAS at `/volume1/CRC1382_NAS01/CRC1382/OFFICIAL/raw_data`. Note that, for the dataset `220504_Simons_Pabst_MolMedizin_scVDJseq`, we re-run the `CellRanger` pipeline and add **YFP** gene sequence to the mm10 reference genome during alignment steps. 

- `220701_etc_biopsies_data_analysis`: data analysis scripts including `mixcr` pipeline and `gctree` pipeline for the old bulk-BCR dataset. 

- `240805_data_analysis`: data analysis scripts for the new bulk-BCR dataset (240805_BSimons). 

- `240826_data_analysis`: data analysis scripts for the new single-cell BCR VDJ dataset (240826_BSimons).

All raw data are stored in NAS and Coscine RWTH. 

## Data analysis

Data analysis scripts for the output from `mixcr` pipeline are stored in `mixcr_data_analysis` folder. These scripts will take data from `220701_etc_biopsies_preprocess` and `240826_data_preprocess` as inputs.

Scripts:
- `01_clone_summary_and_prepare_fasta_for_GCTree.R`: due to some unknown error/bug/whatever reason it is, this script throws some confusing error when being run in `Rstudio`. The script runs perfectly fine in command line inside the `docker` container `tronghieunguyen/scrna_gex_pipeline_ibex_trex_msa`. This generates the output folder `01_output` which contains a sub-folder `CDR3_<BCR sequence dissimilarity threshold>` (currently use 0.15). `Fasta` files are stored in `<mouse ID>` folder. Summary tables for all clones are saved in `.xlsx` files. 

