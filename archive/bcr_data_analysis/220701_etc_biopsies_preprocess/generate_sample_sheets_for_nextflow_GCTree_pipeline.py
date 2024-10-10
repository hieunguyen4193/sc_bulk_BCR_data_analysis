import pandas as pd
import pathlib
import os

'''
Generate Sample Sheet which contains clones and their paths to the FASTQ files.
Clones are sorted in ascending order of clone sizes. 
'''

outdir = "/home/hieunguyen/CRC1382/outdir/molmed_server"
PROJECT = "mixcr_pipeline_output"

maindir = os.path.join(outdir, PROJECT)
path_to_main_output = os.path.join(maindir, "data_analysis")
path_to_05_output = os.path.join(path_to_main_output, "05_output")

thres = 0.15
maindf = pd.DataFrame()
all_count_files = [item for item in pathlib.Path(os.path.join(path_to_05_output)).glob("*/*/*/*count*")]
for file in all_count_files:
    basedir = "/".join(str(file).split("/")[:-1])
    tmpdf = pd.read_csv(file, index_col = [0])
    mouseid = str(file).split("/")[-3]
    analysis_type = str(file).split("/")[-2].replace("_MIDs", "")
    tmpdf["mouseID"] = mouseid
    tmpdf["analysis_type"] = analysis_type
    tmpdf["filename"] = tmpdf["VJ.combi"].apply(lambda x: x.replace("*", "-"))
    tmpdf["path"] = tmpdf["filename"].apply(lambda x: "{}/{}_{}_{}.aln.fasta".format(basedir, mouseid, analysis_type, x))
    
    tmpdf["filename"] = tmpdf[["mouseID", "analysis_type", "filename"]].apply(lambda x: "_".join(x), axis = 1)
    maindf = pd.concat([maindf, tmpdf], axis = 0)

maindf = maindf.sort_values(by = ["count.seq"], ascending=True)
maindf.to_csv("all_fasta_SampleSheet.csv", index = False)

count = 0
for mouseid in maindf["mouseID"].unique():
    for analysis_type in maindf["analysis_type"].unique():
        tmpdf = maindf[(maindf["mouseID"] == mouseid) & (maindf["analysis_type"] == analysis_type)]
        tmpdf.to_csv("./SampleSheets/{}_{}_SampleSheet.csv".format(mouseid, analysis_type), index = False)
        count += tmpdf.shape[0] 
