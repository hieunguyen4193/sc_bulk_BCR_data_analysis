import pandas as pd
import numpy as np
import pathlib
import matplotlib.pyplot as plt
import seaborn as sns
import os
from tqdm import tqdm
from typing import List, Union, Optional, Callable
import pickle
from Bio import AlignIO, SeqIO
from ete3 import Tree, TreeNode
from gctree import CollapsedTree

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import umap
from ete3 import Tree, faces, TreeStyle, NodeStyle, TextFace, SequenceFace, COLOR_SCHEMES, CircleFace
from GCTree_preparation import *
import warnings
import math
import nltk
from matplotlib.patches import Rectangle
warnings.filterwarnings("ignore")

path_to_storage = "/media/hieunguyen/HNHD01/storage/all_BSimons_datasets"
outdir = "/media/hieunguyen/GSHD_HN01/outdir/sc_bulk_BCR_data_analysis_v0.1"
mouseid = "m3"

bulk_project = "241031_BSimons"
sc_project = ["240411_BSimons", "241002_BSimons"]

PROJECTS = f"{bulk_project}_{'_'.join(sc_project)}"

path_to_01_output = f"{outdir}/tree_analysis/{bulk_project}/01_output"
path_to_07_output = f"{outdir}/tree_analysis/07_output/{PROJECTS}/{mouseid}"
os.system(f"mkdir -p {path_to_07_output}")

with open(f"{path_to_01_output}/saveTreeobj.pkl", "rb") as f:
    saveTreeobj = pickle.load(f)

bulk_metadata = pd.read_excel("/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis/preprocessing/241031_BSimons/241031_sample_sheet.xlsx")
bulk_metadata["MID"] = bulk_metadata["MID"].apply(lambda x: f"MID{x}")
bulk_metadata.columns = ["MID", "mouseID", "organ", "YFP", "population"]

maindf = pd.read_csv(f"{path_to_01_output}/tree_summarydf.csv")
maindf = maindf[maindf["mouseID"] == mouseid]

path_to_04_output = os.path.join(outdir, "VDJ_output", "04_output")
thres = 0.85
clonedf = pd.read_csv(os.path.join(path_to_04_output, "full_clonedf_with_mutation_rate.csv"), index_col= [0])
clonedf = clonedf[clonedf['num_mutation'] != "region_not_covered-skip"]

sc_clonedf = clonedf[clonedf['dataset.name'].isin(sc_project)][["barcode", "id", 'V.gene', 'J.gene', 'D.gene', "aaSeqCDR3", "nSeqCDR3"]].reset_index().drop("index", axis = 1)

sc_clonedf["mouseID"] = sc_clonedf["id"].apply(lambda x: "m" + x.replace("M", "").replace("P", ""))
sc_clonedf = sc_clonedf[sc_clonedf["mouseID"] == mouseid]

color_path = f"{bulk_project}_color.csv"
tree_dist_seqdf = pd.DataFrame()

for cloneid in tqdm(maindf.cloneid.unique()):
    os.system(f"mkdir -p {path_to_07_output}/{cloneid}")
    V_gene = cloneid.split("_")[1]
    J_gene = cloneid.split("_")[2]
    cdr3_len = cloneid.split("_")[3]
    mouseid = cloneid.split("_")[0]

    treeobj = saveTreeobj[cloneid]
    seqdf = treeobj.seqdf.copy()
    idmapdf = treeobj.idmapseqdf.copy()
    bulkdf = seqdf.merge(idmapdf, right_on = "seq", left_on = "seq")
    bulkdf["aaSeqCDR3"] = bulkdf["ID"].apply(lambda x: x.split("|")[2].split(":")[1])
    
    scdf = sc_clonedf[(sc_clonedf["V.gene"] == V_gene.split("-0")[0]) & 
                    (sc_clonedf["J.gene"] == J_gene.split("-0")[0])]

    scdf["len"] = scdf["nSeqCDR3"].apply(lambda x: len(x))
    scdf = scdf[scdf["len"] == int(cdr3_len)]

    if scdf.shape[0] > 0:
        for i in range(0, scdf.shape[0]):
            barcode = scdf.iloc[i].barcode
            sampleid = scdf.iloc[i].id
            sc_seq = scdf[(scdf["barcode"] == barcode) & (scdf["id"] == sampleid)]["aaSeqCDR3"].values[0] 

            bulkdf[f"{sampleid}_{barcode}"] = bulkdf["aaSeqCDR3"].apply(lambda x: nltk.edit_distance(x, sc_seq)/len(x))
        bulkdf["min_dist_to_a_cell"] = bulkdf[[item for item in bulkdf.columns if sampleid in item]].apply(lambda x: min(x), axis = 1)
        bulkdf["cloneid"] = cloneid
        bulkdf.to_csv(os.path.join(path_to_07_output, cloneid, f"tree_Seqdf_{cloneid}.csv"))
        bulkdf_simplified = bulkdf[["ID", "seq", "abundance", "MID", "seqid", "min_dist_to_a_cell"]]
        
        tree_dist_seqdf = pd.concat([tree_dist_seqdf, bulkdf_simplified], axis = 0)
        
        ##### Generate heatmap and save

        heatmap_plotdf = bulkdf[["seqid"] + [item for item in bulkdf.columns if "_" in item and item != "min_dist_to_a_cell"]].set_index("seqid")
        plt.figure(figsize= (10, 10))
        sns.heatmap(heatmap_plotdf,
                    cmap="coolwarm", cbar=-1, linewidths=0.5, linecolor='black')
        for tick_label in plt.gca().get_xticklabels():
            tick_text = tick_label.get_text()
            if "M" in tick_text:
                tick_label.set_color('blue')
            elif "P" in tick_text:
                tick_label.set_color('red')
        plt.tight_layout()
        plt.savefig(os.path.join(path_to_07_output, cloneid, f"heatmap_{cloneid}.svg"))
        plt.close()

        plt.figure(figsize= (10, 10))
        coordinates = list(zip(*np.where(heatmap_plotdf.values == 0)))
        ax = sns.heatmap(heatmap_plotdf,
                    cmap="coolwarm", cbar=-1, linewidths=0.5, linecolor='black')
        for tick_label in plt.gca().get_xticklabels():
            tick_text = tick_label.get_text()
            if "M" in tick_text:
                tick_label.set_color('blue')
            elif "P" in tick_text:
                tick_label.set_color('red')

        for i in range(len(coordinates)):
            ax.add_patch(Rectangle((coordinates[i][1], coordinates[i][0]),1,1, fill=False, edgecolor='red', lw=3))
        plt.tight_layout()
        plt.savefig(os.path.join(path_to_07_output, cloneid, f"heatmap_{cloneid}.annotated.svg"))
        plt.close()

        ##### generate tree and save
        treeobj = saveTreeobj[cloneid] 
        avai_mids = treeobj.seqdf["MID"].unique()
        mid_color_pal = pd.read_csv(color_path, index_col = [0]).to_dict()["hex color"]

        ts = treeobj.generate_tree_style(color_path = color_path)

        for input_mid in avai_mids:
            if input_mid == "GL":
                input_mid_col = "gray"
            else:
                input_mid_col = mid_color_pal[input_mid]
            ts.legend.add_face(CircleFace(10, input_mid_col), column = 0)
            ts.legend.add_face(TextFace(bulk_metadata[bulk_metadata["MID"]==input_mid]["population"].values[0]), column = 0)

        treeobj.tree.render(f"{path_to_07_output}/{cloneid}/{cloneid}.tree.svg", tree_style = ts)
    else:
        bulkdf.to_csv(os.path.join(path_to_07_output, f"tree_Seqdf_{cloneid}.csv"))

    scdistdf = pd.DataFrame()
for cloneid in tqdm(maindf.cloneid.unique()):
    V_gene = cloneid.split("_")[1]
    J_gene = cloneid.split("_")[2]
    cdr3_len = cloneid.split("_")[3]
    mouseid = cloneid.split("_")[0]

    treeobj = saveTreeobj[cloneid]
    seqdf = treeobj.seqdf.copy()
    idmapdf = treeobj.idmapseqdf.copy()
    bulkdf = seqdf.merge(idmapdf, right_on = "seq", left_on = "seq")
    bulkdf["aaSeqCDR3"] = bulkdf["ID"].apply(lambda x: x.split("|")[2].split(":")[1])

    scdf = sc_clonedf[(sc_clonedf["V.gene"] == V_gene.split("-0")[0]) & 
                    (sc_clonedf["J.gene"] == J_gene.split("-0")[0])]

    scdf["len"] = scdf["nSeqCDR3"].apply(lambda x: len(x))
    scdf = scdf[scdf["len"] == int(cdr3_len)]

    all_barcodes = [f"{scdf['id'].values[i]}_{scdf['barcode'].values[i]}" for i in range(scdf.shape[0])]

    if scdf.shape[0] > 0:
        for i in range(0, scdf.shape[0]):
            barcode = scdf.iloc[i].barcode
            sampleid = scdf.iloc[i].id
            sc_seq = scdf[(scdf["barcode"] == barcode) & (scdf["id"] == sampleid)]["aaSeqCDR3"].values[0] 

            bulkdf[f"{sampleid}_{barcode}"] = bulkdf["aaSeqCDR3"].apply(lambda x: nltk.edit_distance(x, sc_seq)/len(x))
            
    tmp_scdistdf = bulkdf[[item for item in bulkdf.columns if item in all_barcodes]].min().reset_index()
    tmp_scdistdf.columns = ["barcode", "min_dist_to_a_tree"]
    tmp_scdistdf["treeID"] = cloneid
    scdistdf = pd.concat([scdistdf, tmp_scdistdf], axis = 0)
scdistdf.to_csv(os.path.join(path_to_07_output, "scdistdf.csv"))