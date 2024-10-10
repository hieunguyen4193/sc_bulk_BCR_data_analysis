import pandas as pd
from collections import defaultdict
import os
import sys
import argparse

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from PIL import Image

# headless, no screen, avoid the erorr
# qt.qpa.xcb: could not connect to display 
# qt.qpa.plugin: Could not load the Qt platform plugin "xcb" in "" even though it was found.
# This application failed to start because no Qt platform plugin could be initialized. Reinstalling the application may fix this problem.
# Available platform plugins are: eglfs, linuxfb, minimal, minimalegl, offscreen, vnc, wayland-egl, wayland, wayland-xcomposite-egl, wayland-xcomposite-glx, webgl, xcb.
# Aborted (core dumped)
os.environ['QT_QPA_PLATFORM'] = 'offscreen'
import math
from ete3 import Tree, faces, TreeStyle, NodeStyle, TextFace, SequenceFace, COLOR_SCHEMES
from Bio import AlignIO, SeqIO
import pickle
import numpy as np

def get_parser(argv):
    parser = argparse.ArgumentParser(
        description=(
            " Create a pie chart like node in the GCtree."
        )
    )
    parser.add_argument(
        "--input_fasta",
        type=str,
        help="The original input fasta file. Sequence ID must contain full information of the sequence",
    )
    parser.add_argument(
        "--input_idmap",
        type=str,
        help="ID MAP file for mapping between seq IDs and real sequences, generated from the function deduplicate",
    )
    parser.add_argument(
        "--gctree_inference_file",
        type=str,
        help="The file *gctree.out.inference.1.p output from gctree inference",
    )
    parser.add_argument(
        "--color_path",
        type=str,
        help="Path to color pallete file for each MID",
    )
    parser.add_argument(
            "--output",
            type=str,
            help="Output dir",
        )
    parser.add_argument(
            "--svg_name",
            type=str,
            help="Tree name",
        )

    return parser.parse_args(argv)


def main(args):
    # idmapseqdf = pd.read_csv("/home/hieu/src/BCRTree_release/gctree/example/output/m11_IGHV1-26*01_IGHJ2*01_45/m11_IGHV1-26*01_IGHJ2*01_45.id_map_seq.csv")
    # path_to_orig_fasta = "/home/hieu/src/BCRTree_release/gctree/example/m11_IGHV1-26*01_IGHJ2*01_45.aln.fasta"
    # color_path = "/home/hieu/outdir/mixcr_pipeline_output/data_analysis/01_output/mid_color_pal.csv"
    # path_to_gctree_inference = "/home/hieu/src/BCRTree_release/gctree/example/output/m11_IGHV1-26*01_IGHJ2*01_45/gctree.out.inference.1.p"
    
    idmapseqdf = pd.read_csv(args.input_idmap)
    path_to_orig_fasta = args.input_fasta
    path_to_gctree_inference = args.gctree_inference_file
    color_path = args.color_path
    outputdir = args.output
    svg_name = args.svg_name
    
    idmapseqdf.columns = ["seqid", "seq"]
    with open(path_to_orig_fasta) as fasta_file:  # Will close handle cleanly
        identifiers = []
        seqs = []
        for seq_record in SeqIO.parse(fasta_file, 'fasta'):  # (generator)
            identifiers.append(seq_record.id)
            seqs.append(str(seq_record.seq))
    
    seqdf = pd.DataFrame(data = identifiers, columns = ["ID"])
    seqdf["seq"] = seqs
    seqdf = seqdf[seqdf["ID"] != "GL"]
    seqdf["abundance"] = seqdf["ID"].apply(lambda x: int(x.split("|")[-1].replace("Abundance:", "")))
    seqdf["MID"] = seqdf["ID"].apply(lambda x: str(x.split("|")[0].replace("Sample:", "")))
    seqdf_orig = seqdf.copy()
    seqdf = seqdf.groupby("seq")["abundance"].sum().reset_index().copy()

    # fix color list, MID color.
    avai_mids = seqdf_orig.MID.unique()
    mid_color_pal = pd.read_csv(color_path, index_col = [0]).to_dict()["hex color"]

    # assert(seqdf.shape[0] == idmapseqdf.shape[0])
    # assert(len([item for item in seqdf.seq.values if item in idmapseqdf.seq.values]) == idmapseqdf.shape[0])
    
    seqdf = seqdf.merge(idmapseqdf, right_on = "seq", left_on = "seq")
    
    abund_pct = dict()
    for node_name in seqdf.seqid.unique():
        seq = seqdf[seqdf["seqid"] == node_name].seq.unique()[0]
        total_abund = seqdf[seqdf["seqid"] == node_name].abundance.unique()[0]
        abund_pct[node_name] = dict()
        tmpdf = seqdf_orig[seqdf_orig["seq"] == seq]
        for mid in tmpdf["MID"].unique():
            abund_pct[node_name][mid] = 100*np.sum(tmpdf[(tmpdf["MID"] == mid)]["abundance"].values)/total_abund
    
    def layout(n):
        size = max(1, 10 * math.sqrt(n.abundance))
    
        if n.abundance > 1:
            cols = [mid_color_pal[mid] for mid in abund_pct[n.name].keys()]
            values = [abund_pct[n.name][mid] for mid in abund_pct[n.name].keys()]
            F = faces.PieChartFace(values, colors=cols,
                                   width=size * 2, height=size * 2)
            F.border.width = None
            # F.opacity = 0.6
            faces.add_face_to_node(F, n, 0, position="branch-right")
            ns = NodeStyle()
            ns["size"] = 0
            n.set_style(ns)
    
    with open(path_to_gctree_inference, "rb") as fd:
        p = pickle.load(fd)
    t = p.tree
    
    for n in t.traverse("postorder"):
        if n.abundance > 1:
            cols = [mid_color_pal[mid] for mid in abund_pct[n.name].keys()]
            values = [abund_pct[n.name][mid] for mid in abund_pct[n.name].keys()]
    
    ts = TreeStyle()
    ts.layout_fn = layout
    ts.mode = "r"
    ts.rotation = 90
    ts.show_leaf_name = False
    # t.show(tree_style=ts)
    t.render(os.path.join(outputdir, "{}.png".format(svg_name)), w=1280, tree_style=ts)

    ##### Add legend
    img = mpimg.imread(os.path.join(outputdir, "{}.png".format(svg_name)))
    
    # Create a figure
    fig, ax = plt.subplots()
    ax.imshow(img)
    ax.axis('off')  # Hide axes
    
    # Add legend manually
    legend_labels = avai_mids
    legend_colors = [mid_color_pal[mid] for mid in avai_mids]
    
    # Create legend patches
    import matplotlib.patches as mpatches
    patches = [mpatches.Patch(color=color, label=label) for color, label in zip(legend_colors, legend_labels)]
    plt.legend(handles=patches, loc='lower right')
    
    # Save the final image
    plt.savefig(os.path.join(outputdir, "{}.withLegend.svg".format(svg_name)), bbox_inches='tight', pad_inches=0.1, dpi=300, format="svg")
    plt.show()
    
if __name__ == '__main__':
    main(get_parser(sys.argv[1:]))
