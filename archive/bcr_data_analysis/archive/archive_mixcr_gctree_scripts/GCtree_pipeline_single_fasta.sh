#First we deduplicate the sequences and convert to phylip alignment format,
# and also determine sequence abundances. The deduplicate command writes
# the phylip alignment of unique sequences to stdout (which we redirect to
# a file here). The argument --root indicates the root id.
# The flag --id_abundances can be used to indicate that integer sequence ids
# should be interepreted as abundances. The argument --abundance_file
# indicates that sequence abundance information should be written to the
# specified csv file. The argument --idmapfile allows us to specify
# a filename for the output file containing a mapping of new, unique sequence
# IDs to original sequence IDs from the fasta file.

input_fasta=$1;
outputdir=$2;
mkdir -p ${outputdir};

orig_fasta=${input_fasta}
deduplicate_src="/home/hieu/src/BCRTree_release/gctree/deduplicated.py";
modify_tree_colors="/home/hieu/src/BCRTree_release/gctree/modify_tree_colors.py";
color_path="/home/hieu/src/BCRTree_release/gctree/hex_color.csv";

filename=$(echo $input_fasta | xargs -n 1 basename);
filename=${filename%.aln.fasta*}
mkdir -p ${outputdir}/${filename}

cat ${input_fasta} | sed 's/>.*|Abundance:\([0-9]\+\)/>\1/' > ${outputdir}/${filename}/${filename}.modified.fa

input_fasta=${outputdir}/${filename}/${filename}.modified.fa;

echo -e "Running deduplicate ..."
##### We modify the python function deduplicate so that it allows repeated Sequence ID. 
# Sequence ID = sequence abundance values.
python ${deduplicate_src} \
--input ${input_fasta} \
--root GL \
--frame 0 \
--id_abundances \
--output_name ${filename} \
--output ${outputdir}/${filename}

echo -e "Running mkconfig..."
mkconfig --quick ${outputdir}/${filename}/${filename}.phylip dnapars \
> ${outputdir}/${filename}/${filename}_dnapars.cfg

echo -e "Running dnapars..."
dnapars < ${outputdir}/${filename}/${filename}_dnapars.cfg > ${outputdir}/${filename}/${filename}_dnapars.log

export QT_QPA_PLATFORM=offscreen
export XDG_RUNTIME_DIR=/tmp/runtime-runner
export MPLBACKEND=agg

echo -e "Inferring gctree ..."
gctree infer --verbose --root GL --frame 1 --idlabel outfile ${outputdir}/${filename}/${filename}.abundance.csv

mv outfile ${outputdir}/${filename}/outfile
mv outtree ${outputdir}/${filename}/outtree
mv gctree.out.* ${outputdir}/${filename}

echo -e "Modify tree colors..."
python ${modify_tree_colors} \
--input_fasta ${orig_fasta} \
--input_idmap ${outputdir}/${filename}/${filename}.id_map_seq.csv \
--gctree_inference_file ${outputdir}/${filename}/gctree.out.inference.1.p \
--color_path ${color_path} \
--output ${outputdir}/${filename} \
--svg_name ${filename}.color

echo -e "Finished"
