echo -e "#############################################";
echo -e "Working on sample " ${sampleid} "\n";
echo -e "#############################################";
sampleid=$1;
fastq1=$2;
fastq2=$3;

inputdir="/home/hieu/src/BCRTree/input";
output="/home/hieu/outdir/mixcr_pipeline_output/mid_based_output/"${sampleid};
mkdir -p ${output};

read_pattern="^N{35}(UMI:N{16})N{4:6}(R1:*)\^N{34}(R2:*)"

#####--------------------------------------------------#####
##### UPSTREAM PIPELINE
#####--------------------------------------------------#####

#### ALIGNMENT
echo -e "ALIGNMENT"

mixcr align \
-p generic-amplicon-with-umi \
--rna \
--library imgt.202312-3 \
--species mmu \
--report ${output}/${sampleid}.report.txt \
--json-report ${output}/${sampleid}.align.report.json \
--tag-pattern ${read_pattern} \
--floating-left-alignment-boundary \
--floating-right-alignment-boundary C \
$fastq1 \
$fastq2 \
${output}/${sampleid}.vdjca

##### Refine tags and sort
echo -e "Refine tags and sort"

mixcr refineTagsAndSort \
    --report ${output}/${sampleid}.refine.report.txt \
    --json-report ${output}/${sampleid}.refine.report.json \
    ${output}/${sampleid}.vdjca \
    ${output}/${sampleid}.refined.vdjca

##### Assemble
echo -e "Assemble"

mixcr assemble \
    --report ${output}/${sampleid}.assemble.report.txt \
    --json-report ${output}/${sampleid}.assemble.report.json \
    -OassemblingFeatures={FR1Begin:FR4End} \
    -OseparateByC=true \
    ${output}/${sampleid}.refined.vdjca \
    ${output}/${sampleid}.clns

##### Export CLONES
echo -e "EXRPOT CLONES"
mixcr exportClones -uniqueTagCount Molecule -count \
${output}/${sampleid}.clns ${output}/${sampleid}.tsv

##### Export QUALITY CONTROLS
echo -e "Export QC"
mixcr exportQc align ${output}/*.clns ${output}/alignQc.pdf

#####--------------------------------------------------#####
##### TREES
#####--------------------------------------------------#####
mixcr findAlleles \
    --report ${output}/${sampleid}.findAlleles.report.txt \
    --json-report ${output}/${sampleid}.findAlleles.report.json \
    --export-alleles-mutations ${output}/${sampleid}_alleles.tsv \
    --export-library ${output}/${sampleid}_alleles.json \
    --output-template {file_dir_path}/{file_name}.reassigned.clns \
    ${output}/${sampleid}.clns

mixcr exportClones -uniqueTagCount Molecule -count ${output}/${sampleid}.reassigned.clns ${output}/${sampleid}.reassigned.tsv
