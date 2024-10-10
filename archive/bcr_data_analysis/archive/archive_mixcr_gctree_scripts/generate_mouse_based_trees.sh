inputdir="/home/hieu/outdir/mixcr_pipeline_output/mouse_based_output";
all_mouses=$(ls $inputdir)
for mouse in $all_mouses;do \
echo -e "working on mouse sample " $mouse;

echo -e "find alleles";
mixcr findAlleles \
      --report $inputdir/$mouse/${mouse}.findAlleles.report.txt \
      --export-alleles-mutations $inputdir/$mouse/${mouse}_alleles.tsv \
      --export-library $inputdir/$mouse/${mouse}_alleles.json \
      --output-template {file_dir_path}/{file_name}.reassigned.clns \
      $inputdir/$mouse/*.clns;

#echo -e "export clones";
#mixcr exportClones \ 
#      $inputdir/$mouse/${mouse}.reassigned.clns \
#      $inputdir/$mouse/${mouse}.reassigned_exportClones.tsv;

echo -e "find SHM trees";
mixcr findShmTrees \
      --report $inputdir/$mouse/${mouse}_trees.log \
      $inputdir/$mouse/*.reassigned.clns \
      $inputdir/$mouse/${mouse}.shmt;

echo -e "export SHM trees with nodes";
mixcr exportShmTreesWithNodes \
      $inputdir/$mouse/${mouse}.shmt \
      $inputdir/$mouse/${mouse}_trees.tsv;

echo -e "export SHM Trees NEWICK";
mixcr exportShmTreesNewick\
    $inputdir/$mouse/${mouse}.shmt \
    $inputdir/$mouse/IM_newick;
done
