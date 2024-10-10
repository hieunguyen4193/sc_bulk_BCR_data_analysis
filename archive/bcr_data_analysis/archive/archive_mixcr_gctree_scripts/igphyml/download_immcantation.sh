##### donwload the immcantation database
# Enter these commands in a terminal, not an R session!

# Move to the directory of interest
mkdir germlines

# Download the Immcantation repository
git clone https://bitbucket.org/kleinstein/immcantation

# Run script to obtain IMGT gapped sequences
immcantation/scripts/fetch_imgtdb.sh -o germlines

# View added directories
ls germlines
# human  IMGT.yaml  immcantation  mouse  rabbit  rat  rhesus_monkey