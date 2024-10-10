library(dowser)
library(dplyr)

data(ExampleAirr)

# Read in IMGT-gapped sequences
path.to.project.src <- "/home/hieunguyen/CRC1382/src_2023/bcrtree/igphyml"
references = readIMGT(dir = file.path(path.to.project.src, "germlines", "mouse", "vdj"))

# remove germline alignment columns for this example
db = select(ExampleAirr, -"germline_alignment", 
            -"germline_alignment_d_mask")

# Reconstruct germline sequences
ExampleAirr = createGermlines(db, references, nproc=1)

# Check germline of first row
ExampleAirr$germline_alignment_d_mask[1]