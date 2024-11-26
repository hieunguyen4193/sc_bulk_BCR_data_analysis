library(dplyr)
data("combined_pbmc")

# piping the plot can be nice to read syntactically -
# By default, assigns unique colors to highlights and everything else is gray


# one useful application is to highlight shared clones - beware that the
# clonotype sequences may get extremely long in the legend
shared_aa_clones <- names(getSharedClones(combined_pbmc, clonecall = "aa"))

vizAPOTC(combined_pbmc, clonecall = "aa", verbose = FALSE) %>%
  showCloneHighlight(shared_aa_clones)
