######################################################
# Script for cleaning dataset from Basu, et al. 2016 #
######################################################

# load required packages
library(lfa)
library(alstructure)

# set file path to location of India367 dataset files
data_dir <- "/path/to/dir/India367"

X_basu <- lfa::read.bed(data_dir)

# retrieve and parse .fam info
fam_df <- read_tsv(paste0(data_dir, ".fam"), col_names = FALSE)
fam_split <- strsplit(fam_df[[2]], "-|_")
fam <- matrix(0, nrow = n, ncol = 1)
for(i in 1:n){
  fam[i] <- fam_split[[i]][1]
}

# remove individuals from island populations
island_inds <- which(fam %in% c("JW", "ONG"))
X_basu <- X_basu[, -island_inds]

# write .bed, .bim, .fam files of cleaned dataset
alstructure:::make_bed(X_basu, paste0(data_dir, "_clean"))