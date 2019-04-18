################################
# run all methods on real data #
################################

library(alstructure)
library(lfa)
source("./misc_functions.R")

####################
# Data Directories #
####################
base_dir <- "./data/"
out_dir <- "/path/to/original_fits/"
write_file <- paste0(out_dir, "all_vals.txt")

######################
# Method Directories #
######################
adx_dir <- "/path/to/admixture/software/admixture"
fs_dir <- "/path/to/faststructure/software/structure.py"
ts_dir <- "/path/to/terastructure/software/terastructure"

# get array index
AI <- commandArgs()[length(commandArgs())]
print(AI)

options(scipen = 10) # scientific output (important)

# set the current values of method and dataset for the current array index AI
real_datasets <- c("tgp", "hgdp", "ho", "basu")
methods <- c("ADX", "FS", "TS", "ALS")
parameters <- expand.grid(real_datasets, methods)
dataset <- as.character(parameters[AI,1])
method <- as.character(parameters[AI,2])
  
# parameters
ks <- c(8, 10, 14, 4)
ms <- c(1229310, 550303, 372446, 803570)
ns <- c(1815, 940, 2251, 331)
  
setwd(base_dir)
  
# real if statemets
if(dataset == "tgp"){
  k <- ks[1]
  m <- ms[1]
  n <- ns[1]
  temp_data_str <- paste0("TGP_phase3_chip")
} else if(dataset == "hgdp"){
  k <- ks[2]
  m <- ms[2]
  n <- ns[2]
  temp_data_str <- paste0("HGDP_hwe_940")
} else if(dataset == "ho"){
  k <- ks[3]
  m <- ms[3]
  n <- ns[3]
  temp_data_str <- paste0("Lazaridis")
} else if(dataset == "basu"){
  k <- ks[4]
  m <- ms[4]
  n <- ns[4]
  temp_data_str <- paste0("India367_clean")
}

# a list of the array elements to be recorded
column_names_master <- c("m", "n", "k", "dataset", "method", "run_time")

if(method == "ADX"){
  
  ###############
  # execute ADX #
  ###############
  
  # output and execution strings
  str_ad_out <- paste0(out_dir, "ADX_", dataset, ".out")
  str_ad_exec <- paste0(adx_dir, " ", base_dir, temp_data_str, ".bed ", k, " > ", str_ad_out)
  
  # execute ADX
  print("ADX execution")
  system(str_ad_exec)
  
  # record time
  time <- scrub_adx_time(str_ad_out)
  
  # obtain P and Q
  str_ad_PQ <- paste0(base_dir, temp_data_str, ".", k)
  ad_PQ <- scrub_adx_PQ(str_ad_PQ)
  
  # write_info
  write(c(m, n, k, dataset, "ADX", time), file = write_file, ncolumns = length(column_names_master), append = TRUE, sep = "\t")
  
  ######################################
  # Save the fits for later simulation #
  ######################################
  save(ad_PQ, file = paste0(out_dir, "PQ_", method, "_", dataset, ".Rdata"))
  
  # move the output files into place
  system(paste0("mv ", temp_data_str, ".", k, ".P ", out_dir))
  system(paste0("mv ", temp_data_str, ".", k, ".Q ", out_dir))
}

if (method == "FS"){
  ##############
  # execute FS #
  ##############
  
  # output and execution strings
  str_fs_out <- paste0(out_dir, "FS_", dataset)
  str_fs_exec <- paste0("python2.7 ", fs_dir, " -K ", k, " --input=", temp_data_str, " --out=", str_fs_out)
  
  # execute FS
  print("FS execution")
  system(str_fs_exec)
  
  # record time
  time <- scrub_fs_time(paste0(str_fs_out, ".", k, ".log"))
  
  # obtain P and Q
  str_fs_PQ <- paste0(str_fs_out, ".", k)
  fs_PQ <- scrub_fs_PQ(str_fs_PQ)
  
  # write_info
  write(c(m, n, k, dataset, "FS", time), file = write_file, ncolumns = length(column_names_master), append = TRUE, sep = "\t")
  
  ######################################
  # Save the fits for later simulation #
  ######################################
  save(fs_PQ, file = paste0(out_dir, "PQ_", method, "_", dataset, ".Rdata"))
}

if(method == "TS"){
  ##############
  # execute TS #
  ##############
  
  # write the locations.loc file
  locations <- rep(c(0:(m - 1)), each = 2)
  locations_mat <- matrix(locations, nrow = m, ncol = 2, byrow = TRUE)
  write.table(locations_mat, file = paste0("locations_", method, "_", dataset, ".loc"), sep = "\t", eol = "\n", row.names = FALSE, col.names = FALSE)
  
  # output and execution strings
  str_ts_dir <- paste0(base_dir, "n", n, "-k", k, "-l", m, "-TS-", dataset)
  str_ts_exec <- paste0(ts_dir, " -file ", temp_data_str, ".bed -n ", n, " -l ", m, " -k ", k, " -stochastic -rfreq 10000 -force -label ", "TS-", dataset)
  str_ts_beta_run <- paste0(ts_dir, " -file ../", temp_data_str, ".bed -locations-file ../locations_", method, "_", dataset, ".loc -n ", n, " -l ", m ," -k ", k, " -stochastic -label ", method, "-", dataset, "-beta -compute-beta")
  
  # execute TS
  print("TS execution")
  system(str_ts_exec)
  p_time <- proc.time()
  setwd(str_ts_dir)
  system(str_ts_beta_run)
  setwd(base_dir)
  beta_time <- proc.time() - p_time
  
  # record time
  time <- scrub_ts_time(paste0(str_ts_dir, "/infer.log"), as.numeric(beta_time[3]))
  
  # obtain P and Q
  str_ts_Q <- str_ts_dir
  str_ts_P <- paste0(str_ts_dir, "/n", n, "-k", k, "-l", m, "-", method, "-", dataset, "-beta")
  ts_PQ <- scrub_ts_PQ(str_ts_P, str_ts_Q)
  
  # write_info
  write(c(m, n, k, dataset, "TS", time), file = write_file, ncolumns = length(column_names_master), append = TRUE, sep = "\t")
  
  ######################################
  # Save the fits for later simulation #
  ######################################
  save(ts_PQ, file = paste0(out_dir, "PQ_", method, "_", dataset, ".Rdata"))
  
  # move the terastructure folder to the proper place
  system(paste0("mv ", str_ts_dir, " ", out_dir))
}

if (method == "ALS"){
  ###############
  # execute ALS #
  ###############
  
  bed_file <- paste0(base_dir, temp_data_str)
  X <- lfa::read.bed(bed_file)
  
  # output string and make directory
  str_als_dir <- paste0(out_dir, "ALS_", dataset)
  system(paste0("mkdir ", str_als_dir))
  
  # execute ALS method
  print("ALS execution")
  time_init <- proc.time()
  als_PQ <- alstructure(X, d_hat = k, svd_method = "base", tol = 0.0001, max_iters = 1000)
  time_final <- proc.time()
  time <- time_final - time_init
  time <- time[3]
  names(als_PQ)[1] <- "P"
  names(als_PQ)[2] <- "Q"
  
  # write_info
  write(c(m, n, k, dataset, "ALS", time), file = write_file, ncolumns = length(column_names_master), append = TRUE, sep = "\t")
  
  # write P and Q
  write.table(als_PQ$P, file = paste0(str_als_dir, "/P.txt"), col.names = FALSE, row.names = FALSE, sep = "\t")
  write.table(als_PQ$Q, file = paste0(str_als_dir, "/Q.txt"), col.names = FALSE, row.names = FALSE, sep = "\t")
  
  ######################################
  # Save the fits for later simulation #
  ######################################
  save(als_PQ, file = paste0(out_dir, "PQ_", method, "_", dataset, ".Rdata"))
  
}
