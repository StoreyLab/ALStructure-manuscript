##############################################
# Fit datasets simulated from real data fits #
##############################################

library(alstructure)
library(lfa)
library(clue)
source("./misc_functions.R")

out_dir <- "/path/to/simulated_fits/"
true_dir <- "/path/to/original_fits/"
data_dir <- "/path/to/simulated/datasets/"
write_file <- paste0(out_dir, "all_vals.txt")

######################
# Method Directories #
######################
adx_dir <- "/path/to/admixture/software/admixture"
fs_dir <- "/path/to/faststructure/software/structure.py"
ts_dir <- "/path/to/terastructure/software/terastructure"
  
# k, m, and n values
ks <- c(8, 10, 14)
ms <- c(1229310, 550303, 372446)
ns <- c(1815, 940, 2251)
  
# get array index
AI <- commandArgs()[length(commandArgs())]
print(AI)
options(scipen = 10) # scientific output (important)

# set the current values of method and dataset for the current array index AI
real_datasets <- c("tgp", "hgdp", "ho")
fit1_methods <- c("adx", "fs", "ts", "als")
fit2_methods <- c("adx", "fs", "ts","als")
repetition <- c("1", "2", "3", "4")
parameters <- expand.grid(real_datasets, fit1_methods, fit2_methods, repetition)

# set simulation vals
dataset <- as.character(parameters[AI,1])
fit1_method <- as.character(parameters[AI,2])
fit2_method <- as.character(parameters[AI,3])
repetition <- as.character(parameters[AI, 4])

if(dataset == "tgp"){
  k <- ks[1]
  m <- ms[1]
  n <- ns[1]
} else if(dataset == "hgdp"){
  k <- ks[2]
  m <- ms[2]
  n <- ns[2]
} else if(dataset == "ho"){
  k <- ks[3]
  m <- ms[3]
  n <- ns[3]
}

# a list of the array elements to be recorded
column_names_master <- c("AI", "m", "n", "k", "fit1_method", "fit2_method", "dataset",
  "Q_ABS", "Q_RMSE", "P_ABS", "P_RMSE", "binomial_likelihood1", "binomial_llikelihood1", "binomial_likelihood2", "binomial_llikelihood2", "run_time")

# the name of the input dataset created by simulate_datasets.R
temp_data_str <- paste0("X_", fit1_method, "_", dataset, "_", repetition)

# read in data (for determining binomial likelihood later)
X <- lfa::read.bed(paste0(data_dir, temp_data_str))

if(fit2_method == "adx"){
  
  ###############
  # execute ADX #
  ###############
  
  # load "true" P and Q matrices
  if(fit1_method == "adx"){
    load(paste0(true_dir, "PQ_ADX_", dataset, ".Rdata"))
    true_PQ <- ad_PQ
  } else if(fit1_method == "fs"){
    load(paste0(true_dir, "PQ_FS_", dataset, ".Rdata"))
    true_PQ <- fs_PQ
  } else if(fit1_method == "ts"){
    load(paste0(true_dir, "PQ_TS_", dataset, ".Rdata"))
    true_PQ <- ts_PQ
  } else if(fit1_method == "als"){
    load(paste0(true_dir, "PQ_ALS_", dataset, ".Rdata"))
    true_PQ <- als_PQ
  }
  
  # output and execution strings
  fit_str <- paste0(out_dir, "ADX_", fit1_method, "_", fit2_method, "_", dataset, "_", repetition)
  str_ad_out <- paste0(fit_str, ".out")
  str_ad_exec <- paste0(adx_dir, " ", data_dir, temp_data_str, ".bed ", k, " > ", str_ad_out)
  
  # execute ADX
  print("ADX execution")
  system(str_ad_exec)
  
  # move the output files into place
  print(paste0("mv ", "./", temp_data_str, ".", k, ".P ", fit_str, ".P"))
  print(paste0("mv ", "./", temp_data_str, ".", k, ".Q ", fit_str, ".Q"))
  system(paste0("mv ", "./", temp_data_str, ".", k, ".P ", fit_str, ".P"))
  system(paste0("mv ", "./", temp_data_str, ".", k, ".Q ", fit_str, ".Q"))
  
  # record time
  time <- scrub_adx_time(str_ad_out)
  
  # obtain P and Q
  ad_PQ <- scrub_adx_PQ(fit_str)
  
  # find best rotation
  ad_perm <- naive_perm(true_PQ$Q, ad_PQ$Q)
  ad_PQ$Q <- ad_perm %*% ad_PQ$Q; ad_PQ$P <- ad_PQ$P %*% t(ad_perm)
  
  # compute errors
  Q_ABS_vals <- alstructure:::mean_ABS(ad_PQ$Q, true_PQ$Q)
  Q_RMSE_vals <- alstructure:::RMSE(ad_PQ$Q, true_PQ$Q)
  P_ABS_vals <- alstructure:::mean_ABS(ad_PQ$P, true_PQ$P)
  P_RMSE_vals <- alstructure:::RMSE(ad_PQ$P, true_PQ$P)
  liks1 <- alstructure:::binomial_likelihood(X, ad_PQ$P %*% ad_PQ$Q)
  liks2 <- alstructure:::binomial_likelihood(2 - X, ad_PQ$P %*% ad_PQ$Q)
  LL1 <- liks1$L
  ll1 <- liks1$l
  LL2 <- liks2$L
  ll2 <- liks2$l
  
  # write_info
  write(c(AI, m, n, k, fit1_method, fit2_method, dataset, Q_ABS_vals, Q_RMSE_vals, P_ABS_vals, P_RMSE_vals, LL1, ll1, LL2, ll2, time), file = write_file, ncolumns = length(column_names_master), append = TRUE, sep = "\t")
  
  ######################################
  # Save the fits for later simulation #
  ######################################
  save(ad_PQ, file = paste0(out_dir, "PQ_", fit1_method, "_", fit2_method, "_", dataset, "_", repetition, ".Rdata"))
}

if (fit2_method == "fs"){
  ##############
  # execute FS #
  ##############
  
  # load "true" P and Q matrices
  if(fit1_method == "adx"){
    load(paste0(true_dir, "PQ_ADX_", dataset, ".Rdata"))
    true_PQ <- ad_PQ
  } else if(fit1_method == "fs"){
    load(paste0(true_dir, "PQ_FS_", dataset, ".Rdata"))
    true_PQ <- fs_PQ
  } else if(fit1_method == "ts"){
    load(paste0(true_dir, "PQ_TS_", dataset, ".Rdata"))
    true_PQ <- ts_PQ
  } else if(fit1_method == "als"){
    load(paste0(true_dir, "PQ_ALS_", dataset, ".Rdata"))
    true_PQ <- als_PQ
  }
  
  # output and execution strings
  str_fs_out <- paste0(out_dir, "FS_", fit1_method, "_", fit2_method, "_", dataset, "_", repetition)
  str_fs_exec <- paste0("python2.7 ", fs_dir, " -K ", k, " --input=", data_dir, temp_data_str, " --out=", str_fs_out)
  print(str_fs_exec)
  
  # execute FS
  system(str_fs_exec)
  
  # record time
  time <- scrub_fs_time(paste0(str_fs_out, ".", k, ".log"))
  
  # obtain P and Q
  str_fs_PQ <- paste0(str_fs_out, ".", k)
  fs_PQ <- scrub_fs_PQ(str_fs_PQ)
  
  # find best rotation
  fs_perm <- naive_perm(true_PQ$Q, fs_PQ$Q)
  fs_PQ$Q <- fs_perm %*% fs_PQ$Q; fs_PQ$P <- fs_PQ$P %*% t(fs_perm)
  
  # compute errors
  Q_ABS_vals <- alstructure:::mean_ABS(fs_PQ$Q, true_PQ$Q)
  Q_RMSE_vals <- alstructure:::RMSE(fs_PQ$Q, true_PQ$Q)
  P_ABS_vals <- alstructure:::mean_ABS(fs_PQ$P, true_PQ$P)
  P_RMSE_vals <- alstructure:::RMSE(fs_PQ$P, true_PQ$P)
  liks1 <- alstructure:::binomial_likelihood(X, fs_PQ$P %*% fs_PQ$Q)
  liks2 <- alstructure:::binomial_likelihood(2 - X, fs_PQ$P %*% fs_PQ$Q)
  LL1 <- liks1$L
  ll1 <- liks1$l
  LL2 <- liks2$L
  ll2 <- liks2$l
  
  # write_info
  write(c(AI, m, n, k, fit1_method, fit2_method, dataset, Q_ABS_vals, Q_RMSE_vals, P_ABS_vals, P_RMSE_vals, LL1, ll1, LL2, ll2, time), file = write_file, ncolumns = length(column_names_master), append = TRUE, sep = "\t")
  
  ######################################
  # Save the fits for later simulation #
  ######################################
  save(fs_PQ, file = paste0(out_dir, "PQ_", fit1_method, "_", fit2_method, "_", dataset, "_", repetition, ".Rdata"))
}

if(fit2_method == "ts"){
  ##############
  # execute TS #
  ##############
  
  # load "true" P and Q matrices
  if(fit1_method == "adx"){
    load(paste0(true_dir, "PQ_ADX_", dataset, ".Rdata"))
    true_PQ <- ad_PQ
  } else if(fit1_method == "fs"){
    load(paste0(true_dir, "PQ_FS_", dataset, ".Rdata"))
    true_PQ <- fs_PQ
  } else if(fit1_method == "ts"){
    load(paste0(true_dir, "PQ_TS_", dataset, ".Rdata"))
    true_PQ <- ts_PQ
  } else if(fit1_method == "als"){
    load(paste0(true_dir, "PQ_ALS_", dataset, ".Rdata"))
    true_PQ <- als_PQ
  }
  
  # write the locations.loc file
  locations <- rep(c(0:(m - 1)), each = 2)
  locations_dir <- paste0(data_dir, "locations_", fit1_method, "_", fit2_method, "_", dataset, "_", repetition, ".loc")
  locations_mat <- matrix(locations, nrow = m, ncol = 2, byrow = TRUE)
  write.table(locations_mat, file = locations_dir, sep = "\t", eol = "\n", row.names = FALSE, col.names = FALSE)
  
  # output and execution strings
  str_ts_dir <- paste0("./n", n, "-k", k, "-l", m, "-", fit1_method, "-", fit2_method, "-", dataset, "-", repetition)
  str_ts_exec <- paste0(ts_dir, " -file ", data_dir, temp_data_str, ".bed -n ", n, " -l ", m, " -k ", k, " -stochastic -rfreq 10000 -force -label ", fit1_method, "-", fit2_method, "-", dataset, "-", repetition)
  str_ts_beta_run <- paste0(ts_dir, " -file ", data_dir, temp_data_str, ".bed -locations-file ", locations_dir, " -n ", n, " -l ", m ," -k ", k, " -stochastic -label ", fit1_method, "-", fit2_method, "-", dataset, "-", repetition, "-beta -compute-beta")
  
  # execute TS
  print("TS execution")
  system(str_ts_exec)
  p_time <- proc.time()
  print(str_ts_dir)
  setwd(str_ts_dir)
  system(str_ts_beta_run)
  beta_time <- proc.time() - p_time
  
  # record time
  time <- scrub_ts_time(paste0(str_ts_dir, "/infer.log"), as.numeric(beta_time[3]))
  
  # obtain P and Q
  str_ts_Q <- str_ts_dir
  str_ts_P <- paste0(str_ts_dir, "/n", n, "-k", k, "-l", m, "-", fit1_method, "-", fit2_method, "-", dataset, "-", repetition, "-beta")
  ts_PQ <- scrub_ts_PQ(str_ts_P, str_ts_Q)
  
  # find best rotation
  ts_perm <- naive_perm(true_PQ$Q, ts_PQ$Q)
  ts_PQ$Q <- ts_perm %*% ts_PQ$Q; ts_PQ$P <- ts_PQ$P %*% t(ts_perm)
  
  # compute errors
  Q_ABS_vals <- alstructure:::mean_ABS(ts_PQ$Q, true_PQ$Q)
  Q_RMSE_vals <- alstructure:::RMSE(ts_PQ$Q, true_PQ$Q)
  P_ABS_vals <- alstructure:::mean_ABS(ts_PQ$P, true_PQ$P)
  P_RMSE_vals <- alstructure:::RMSE(ts_PQ$P, true_PQ$P)
  liks1 <- alstructure:::binomial_likelihood(X, ts_PQ$P %*% ts_PQ$Q)
  liks2 <- alstructure:::binomial_likelihood(2 - X, ts_PQ$P %*% ts_PQ$Q)
  LL1 <- liks1$L
  ll1 <- liks1$l
  LL2 <- liks2$L
  ll2 <- liks2$l
  
  write(c(AI, m, n, k, fit1_method, fit2_method, dataset, Q_ABS_vals, Q_RMSE_vals, P_ABS_vals, P_RMSE_vals, LL1, ll1, LL2, ll2, time), file = write_file, ncolumns = length(column_names_master), append = TRUE, sep = "\t")
  
  ######################################
  # Save the fits for later simulation #
  ######################################
  save(ts_PQ, file = paste0(out_dir, "PQ_", fit1_method, "_", fit2_method, "_", dataset, "_", repetition, ".Rdata"))
  
  # move the files into place
  system(paste0("mv ", str_ts_dir, " ", out_dir))
}

if (fit2_method == "als"){
  ###############
  # execute ALS #
  ###############
  
  # load "true" P and Q matrices
  if(fit1_method == "adx"){
    load(paste0(true_dir, "PQ_ADX_", dataset, ".Rdata"))
    true_PQ <- ad_PQ
  } else if(fit1_method == "fs"){
    load(paste0(true_dir, "PQ_FS_", dataset, ".Rdata"))
    true_PQ <- fs_PQ
  } else if(fit1_method == "ts"){
    load(paste0(true_dir, "PQ_TS_", dataset, ".Rdata"))
    true_PQ <- ts_PQ
  } else if(fit1_method == "als"){
    load(paste0(true_dir, "PQ_ALS_", dataset, ".Rdata"))
    true_PQ <- als_PQ
  }
  
  # output string
  str_als_dir <- paste0(out_dir, fit1_method, "_", fit2_method, "_", dataset, "_", repetition)
  system(paste0("mkdir ", str_als_dir))
  
  # execute ALS method
  print("ALS execution")
  time_init <- proc.time()
  als_PQ <- alstructure(X, d_hat = k, svd_method = "base", tol = 0.00001, max_iters = 10000)
  time_fin <- proc.time()
  time <- (time_fin - time_init)[3]
  
  # find best rotation
  als_perm <- naive_perm(true_PQ$Q, als_PQ$Q)
  als_PQ$Q <- als_perm %*% als_PQ$Q; als_PQ$P <- als_PQ$P %*% t(als_perm)
  
  # compute errors
  Q_ABS_vals <- alstructure:::mean_ABS(als_PQ$Q, true_PQ$Q)
  Q_RMSE_vals <- alstructure:::RMSE(als_PQ$Q, true_PQ$Q)
  P_ABS_vals <- alstructure:::mean_ABS(als_PQ$P, true_PQ$P)
  P_RMSE_vals <- alstructure:::RMSE(als_PQ$P, true_PQ$P)
  liks1 <- alstructure:::binomial_likelihood(X, als_PQ$P %*% als_PQ$Q)
  liks2 <- alstructure:::binomial_likelihood(2 - X, als_PQ$P %*% als_PQ$Q)
  LL1 <- liks1$L
  ll1 <- liks1$l
  LL2 <- liks2$L
  ll2 <- liks2$l
  
  # write_info
  write(c(AI, m, n, k, fit1_method, fit2_method, dataset, Q_ABS_vals, Q_RMSE_vals, P_ABS_vals, P_RMSE_vals, LL1, ll1, LL2, ll2, time), file = write_file, ncolumns = length(column_names_master), append = TRUE, sep = "\t")
  
  # write P and Q
  write.table(als_PQ$P, file = paste0(str_als_dir, "/P.txt"), col.names = FALSE, row.names = FALSE, sep = "\t")
  write.table(als_PQ$Q, file = paste0(str_als_dir, "/Q.txt"), col.names = FALSE, row.names = FALSE, sep = "\t")
  
  ######################################
  # Save the fits for later simulation #
  ######################################
  save(als_PQ, file = paste0(out_dir, "PQ_", fit1_method, "_", fit2_method, "_", dataset, "_", repetition, ".Rdata"))
  
}
