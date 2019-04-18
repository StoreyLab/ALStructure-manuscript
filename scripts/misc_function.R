###########################################
# Utility functions used in data analyses #
###########################################

scrub_adx_time <- function(file){
  # function to scrub time of evaluation from ADX .out file
  
  adx_out <- readLines(file)
  time_index <- grep("^Converged", adx_out)
  time_line <- strsplit(adx_out[time_index], split = " ")
  time_raw <- strsplit(time_line[[1]][5], split = "")[[1]][-1]
  time_adx <- as.numeric(paste(time_raw, collapse = ""))
}

scrub_fs_time <- function(file){
  # function to scrub time of evaluation from ADX .out file
  
  fs_out <- readLines(file)
  time_index <- grep("^Total time", fs_out)
  time_line <- strsplit(fs_out[time_index], split = " ")
  time_fs <- as.numeric(time_line[[1]][4])
}

scrub_ts_time <- function(file, p_time){
  # function to scrub the time of evaluation from TS infer.log file
  
  ts_out <- readLines(file)
  time_line <- strsplit(ts_out[length(ts_out) - 1], split = " ")
  time_ts <- as.numeric(time_line[[1]][length(time_line[[1]]) - 1])
  time_ts <- time_ts + p_time
}


scrub_als_time <- function(file){
  # function to scrub the time of evaluation from ALStructure summary.txt
  
  als_out <- readLines(file)
  time_line_ind <- grep("^time:", als_out)
  time_line_split <- strsplit(als_out[time_line_ind], split = " ")
  time_als <- as.numeric(time_line_split[[1]][2])
}

scrub_adx_PQ <- function(file){
  # function to scrub the P and Q matrices from ADX output
  file_P <- paste0(file, ".P")
  file_Q <- paste0(file, ".Q")
  P <- as.matrix(read.table(file_P))
  Q <- t(as.matrix(read.table(file_Q)))
  rownames(P) <- NULL; rownames(Q) <- NULL
  
  final_result = list(P = P, Q = Q)
  return(final_result)
}

scrub_fs_PQ <- function(file){
  # function to scrub the P and Q matrices from fs output
  file_P <- paste0(file, ".meanP")
  file_Q <- paste0(file, ".meanQ")
  P <- as.matrix(read.table(file_P, row.names = NULL))
  Q <- t(as.matrix(read.table(file_Q, row.names = NULL)))
  rownames(P) <- NULL; rownames(Q) <- NULL
  
  final_result = list(P = P, Q = Q)
  return(final_result)
}

scrub_ts_PQ <- function(P_folder, Q_folder){
  # function to scrub the P and Q matrices from fs output
  file_P <- paste0(P_folder, "/beta.txt")
  file_Q <- paste0(Q_folder, "/theta.txt")
  P <- as.matrix(read.table(file_P, row.names = NULL))
  Q <- t(as.matrix(read.table(file_Q, row.names = NULL)))
  P <- P[, 2:dim(P)[2]]; class(P) <- "numeric"
  Q <- (Q[3:(dim(Q)[1] - 1), ]); class(Q) <- "numeric"
  rownames(P) <- NULL; rownames(Q) <- NULL
  
  final_result = list(P = P, Q = Q)
  return(final_result)
}

mash_txts <- function(){
  # function to mash together all .txt files in all_vals folder and combine them
  # into a single .txt master file. mash_txts() will also save all of this into
  # a data.frame
  
  setwd("./")
  str_vals_dir <- "./all_vals"
  file.names <- dir(str_vals_dir, pattern =".txt")
  
  # initialize the main master data frame
  df_master <- data.frame()
  
  # loop through file names and write 
  for (i in 1:length(file.names)){
    master_vec <- read.table(paste0("./all_vals/",file.names[i]), header = FALSE, sep = "\t")
    df_master <- rbind(df_master, master_vec)
  }
  
  all_names <- c("AI", "m", "n", "k", "alpha_proto", "alpha_rep", "method",
    "Q_KL", "Q_ABS", "Q_RMSE", "P_KL", "P_ABS", "P_RMSE", "run_time")
  names(df_master) <- all_names
  df_master <- df_master[order(df_master$AI),]
  save(df_master, file = "df_master.RData")
  write.table(df_master, file = "master_table", row.names = FALSE)
}

naive_perm <- function(A, B){
  # function finds best permutation matrix P such that PB is close to A.
  # expects a short, fat matrix. This is a greedy algorithm
  
  d <- dim(A)[1]
  norms <- matrix(0, nrow = d, ncol = d, dimnames = list(1:d, 1:d))
  perm_mat <- matrix(0, nrow = d, ncol = d)
  
  for(i in 1:d){
    for(j in 1:d){
      norms[i, j] <- sum((A[i, ] - B[j, ])^2)
    }
  }
  
  # solve the asymmetric travelling salesman problem)
  assignment <- as.integer(solve_LSAP(norms))
  for(i in 1:d){
    perm_mat[i, assignment[i]] <- 1
  }
  
  perm_mat
  
}