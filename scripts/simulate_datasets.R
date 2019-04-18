##########################################
# Simulate data from from real data fits #
##########################################

library(ALStructure)
source("./misc_functions.R")
orig_fits_dir <- "/path/to/orig_fits/"
out_dir <- "/path/to/simulated/datasets/"
setwd(base_dir)

AI <- as.numeric(commandArgs()[length(commandArgs())])
print(AI)

# load in fits from original_data_fits.R (TGP)
load(paste0(orig_fits_dir, "PQ_ADX_tgp.Rdata"))
load(paste0(orig_fits_dir, "PQ_FS_tgp.Rdata"))
load(paste0(orig_fits_dir, "PQ_TS_tgp.Rdata"))
load(paste0(orig_fits_dir, "PQ_ALS_tgp.Rdata"))

# assign variables
P_adx_tgp <- ad_PQ$P
P_fs_tgp <- fs_PQ$P
P_ts_tgp <- ts_PQ$P
P_als_tgp <- als_PQ$P

Q_adx_tgp <- ad_PQ$Q
Q_fs_tgp <- fs_PQ$Q
Q_ts_tgp <- ts_PQ$Q
Q_als_tgp <- als_PQ$Q

m_tgp <- dim(P_adx_tgp)[1]
n_tgp <- dim(Q_adx_tgp)[2]
k_tgp <- 8;

# load in fits from poriginal_data_fits.R (HGDP)
load(paste0(orig_fits_dir, "PQ_ADX_hgdp.Rdata"))
load(paste0(orig_fits_dir, "PQ_FS_hgdp.Rdata"))
load(paste0(orig_fits_dir, "PQ_TS_hgdp.Rdata"))
load(paste0(orig_fits_dir, "PQ_ALS_hgdp.Rdata"))
P_adx_hgdp <- ad_PQ$P
P_fs_hgdp <- fs_PQ$P
P_ts_hgpd <- ts_PQ$P
P_als_hgdp <- als_PQ$P

Q_adx_hgdp <- ad_PQ$Q
Q_fs_hgdp <- fs_PQ$Q
Q_ts_hgdp <- ts_PQ$Q
Q_als_hgdp <- als_PQ$Q

m_hgdp <- dim(P_adx_hgdp)[1]
n_hgdp <- dim(Q_adx_hgdp)[2]
k_hgdp <- 10;

# load in fits from original_data_fits.R (HO)
load(paste0(orig_fits_dir, "PQ_ADX_ho.Rdata"))
load(paste0(orig_fits_dir, "PQ_FS_ho.Rdata"))
load(paste0(orig_fits_dir, "PQ_TS_ho.Rdata"))
load(paste0(orig_fits_dir, "PQ_ALS_ho.Rdata"))
P_adx_ho <- ad_PQ$P
P_fs_ho <- fs_PQ$P
P_ts_ho <- ts_PQ$P
P_als_ho <- als_PQ$P

Q_adx_ho <- ad_PQ$Q
Q_fs_ho <- fs_PQ$Q
Q_ts_ho <- ts_PQ$Q
Q_als_ho <- als_PQ$Q

m_ho <- dim(P_adx_ho)[1]
n_ho <- dim(Q_adx_ho)[2]
k_ho <- 10;

# Generate data (TGP)
F_adx_tgp <- P_adx_tgp %*% Q_adx_tgp;
F_fs_tgp <- P_fs_tgp %*% Q_fs_tgp;
F_ts_tgp <- P_ts_tgp %*% Q_ts_tgp;
F_als_tgp <- P_als_tgp %*% Q_als_tgp;
F_als_tgp[F_als_tgp > 1] <- 1 

set.seed(1234 * AI)
X_adx_tgp <- matrix(rbinom(n_tgp * m_tgp, 2, c(F_adx_tgp)), ncol = n_tgp, nrow = m_tgp)
X_fs_tgp <- matrix(rbinom(n_tgp * m_tgp, 2, c(F_fs_tgp)), ncol = n_tgp, nrow = m_tgp)
X_ts_tgp <- matrix(rbinom(n_tgp * m_tgp, 2, c(F_ts_tgp)), ncol = n_tgp, nrow = m_tgp)
X_als_tgp <- matrix(rbinom(n_tgp * m_tgp, 2, c(F_als_tgp)), ncol = n_tgp, nrow = m_tgp)

# Generate data (HGDP)
F_adx_hgdp <- P_adx_hgdp %*% Q_adx_hgdp;
F_fs_hgdp <- P_fs_hgdp %*% Q_fs_hgdp;
F_ts_hgdp <- P_ts_hgpd %*% Q_ts_hgdp;
F_als_hgdp <- P_als_hgdp %*% Q_als_hgdp;
F_als_hgdp[F_als_hgdp > 1] <- 1 

set.seed(1234 * AI)
X_adx_hgdp <- matrix(rbinom(n_hgdp * m_hgdp, 2, c(F_adx_hgdp)), ncol = n_hgdp, nrow = m_hgdp)
X_fs_hgdp <- matrix(rbinom(n_hgdp * m_hgdp, 2, c(F_fs_hgdp)), ncol = n_hgdp, nrow = m_hgdp)
X_ts_hgdp <- matrix(rbinom(n_hgdp * m_hgdp, 2, c(F_ts_hgdp)), ncol = n_hgdp, nrow = m_hgdp)
X_als_hgdp <- matrix(rbinom(n_hgdp * m_hgdp, 2, c(F_als_hgdp)), ncol = n_hgdp, nrow = m_hgdp)

# Generate data (HO)
F_adx_ho <- P_adx_ho %*% Q_adx_ho;
F_fs_ho <- P_fs_ho %*% Q_fs_ho;
F_ts_ho <- P_ts_ho %*% Q_ts_ho;
F_als_ho <- P_als_ho %*% Q_als_ho;
F_als_ho[F_als_ho > 1] <- 1 

set.seed(1234 * AI)
X_adx_ho <- matrix(rbinom(n_ho * m_ho, 2, c(F_adx_ho)), ncol = n_ho, nrow = m_ho)
X_fs_ho <- matrix(rbinom(n_ho * m_ho, 2, c(F_fs_ho)), ncol = n_ho, nrow = m_ho)
X_ts_ho <- matrix(rbinom(n_ho * m_ho, 2, c(F_ts_ho)), ncol = n_ho, nrow = m_ho)
X_als_ho <- matrix(rbinom(n_ho * m_ho, 2, c(F_als_ho)), ncol = n_ho, nrow = m_ho)

# save .Rdata files
save(X_adx_tgp, file = paste0(out_dir, "X_adx_tgp_", AI, ".Rdata"))
save(X_fs_tgp, file = paste0(out_dir, "X_fs_tgp_", AI, ".Rdata"))
save(X_ts_tgp, file = paste0(out_dir, "X_ts_tgp_", AI, ".Rdata"))
save(X_als_tgp, file = paste0(out_dir, "X_als_tgp_", AI, ".Rdata"))

save(X_adx_hgdp, file = paste0(out_dir, "X_adx_hgdp_", AI, ".Rdata"))
save(X_fs_hgdp, file = paste0(out_dir, "X_fs_hgdp_", AI, ".Rdata"))
save(X_ts_hgdp, file = paste0(out_dir, "X_ts_hgdp_", AI, ".Rdata"))
save(X_als_hgdp, file = paste0(out_dir, "X_als_hgdp_", AI, ".Rdata"))

save(X_adx_ho, file = paste0(out_dir, "X_adx_ho_", AI, ".Rdata"))
save(X_fs_ho, file = paste0(out_dir, "X_fs_ho_", AI, ".Rdata"))
save(X_ts_ho, file = paste0(out_dir, "X_ts_ho_", AI, ".Rdata"))
save(X_als_ho, file = paste0(out_dir, "X_als_ho_", AI, ".Rdata"))

# save .bed, .bim, .fam files
make_bed(X = X_adx_tgp, out_file = paste0(out_dir, "X_adx_tgp_", AI))
make_bed(X = X_fs_tgp, out_file = paste0(out_dir, "X_fs_tgp_", AI))
make_bed(X = X_ts_tgp, out_file = paste0(out_dir, "X_ts_tgp_", AI))
make_bed(X = X_als_tgp, out_file = paste0(out_dir, "X_als_tgp_", AI))

make_bed(X = X_adx_hgdp, out_file = paste0(out_dir, "X_adx_hgdp_", AI))
make_bed(X = X_fs_hgdp, out_file = paste0(out_dir, "X_fs_hgdp_", AI))
make_bed(X = X_ts_hgdp, out_file = paste0(out_dir, "X_ts_hgdp_", AI))
make_bed(X = X_als_hgdp, out_file = paste0(out_dir, "X_als_hgdp_", AI))

make_bed(X = X_adx_ho, out_file = paste0(out_dir, "X_adx_ho_", AI))
make_bed(X = X_fs_ho, out_file = paste0(out_dir, "X_fs_ho_", AI))
make_bed(X = X_ts_ho, out_file = paste0(out_dir, "X_ts_ho_", AI))
make_bed(X = X_als_ho, out_file = paste0(out_dir, "X_als_ho_", AI))
