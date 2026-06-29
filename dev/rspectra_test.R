# %%
library(Matrix)
library(RSpectra)
library(dplyr)
library(tibble)


options(scipen = 9999)
# %%
# ==============================================================================
# 1. Construct a Large-Scale (500x500) Chromatin Domain Matrix
# ==============================================================================
set.seed(123)
n_bins <- 9000

# Step A: Initialize the 1D distance-decay background (Toeplitz structure)
# Mimics normal polymer diffusion where close bins interact more than distal ones
decay_rates <- 10 / (1:n_bins)^0.7
base_matrix <- toeplitz(decay_rates)

# Step B: Superimpose two massive, distinct structural domains (Sub-TADs)
# Domain 1: Bins 1 to 250 (High interaction block)
# Domain 2: Bins 251 to 500 (High interaction block)
base_matrix[1:250, 1:250]     <- base_matrix[1:250, 1:250] + matrix(runif(250*250, 2, 5), nrow = 250)
base_matrix[251:500, 251:500] <- base_matrix[251:500, 251:500] + matrix(runif(250*250, 2, 5), nrow = 250)

# Step C: Enforce perfect matrix symmetry and clean the diagonal
base_matrix[lower.tri(base_matrix)] <- t(base_matrix)[lower.tri(base_matrix)]
diag(base_matrix) <- 0

# Assign sequential genomic coordinate headers (e.g., 5kb bins starting at 10,000,000)
bin_coordinates <- seq(10000000, by = 5000, length.out = n_bins)
rownames(base_matrix) <- colnames(base_matrix) <- bin_coordinates

# Convert into the package's production-grade sparse matrix format
M_large_sparse <- Matrix::Matrix(base_matrix, sparse = TRUE)
# %%
# ==============================================================================
# 2. Benchmark Both Paths over the 500-Bin Space
# ==============================================================================

# Execution Path A: Direct dense calculation (Standard baseline)
time_eigen <- system.time({
  Dinv_e <- Matrix::Diagonal(nrow(M_large_sparse), 1 / Matrix::rowSums(M_large_sparse))
  lp_matrix_e <- Matrix::Diagonal(nrow(M_large_sparse), 1) - Dinv_e %*% M_large_sparse
  
  res_eigen <- eigen(lp_matrix_e)
  
  idx_fiedler <- n_bins - 1
  idx_zero    <- n_bins
  
  tbl_eigen <- tibble::as_tibble(res_eigen[["vectors"]][, c(idx_fiedler, idx_zero)], .name_repair = "minimal")
  colnames(tbl_eigen) <- c("fiedler", "zero")
  tbl_eigen <- tbl_eigen  |> # Standardize vector sign
    dplyr::mutate(bins = as.numeric(rownames(M_large_sparse)), path = "eigen")
})
# %%
# Execution Path B: Sparse shift-invert calculation (RSpectra scaling engine)
time_rspectra <- system.time({
  Dinv_r <- Matrix::Diagonal(nrow(M_large_sparse), 1 / Matrix::rowSums(M_large_sparse))
  lp_matrix_r <- Matrix::Diagonal(nrow(M_large_sparse), 1) - Dinv_r %*% M_large_sparse
  
  res_rspectra <- RSpectra::eigs(lp_matrix_r, k = 2, sigma = 1e-10, which = "LM", maxitr = 10000)
  
  tbl_rspectra <- tibble::as_tibble(res_rspectra[["vectors"]], .name_repair = "minimal")

  colnames(tbl_rspectra) <- c("fiedler", "zero")
  tbl_rspectra <- tbl_rspectra |>
    dplyr::mutate(bins = as.integer(rownames(M_large_sparse)), path = "rspectra")
})

plot(Re(tbl_rspectra[[1]]),tbl_eigen[[1]])


# %%

time_rspectra_sym <- system.time({
Dinv_sqrt_r <- Matrix::Diagonal(nrow(M_large_sparse), 1 / sqrt(Matrix::rowSums(M_large_sparse)))
# Construct L_sym (Stable, perfectly symmetric)
I_mat <- Matrix::Diagonal(n_bins, 1)
lp_sym <- I_mat - (Dinv_sqrt_r %*% M_large_sparse %*% Dinv_sqrt_r)
#temp <- RSpectra::eigs_sym(lp_sym, k = 2, sigma = 1e-10,which = "LM", maxitr = 10000)
temp <- RSpectra::eigs_sym(lp_sym, k = 2,which = "SM", maxitr = 10000)
# CRITICAL: Transform u_sym back to the recommended v_rw (Random-Walk) space
rw_vectors <- as.matrix(Dinv_sqrt_r %*% temp[["vectors"]])
tmp_tbl <- tibble::as_tibble(rw_vectors, .name_repair = "minimal")
colnames(tmp_tbl) <- c("fiedler", "zero")
tmp_tbl |> dplyr::mutate(bins = as.numeric(rownames(M_large_sparse)))
})
# %%


plot(tmp_tbl[[1]],tbl_rspectra[[1]])

# ==============================================================================
# 3. Print Performance and Alignment Verification
# ==============================================================================

message("\n=== PERFORMANCE PROFILING ===")
cat("Dense eigen() calculation time:    ", time_eigen["elapsed"], " seconds\n")
cat("Sparse RSpectra calculation time: ", time_rspectra_sym["elapsed"], " seconds\n")

# Verify that the calculated bisections split at the identical geometric midpoint
split_point_eigen    <- which(diff(sign(tbl_eigen$fiedler)) != 0)
split_point_rspectra <- which(diff(sign(Re(tmp_tbl$fiedler))) != 0)

message("\n=== GEOMETRIC SPLIT VERIFICATION ===")
cat("Domain cut index identified by eigen():    ", split_point_eigen, "\n")
cat("Domain cut index identified by RSpectra:  ", split_point_rspectra, "\n")
# %%
