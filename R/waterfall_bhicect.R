# %%
library(igraph)
library(hictkR)
library(dplyr)
library(MASS)
# %%

#' Generate a Cleaned Adjacency Matrix from an igraph Object
#'
#' This function takes a chromosomal network (igraph object) and converts it into 
#' a weighted adjacency matrix. It ensures the matrix has no self-loops by 
#' zeroing the diagonal and removes any "orphan" nodes (nodes with no connections).
#'
#' @param g_chr1 An \code{igraph} object, typically representing a chromosomal 
#'   interaction network.
#'
#' @return A weighted \code{dgCMatrix} (sparse matrix format) with zero 
#'   diagonals and no zero-sum columns/rows.
#' 
#' @export
#'
#' @importFrom igraph as_adjacency_matrix
#' @importFrom Matrix colSums
#'
#' @examples
#' \dontrun{
#' # Assuming g is a pre-existing igraph object
#' adj_matrix <- get_adj_mat_fn(g)
#' }
get_adj_mat_fn <- function(g_chr1) {
# Convert to weighted adjacency matrix
  chr_mat <- igraph::as_adjacency_matrix(g_chr1, type = "both", attr = "weight")
# Ensure no self-loops (diagonal = 0) 
  diag(chr_mat) <- 0
# Identify and remove orphan nodes (zero-sum columns/rows)
  if (any(Matrix::colSums(chr_mat) == 0)) {
    out <- which(Matrix::colSums(chr_mat) == 0)
    chr_mat <- chr_mat[-out, ]
    chr_mat <- chr_mat[, -out]
  }
  return(chr_mat)
}

#' Compute the Normalized Laplacian and Fiedler Vector
#'
#' This function calculates the normalized Graph Laplacian of an adjacency matrix 
#' and extracts the Fiedler vector (the eigenvector corresponding to the second 
#' smallest eigenvalue). It switches between base \code{eigen} and 
#' \code{RSpectra::eigs_sym} based on matrix size for computational efficiency.
#'
#' @param x A square \code{dgCMatrix} or \code{matrix} representing a 
#'   weighted adjacency matrix.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{vectors}: A \code{tibble} with the Fiedler vector, the zero 
#'   eigenvector, and the original row indices (bins).
#'   \item \code{values}: The corresponding eigenvalues.
#' }
#' 
#' @details 
#' The normalized Laplacian is calculated as:
#' \eqn{L = I - D^{-1}A}
#' where \eqn{I} is the identity matrix, \eqn{D} is the degree matrix, 
#' and \eqn{A} is the adjacency matrix.
#'
#' @importFrom Matrix Diagonal rowSums
#' @importFrom tibble as_tibble
#' @importFrom dplyr mutate
#' @importFrom RSpectra eigs_sym
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming adj is a sparse adjacency matrix from get_adj_mat_fn
#' result <- lp_fn(adj)
#' plot(result$vectors$fiedler)
#' }
lp_fn <- function(x) {
  Dinv <- Matrix::Diagonal(nrow(x), 1 / Matrix::rowSums(x))

  lp_chr1 <- Matrix::Diagonal(nrow(x), 1) - Dinv %*% x
  if (dim(lp_chr1)[1] > 10000) {
    temp <- RSpectra::eigs_sym(lp_chr1, k = 2, sigma = .Machine$double.xmin, which = "LM", maxitr = 10000)
    tmp_tbl <- tibble::as_tibble(temp[["vectors"]],.name_repair = make.names)
    colnames(tmp_tbl) <- c("fiedler", "zero")
    tmp_tbl <- tmp_tbl |>
      dplyr::mutate(bins = as.integer(rownames(x)))

    return(list(vectors = tmp_tbl, values = temp[["values"]]))
  } else {
    temp <- eigen(lp_chr1)
    tmp_tbl <- tibble::as_tibble(temp[["vectors"]][, c(length(temp$values) - 1, length(temp$values))],.name_repair = make.names)
    colnames(tmp_tbl) <- c("fiedler", "zero")
    tmp_tbl <- tmp_tbl %>%
      dplyr::mutate(bins = as.integer(rownames(x)))
    return(list(vectors = tmp_tbl, values = temp[["values"]][c(length(temp$values) - 1, length(temp$values))]))
  }
}

#' Calculate the Total Sum of Squares (SS)
#'
#' This is a helper function that calculates the sum of squared deviations 
#' from the mean. It is equivalent to calculating the squared Euclidean 
#' norm of a centered vector.
#'
#' @param x A numeric vector or matrix.
#'
#' @return A numeric value representing the sum of squares.
#' 
#' @details 
#' The function uses \code{scale(x, scale = FALSE)} to center the data 
#' around the mean before squaring, which is numerically more stable 
#' than the raw formula \eqn{\sum x^2 - (\sum x)^2/n}.
#'
#' @export
#'
#' @examples
#' data <- c(1, 2, 3, 4, 5)
#' ss(data)
#' 
#' # Example in a spatial context: variance of a gene expression vector
#' gene_expr <- c(0.1, 1.2, 0.5, 2.1)
#' total_ss <- ss(gene_expr)
ss <- function(x) {

  sum(scale(x, scale = FALSE)^2)
}

#' Perform Optimal Bi-partitioning on a Fiedler Vector
#'
#' This function identifies the optimal threshold for splitting a network into 
#' two clusters based on the Fiedler vector. It iterates through all possible 
#' split points and selects the one that maximizes the ratio of the 
#' between-cluster sum of squares to the total sum of squares.
#'
#' @param lp_res A list containing the results from \code{lp_fn}, specifically 
#'   a \code{vectors} tibble containing the "fiedler" column.
#' @param tmp_res A numeric or character value indicating the current 
#'   resolution or iteration identifier to be appended to the output.
#'
#' @return A \code{tibble} containing:
#' \itemize{
#'   \item \code{real_fiedler}: The real component of the Fiedler vector.
#'   \item \code{stat}: The variance-explained ratio for that specific split point.
#'   \item \code{smpl.cl}: The resulting cluster assignment (1 or 2).
#'   \item \code{res}: The provided \code{tmp_res} identifier.
#' }
#' 
#' @details 
#' The function calculates a "stat" for every value in the Fiedler vector. 
#' This stat represents the Sum of Squares of the predicted means divided by 
#' the total Sum of Squares. The split that maximizes this value is chosen 
#' as the partition point.
#'
#' @importFrom dplyr mutate slice_max pull
#' @importFrom purrr map_dbl
#' @importFrom stats Re
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming lp_results is the output from lp_fn()
#' partition <- simple_partition_tbl_fn(lp_results, tmp_res = 1)
#' }
simple_partition_tbl_fn <- function(lp_res, tmp_res) {
  smpl_thresh_tbl <- lp_res$vectors |>
    dplyr::mutate(real_fiedler = Re(fiedler))|>
    dplyr::mutate(stat = purrr::map_dbl(real_fiedler, function(x) {
      cl_a <- real_fiedler[which(real_fiedler <= x)]
      cl_b <- real_fiedler[which(real_fiedler > x)]
      return(ss(rep(c(mean(cl_a), mean(cl_b)), c(length(cl_a), length(cl_b)))) / ss(real_fiedler))
    }))
  smpl_thresh <- smpl_thresh_tbl |>
    dplyr::slice_max(stat,n = 1, with_ties = FALSE) |>
    dplyr::pull(real_fiedler)
  smpl_thresh_tbl <- smpl_thresh_tbl |>
    dplyr::mutate(
      smpl.cl = ifelse(real_fiedler <= smpl_thresh, 1, 2),
      res = tmp_res
    )
  smpl_thresh_tbl
}

# %%

#' Execute Spectral Bi-partitioning on a Data Frame
#'
#' This high-level wrapper orchestrates the full spectral clustering pipeline: 
#' it constructs a graph from a data frame, cleans self-loops, generates an 
#' adjacency matrix, computes the Laplacian, and identifies the optimal 
#' bi-partition.
#'
#' @param data_tbl A data frame (or tibble) with at least two columns 
#'   representing edges between nodes, and an optional column for weights.
#' @param res_level A numeric or character identifier for the current 
#'   resolution
#'
#' @return A \code{tibble} containing the Fiedler vector, partition statistics, 
#'   cluster assignments, and the resolution level.
#' 
#' @details 
#' The function follows a four-step transversal workflow:
#' \enumerate{
#'   \item Graph construction via \code{igraph}.
#'   \item Removal of self-loops.
#'   \item Calculation of the normalized Laplacian.
#'   \item Variance-based bi-partitioning.
#' }
#'
#' @importFrom igraph graph_from_data_frame delete_edges E which_loop
#' @export
#'
#' @examples
#' \dontrun{
#' # Example usage with a spatial interaction data frame
#' results <- spectral_bipartition(spatial_data, res_level = "Level_1")
#' }
spectral_bipartition <- function(data_tbl, res_level) {

  g_chr1 <- igraph::graph_from_data_frame(data_tbl, directed = FALSE)
  # eliminate self loop

  g_chr1 <- igraph::delete_edges(g_chr1, igraph::E(g_chr1)[which(igraph::which_loop(g_chr1))])
  chr_mat <- get_adj_mat_fn(g_chr1)
  # whole chromosome laplacian
  lpe_chr1 <- lp_fn(chr_mat)
  # spectral clusters
  smpl_thresh_tbl <- simple_partition_tbl_fn(lpe_chr1, res_level)

  return(smpl_thresh_tbl)
}
# %%

#' Create a New flat_node Object
#'
#' A constructor function for the \code{rc_node} S3 class. This object stores 
#' metadata, performance metrics, and structural information for a single node 
#' within a hierarchical spatial partition.
#'
#' @param id A unique character or numeric identifier for the node.
#' @param type A character string indicating the node type (e.g., 'leaf' or 'internal').
#' @param size Integer. The number of bins or observations contained within this node.
#' @param res_level Numeric. The resolution level at which this node was created.
#' @param depth Integer. The depth of the node within the hierarchy.
#' @param perf Numeric. Performance metric associated with the partition (e.g., variance explained).
#' @param parent_id The identifier of the parent node. Use \code{NA} for the root.
#' @param ids A vector of identifiers (e.g., bin IDs) belonging to this node.
#' @param null_mu Numeric. The mean of the null distribution for significance testing.
#' @param null_sd Numeric. The standard deviation of the null distribution.
#' @param obs Numeric. The observed value/statistic for this node.
#' @param df Integer. Degrees of freedom associated with the node's statistical test.
#'
#' @return An object of class \code{rc_node}, which is a structured list.
#' 
#' @export
#'
#' @examples
#' node <- new_flat_node(
#'   id = "node_1", type = "internal", size = 100, 
#'   res_level = 1, depth = 0, perf = 0.8, 
#'   parent_id = NA, ids = c(1:100), 
#'   null_mu = 0.5, null_sd = 0.1, obs = 0.7, df = 99
#' )

new_flat_node <- function(id, type, size, res_level, depth, perf, parent_id, ids, null_mu,null_sd,obs,df) {
  structure(list(id = id, type = type, size = size, res_level = res_level, depth = depth, perf = perf, parent_id= parent_id, ids = ids, null_mu = null_mu,null_sd = null_sd, obs = obs,df = df), class = "rc_node")
}
# %%

#' Fetch Nested High-Resolution Locations
#'
#' This function performs a recursive-style lookup to find high-resolution 
#' genomic bins that are contained within (nested) a set of lower-resolution bins. 
#' It is designed to work with multi-resolution files (e.g., .mcool) via the 
#' \code{hictkR} package.
#'
#' @param chrom A character string specifying the chromosome (e.g., "chr1").
#' @param ids A numeric vector of starting positions (bin IDs) at the current resolution.
#' @param res_level Integer. The current resolution index (1 being the lowest/coarsest).
#' @param high_res_level Integer. The target higher resolution index.
#' @param MresFile A list or object containing multi-resolution file metadata, 
#'   including \code{path}, \code{resolutions}, and \code{chromosomes} info.
#'
#' @return A numeric vector of starting positions representing the high-resolution 
#'   bins nested within the input locations.
#' 
#' @details 
#' The function uses \code{GenomicRanges::countOverlaps} to identify which high-resolution 
#' bins physically intersect with the range defined by the current resolution IDs. 
#' If the \code{res_level} exceeds the available resolutions in \code{MresFile}, 
#' the original \code{ids} are returned.
#'
#' @importFrom GenomicRanges GRanges start countOverlaps
#' @importFrom IRanges IRanges
#' @importFrom dplyr filter pull
#' @importFrom hictkR File fetch
#' @importFrom tibble tibble
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming mres_metadata contains the path to an .mcool file
#' high_res_bins <- fetch_nested_locations(
#'   chrom = "chr1", 
#'   ids = c(1000000, 2000000), 
#'   res_level = 1, 
#'   high_res_level = 2, 
#'   MresFile = mres_metadata
#' )
#' }
fetch_nested_locations <- function(chrom,ids,res_level,high_res_level, MresFile) {
  if (res_level >= length(MresFile$resolutions)){return(ids)} else {
  current_bin_size <- rev(MresFile$resolutions)[res_level]
  new_bin_size <- rev(MresFile$resolutions)[high_res_level]
  tmp_end <- min(c(MresFile$chromosomes|>dplyr::filter(name == chrom)|>dplyr::pull(size),max(ids)+current_bin_size))

  gr_current <-GenomicRanges::GRanges(seqnames = chrom, ranges = IRanges::IRanges(start = ids, width = current_bin_size - 1))

  new_res_hic <- hictkR::fetch(hictkR::File(MresFile$path,resolution = new_bin_size), paste0(chrom,':',min(ids),'-',tmp_end),join=TRUE)

  gr_new <- unique(GenomicRanges::GRanges(seqnames = chrom, ranges = IRanges::IRanges(start = c(new_res_hic$start1,new_res_hic$start2), end = c(new_res_hic$end1 - 1,new_res_hic$end2 -1))))


  new_bin_in_current_data <- GenomicRanges::countOverlaps(gr_new,gr_current)
  high_res_ids <- tibble::tibble(ids = GenomicRanges::start(gr_new), io = new_bin_in_current_data) |>
    dplyr::filter(io > 0) |>
    dplyr::pull(ids)
  return(high_res_ids)
  }
}

# %%

#' Fetch and Transform Interaction Matrix Data
#'
#' This function retrieves genomic interaction data from a multi-resolution file 
#' (e.g., .mcool), filters for interactions within specified genomic bins, 
#' and applies a Box-Cox transformation to the interaction counts to 
#' stabilize variance and normalize the data for downstream network analysis.
#'
#' @param chrom A character string specifying the chromosome (e.g., "chr1").
#' @param ids A numeric vector of bin start positions at the specified resolution.
#' @param res_level Integer. The index of the resolution to fetch from 
#'   (referencing \code{MresFile$resolutions}).
#' @param MresFile A list or object containing multi-resolution file metadata, 
#'   including \code{path}, \code{resolutions}, and \code{chromosomes}.
#'
#' @return A \code{tibble} containing three columns:
#' \itemize{
#'   \item \code{bin1_id}: Start position of the first bin.
#'   \item \code{bin2_id}: Start position of the second bin.
#'   \item \code{weight}: Box-Cox transformed and shifted interaction count.
#' }
#' 
#' @details 
#' The transformation process:
#' \enumerate{
#'   \item Data is fetched using \code{hictkR} with KR (Knight-Ruiz) normalization.
#'   \item Interactions are filtered to include only those where both anchors 
#'         overlap the provided \code{ids}.
#'   \item A Box-Cox power transformation is optimized using \code{MASS::boxcox}.
#'   \item Transformed counts are shifted by \code{-min(x) + 0.001} to ensure 
#'         all weights are strictly positive for graph construction.
#' }
#'
#' @importFrom hictkR File fetch
#' @importFrom GenomicRanges GRanges countOverlaps
#' @importFrom IRanges IRanges
#' @importFrom dplyr filter mutate select pull
#' @importFrom MASS boxcox
#' @importFrom stats na.omit
#' @export
#'
#' @examples
#' \dontrun{
#' # Fetch interactions for a specific region
#' interactions <- get_interaction_matrix(
#'   chrom = "chr1", 
#'   ids = seq(1000000, 2000000, by = 10000), 
#'   res_level = 2, 
#'   MresFile = mres_metadata
#' )
#' }
get_interaction_matrix <- function(chrom,ids,res_level, MresFile) {
  bin_size <-  rev(MresFile$resolutions)[res_level]
  tmp_end <- min(c(MresFile$chromosomes|>dplyr::filter(name == chrom)|>dplyr::pull(size),max(ids)+bin_size))
  data_tbl <- hictkR::fetch(hictkR::File(MresFile$path,resolution = bin_size), paste0(chrom,':',min(ids),'-',tmp_end),normalization='KR',join=TRUE)
  gr_current <-GenomicRanges::GRanges(seqnames = chrom, ranges = IRanges::IRanges(start = ids, width = bin_size))

  gr_a <- GenomicRanges::GRanges(seqnames = chrom, ranges = IRanges::IRanges(start = data_tbl$start1, end = data_tbl$end1 - 1))
  gr_b <- GenomicRanges::GRanges(seqnames = chrom, ranges = IRanges::IRanges(start = data_tbl$start2, end = data_tbl$end2 - 1))

  data_tbl <- data_tbl |>
    dplyr::mutate(a_in = GenomicRanges::countOverlaps(gr_a,gr_current))|>
    dplyr::mutate(b_in = GenomicRanges::countOverlaps(gr_b,gr_current))|>
    dplyr::filter(a_in > 0 & b_in > 0)
  data_tbl <-na.omit(data_tbl)
  if (nrow(data_tbl) > 1) {
  box_obj <- MASS::boxcox(data_tbl$count ~ 1, lambda = seq(-2, 2, 1/10),plotit=FALSE)
  exact_lambda <- box_obj$x[which.max(box_obj$y)]

  if (exact_lambda == 0) {
    transformed_data <- log(data_tbl$count)
  } else {
    transformed_data <- (data_tbl$count^exact_lambda - 1) / exact_lambda
  }
  data_tbl <- data_tbl|>
    dplyr::mutate(weight = transformed_data - min(transformed_data) + 0.001)|>
    dplyr::mutate(bin1_id = start1, bin2_id = start2)|>
    dplyr::select(bin1_id,bin2_id,weight)
  return(data_tbl)
  } else {

  return(NULL)
  }
  
}

# %%
recursive_spectral <- function(MresFile,chrom,ids, res_level = 1, depth = 1, parent_id = NA, threshold,node_registry) {
  # --- 1. Initialization and ID Generation ---
  # Infer maximum resolution from MresFile
  max_res <- length(MresFile$resolutions)
  # Generate unique ID for this attempt
  current_id <- paste0("D", depth, "_R", res_level, "_", length(ids), "_", min(ids), "_", max(ids))
  n_locs <- length(ids)

  # --- BRANCH 1: Less than 3 locations -> Go to higher resolution ---
  if (n_locs < 4 && res_level < max_res) {
    print(paste0(current_id,": higher res because too few bins"))
    # logic to select next higher and divisible resolution
    candidate_res_lvls <- (res_level + 1):max_res
    current_val <- rev(MresFile$resolutions)[res_level]
    candidate_vals <- rev(MresFile$resolutions)[candidate_res_lvls]
    is_aligned <- (current_val %% candidate_vals == 0)
    high_res_level <- candidate_res_lvls[which(is_aligned)[1]]
    high_res_ids <- fetch_nested_locations(chrom,ids, res_level,high_res_level,MresFile)
    return(recursive_spectral(MresFile,chrom,high_res_ids, high_res_level, depth, parent_id,threshold,node_registry))
  }
  # --- 3. Branch: Terminal Leaf (Smallest resolution) ---
  if (res_level == max_res && n_locs < 4){
    print(paste0(current_id,": Found to be leaf because too small"))
    node_obj <- new_flat_node(current_id, "Leaf", n_locs, res_level, depth, NA, parent_id, ids, NA,NA,NA,NA)
    assign(current_id, node_obj, envir = node_registry)
    return(if(!is.na(parent_id)) data.frame(from=parent_id, to=current_id) else data.frame())
    }
  # --- Spectral Clustering Step ---
  # Slice the adjacency matrix for the current resolution and IDs
  print(paste0(current_id,": Spectral clustering"))
  data_tbl <- get_interaction_matrix(chrom,ids,res_level,MresFile)
  if (all(data_tbl$bin1_id == data_tbl$bin2_id) || is.null(data_tbl)){

    print(paste0(current_id,": Found to be leaf because no interaction data"))
    node_obj <- new_flat_node(current_id, "Leaf", n_locs, res_level, depth, NA, parent_id, ids,NA,NA,NA,NA)
    assign(current_id, node_obj, envir = node_registry)
    return(if(!is.na(parent_id)) data.frame(from=parent_id, to=current_id) else data.frame())


  }
  spec_res <- spectral_bipartition(data_tbl,res_level)
  # Evaluate bootstrapped separability
  if (max(spec_res$stat) == 1) {
    is_separable <- TRUE
    perf <- NA
    mu_null <- NA
    sd_null <- NA
    df_null <- NA
    } else {

  boot_stat_vec <- replicate(10, {
      tryCatch({
      tmp_tbl <- spectral_bipartition(data_tbl |> mutate(weight = sample(weight)), res_level)
      max(tmp_tbl$stat)
      },error = function(e){return(NA)})
      })
  mu_null1 <- mean(boot_stat_vec,na.rm=TRUE)
  sd_null1 <- sd(boot_stat_vec,na.rm=TRUE)
  tmp_df1 <- sum(!is.na(boot_stat_vec)) - 1
  t_stat1 <- (max(spec_res$stat,na.rm=TRUE) - mu_null1) / sd_null1
  perf1 <- pt(t_stat1, df=tmp_df1, lower.tail = FALSE)

# Check for "Ambiguity Zone" (e.g., p-value between 0.01 and 0.2)
  is_ambiguous <- perf1 > 0.45 && perf1 < 0.55
  print(is_ambiguous)
  if (is_ambiguous || is.nan(is_ambiguous) || is.na(is_ambiguous)) {
    print(paste0(current_id,": ambiguous separability so running larger bootstrap"))
    # --- STAGE 2: Refinement Batch ---
    high_boot_stat_vec <- replicate(40, {
      tryCatch({
      tmp_tbl <- spectral_bipartition(data_tbl |> mutate(weight = sample(weight)), res_level)
      max(tmp_tbl$stat)
      },error = function(e){return(NA)})
      })

    # Combine results for a high-precision estimate
    final_boot <- c(boot_stat_vec, high_boot_stat_vec)
    mu_null   <- mean(final_boot, na.rm = TRUE)
    sd_null   <- sd(final_boot, na.rm = TRUE)
    df_null   <- sum(!is.na(final_boot)) - 1
    t_stat   <- (max(spec_res$stat)- mu_null) / sd_null
    perf       <- pt(t_stat, df = df_null, lower.tail = FALSE)
  } else {
    perf <- perf1
    mu_null <- mu_null1
    sd_null <- sd_null1
    df_null <- tmp_df1
  }
  is_separable <- (perf < threshold)
  }
  print(paste0(threshold,":",perf))
  # --- BRANCH 2: Is the cluster separable? ---
  # We use the 'is_separable' logic (performance < threshold)
  # --- BRANCH 3 & 4: Logic for Not Separable ---
  if (!is_separable) {
    if (res_level < max_res) {
      # PIVOT: Higher resolution
    print(paste0(current_id,": Higher resolution because not separable"))
    candidate_res_lvls <- (res_level + 1):max_res
    current_val <- rev(MresFile$resolutions)[res_level]
    candidate_vals <- rev(MresFile$resolutions)[candidate_res_lvls]
    is_aligned <- (current_val %% candidate_vals == 0)
    high_res_level <- candidate_res_lvls[which(is_aligned)[1]]
    high_res_ids <- fetch_nested_locations(chrom,ids, res_level,high_res_level,MresFile)
      return(recursive_spectral(MresFile,chrom,high_res_ids, high_res_level, depth, parent_id,threshold,node_registry))
    } else {
      # TERMINATE: Lead Node (Leaf)
      # id, type, size, perf, "res_level, depth, ids
      print(paste0(current_id,": Found to be leaf"))
      node_obj <- new_flat_node(current_id, "Leaf", n_locs, res_level, depth, perf, parent_id, ids, mu_null, sd_null, max(spec_res$stat),df_null)
      assign(current_id, node_obj, envir = node_registry)
      return(if(!is.na(parent_id)) data.frame(from=parent_id, to=current_id) else data.frame())
    }
  }

  # --- BRANCH 5: Separable -> Continue with Children ---
  node_obj <- new_flat_node(current_id, "Cluster", n_locs, res_level, depth, perf, parent_id,ids, mu_null, sd_null, max(spec_res$stat),df_null)
  assign(current_id, node_obj, envir = node_registry)
  # Track edge for this node
  current_edge <- if(!is.na(parent_id)) data.frame(from=parent_id, to=current_id) else data.frame()
  # Split IDs and Recurse to next depth
  print(paste0(current_id,": split into children clusters to process"))
  split_ids <-spec_res%>%group_by(smpl.cl)%>%summarise(ids=list(bins))%>%pull(ids)
  child_edges <- do.call(rbind, lapply(split_ids, function(child_ids) {
    recursive_spectral(MresFile, chrom, child_ids, res_level, depth + 1, current_id,threshold,node_registry)
  }))

  return(rbind(current_edge, child_edges))
}

#' Convert a BHiCect Node to a Tibble Row
#'
#' @description
#' A helper function that flattens a custom node object (S3 or list with attributes) 
#' into a single-row \code{tibble}. It automatically detects vectors or complex 
#' structures and wraps them in list-columns to ensure the resulting data frame 
#' remains flat and compliant with tidy data principles.
#'
#' @param node An object of class \code{bhicect_node} or a named list. 
#'   Expected to contain genomic bin IDs and clustering metadata.
#'
#' @return A \code{\link[tibble]{tibble}} with one row, where complex attributes 
#'   (like genomic IDs) are stored as list-columns.
#' 
#' @details 
#' This function is essential for collapsing the hierarchical \code{node_registry} 
#' (an environment) into a master registry table. It uses \code{unclass} to 
#' strip S3 methods while preserving the underlying data.
#'
#' @export
#'
#' @importFrom tibble as_tibble
#' @importFrom purrr map
#'
#' @examples
#' \dontrun{
#' node <- list(id = "D1_R1", ids = c(1, 2, 3), depth = 1)
#' node_to_row(node)
node_to_row <- function(node) {
  # 1. Capture the underlying list data
  # 2. Capture the metadata attributes
  # 3. Combine them into one flat list
  all_data <- c(unclass(node), attributes(node))
  # 5. Handle list-columns
  # If an attribute is a vector (like bin IDs), wrap it in list() 
  # so it occupies exactly one cell in the tibble row.
  tidy_data <- lapply(all_data, function(x) {
    if (length(x) > 1) return(list(x))
    return(x)
  })
 if ("ids" %in% names(tidy_data)) {
    tidy_data$ids <- list(tidy_data$ids)
  } 
  return(as_tibble(tidy_data))
}

#' Run BHiCect Hierarchical Spectral Clustering
#'
#' @description
#' The primary entry point for the BHiCect2 pipeline. This function initializes 
#' the clustering process on a specified chromosome, manages the recursive 
#' spectral decomposition, and assembles the results into a tidy format.
#'
#' @param MresFile A named list containing the \code{path} to a multi-resolution 
#'   Hi-C file (e.g., .mcool or .hic) and a numeric vector of \code{resolutions}.
#' @param chrom Character. The chromosome name to process (e.g., "chr1").
#' @param threshold Numeric. The significance threshold (p-value) for 
#'   determining cluster separability. Defaults to \code{0.5}.
#'
#' @return A named \code{list} of class \code{bhicect_res} containing:
#' \itemize{
#'   \item \code{nodes}: A \code{\link[tibble]{tibble}} registry of all identified 
#'     clusters and leaf nodes with their metadata.
#'   \item \code{edges}: A \code{data.frame} defining the parent-child 
#'     relationships in the hierarchy.
#' }
#'
#' @details 
#' The function uses an internal hash-table (environment) to store node data 
#' efficiently during recursion. It leverages \code{hictkR} for fast access 
#' to genomic bin data.
#'
#' @export
#'
#' @importFrom hictkR File
#' @importFrom dplyr filter pull
#' @importFrom cli cli_h1 cli_alert_success cli_progress_step
#' @importFrom tibble as_tibble
#' @importFrom dplyr bind_rows
#'
#' @examples
#' \dontrun{
#' mres <- list(path = "data.mcool", resolutions = c(1000, 5000, 10000))
#' results <- waterfall_bhicet(mres, chrom = "chr19", threshold = 0.05)
#' }
waterfall_bhicet <- function(MresFile,chrom,threshold=0.5){

  f <- hictkR::File(MresFile$path, resolution = rev(MresFile$resolutions)[1])
  init_ids <- f$bins|>dplyr::filter(chrom ==chrom)|>dplyr::pull(start)
  node_registry <- new.env(hash = TRUE, parent = emptyenv())

  spec_res <- recursive_spectral(MresFile,chrom,init_ids, res_level = 1, depth = 1, parent_id = NA, threshold = threshold,node_registry)
  cluster_tibble <- do.call(bind_rows,eapply(node_registry,function(x) node_to_row(x)))

  return(list(
      nodes = cluster_tibble,
      edges = spec_res
    ))
}
# %%
options(scipen=9999)
path <- "/home/vipink/Documents/BHiCeCT2/data/HCT116/4DNFIP8RKGDG.mcool"
current_MresFile <- hictkR::MultiResFile(path)
res_obj <- waterfall_bhicet(current_MresFile, "chr10", threshold = 0.5)
# %%
saveRDS(res_obj, file = "~/Documents/BHiCeCT2/data/chr10_res_obj.rds")
# %%

library(igraph)
library(hictkR)
library(MASS)
library(dplyr)
library(ggplot2)
# %%
path <- "/home/vipink/Documents/BHiCeCT2/data/HCT116/4DNFIP8RKGDG.mcool"
current_MresFile <- hictkR::MultiResFile(path)
res_obj <- readRDS("~/Documents/BHiCeCT2/data/chr10_res_obj.rds")
summary_tbl <- res_obj$nodes
summary_tbl <- summary_tbl|>mutate(parent_res_level = stringr::str_split(parent_id,'_')[[1]][2]) |> mutate(parent_res_level = as.integer(stringr::str_sub(parent_res_level,2,-1)))
# %%
summary_tbl |>
  ggplot(aes(perf,col=type))+
  geom_density()
# %%
summary_tbl |>
  ggplot(aes(depth,size))+
  geom_point()+
  scale_y_log10()

# %%

summary_tbl |>
  ggplot(aes(depth,res_level))+
  geom_point()+
  scale_y_log10()

# %%

summary_tbl |>
  mutate(lbin = log10(size))|>
  ggplot(aes(lbin,perf))+
  geom_point()+
  scale_y_log10()
# %%
# Calculate stats from your bootstrap 'null_scores'
plot_node_density_boot <- function(node_id, summary_tbl) {
# 1. Define range up to 1.0 only
  x_range <- seq(0, 1, length.out = 200)
  node <- summary_tbl |> filter(id == node_id)
  # 2. Calculate standard density for the valid range
  density_values <- (1/node$null_sd) * dt((x_range - node$null_mu)/node$null_sd, df = node$df)
  
  # 3. Calculate the "Accumulated Tail" (Area from 1 to Infinity)
  # We use pt() with lower.tail = FALSE to get the probability mass > 1
  t_stat_at_1 <- (1 - node$null_mu) / node$null_sd
  accumulated_mass <- pt(t_stat_at_1, df = node$df, lower.tail = FALSE)
  plt_df <- tibble::tibble(score = x_range, dens = density_values)
  gg <- plt_df |>
    ggplot2::ggplot(ggplot2::aes(score, dens)) +
    ggplot2::geom_line() +
    # Add a visual "spike" at 1.0 to show the accumulated mass
    ggplot2::geom_segment(ggplot2::aes(x = 1, xend = 1, y = 0, yend = max(dens)),
                          linetype = "dotted", color = "blue") +
    ggplot2::geom_point(ggplot2::aes(x = 1, y = max(dens)), color = "blue") +
    # Mark the observed score (capped at 1 for the intercept)
    ggplot2::geom_vline(xintercept = min(node$obs, 1), color = "red", linewidth = 1) +
    ggplot2::labs(
      title = paste("Censored Null Distribution: Node", node$id),
      subtitle = paste("Probability of null greater than observed:", node$perf),
      x = "Separability Score",
      y = "Density"
    ) +
    ggplot2::xlim(0, 1.1) +
    ggplot2::theme_minimal()
  return(gg)
}
# %%
find_genomic_runs <- function(ids, res_level) {
  if (length(ids) == 0) return(NULL)
  ids <- sort(unique(ids))
  thresh <- res_level
  # Find indices where the jump between positions exceeds the threshold
  breaks <- c(0, which(diff(ids) > thresh), length(ids))
  
  purrr::map(1:(length(breaks) - 1), ~{
    start_idx <- breaks[.x] + 1
    end_idx   <- breaks[.x + 1]
    # We define the rectangle from start of first bin to END of last bin
    return(list(start = ids[start_idx], end = ids[end_idx] + res_level))
  })
}
# %%
get_cl_plot_rectangles <- function(x, current_MresFile) {
  tmp_contiguous_blocks <- find_genomic_runs(unlist(x|>pull(ids)),rev(current_MresFile$resolutions)[x|>pull(res_level)])
      # Generate all pairwise combinations of runs within this cluster

  rectangle_set_df <- expand.grid(r1 = seq_along(tmp_contiguous_blocks), r2 = seq_along(tmp_contiguous_blocks)) %>%
  purrr::pmap_dfr(function(r1, r2){
  run_i <- tmp_contiguous_blocks[[r1]]
          run_j <- tmp_contiguous_blocks[[r2]]
          tibble(
            xmin = run_i$start, xmax = run_i$end,
            ymin = run_j$start, ymax = run_j$end,
          )
        })
  return(rectangle_set_df |> mutate(depth = x |> pull(depth)))
}
# %%
all_rect_tbl <- purrr::map_dfr(seq_len(summary_tbl |> nrow()), function(idx) {
        get_cl_plot_rectangles(summary_tbl |> dplyr::slice(idx),current_MresFile)
      })

# %%
ggplot(all_rect_tbl |> arrange(depth)) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = depth),
            color = NA,     # CRITICAL: No borders for 5.5k blocks
            alpha = 0.75) +   # Transparency lets nested TADs show through
  theme_minimal() +             # Removes background noise
  scale_fill_viridis_c(direction = -1, option = "magma")+
#  xlim(0,7e7)+
#  ylim(0,7e7)+
  coord_fixed()
# %%
# Utility function to convert a single cluster into an InteractionSet for integration with existing
# ecosystem
convert_to_InteractionSet <- function(node_id,summary_tbl,MresFile){




}

detect_overlap_with_feature <- function(GRange_of_interest,summary_tbl,MresFile){

}

# %%
edge_tbl <- res_obj$edge_tbl
# %%
library(ggraph)
tree_graph <- graph_from_data_frame(edge_tbl, directed = TRUE)
leaf_indices <- which(degree(tree_graph, mode = "out") == 0)
leaf_dist_mat <- distances(
  graph = tree_graph,
  v = names(leaf_indices),
  to = names(leaf_indices),
  mode = "all"
)
# %%
ggraph(tree_graph, layout = "dendrogram", circular = TRUE) +
  # Use very thin, transparent lines for massive trees
  geom_edge_diagonal(alpha = 1, width = 0.1, color = "grey50") +
  # Rasterize or shrink nodes
  geom_node_point(aes(size = 0.01,alpha=0.01)) +
  # If you want to see specific "important" labels
  # geom_node_text(aes(label = name), repel = TRUE, size = 1) +
  theme_graph() +
  theme(legend.position = "none")

# %%
long_df <- as.data.frame.table(leaf_dist_mat, responseName = "tree_distance")
