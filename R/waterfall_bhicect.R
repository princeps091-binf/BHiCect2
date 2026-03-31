# %%
library(igraph)
library(hictkR)
library(dplyr)
library(MASS)
# %%
get_adj_mat_fn <- function(g_chr1) {
  chr_mat <- igraph::as_adjacency_matrix(g_chr1, type = "both", attr = "weight")
  diag(chr_mat) <- 0
  if (any(Matrix::colSums(chr_mat) == 0)) {
    out <- which(Matrix::colSums(chr_mat) == 0)
    chr_mat <- chr_mat[-out, ]
    chr_mat <- chr_mat[, -out]
  }
  return(chr_mat)
}

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

ss <- function(x) {

  sum(scale(x, scale = FALSE)^2)
}

simple_partition_tbl_fn <- function(lp_res, tmp_res) {
  smpl_thresh_tbl <- lp_res$vectors |>
    dplyr::mutate(real_fiedler = Re(fiedler))|>
    dplyr::mutate(stat = purrr::map_dbl(real_fiedler, function(x) {
      cl_a <- real_fiedler[which(real_fiedler <= x)]
      cl_b <- real_fiedler[which(real_fiedler > x)]
      return(ss(rep(c(mean(cl_a), mean(cl_b)), c(length(cl_a), length(cl_b)))) / ss(real_fiedler))
    }))
  smpl_thresh <- smpl_thresh_tbl |>
    dplyr::slice_max(stat) |>
    dplyr::pull(real_fiedler)
  smpl_thresh_tbl <- smpl_thresh_tbl |>
    dplyr::mutate(
      smpl.cl = ifelse(real_fiedler <= smpl_thresh, 1, 2),
      res = tmp_res
    )
  smpl_thresh_tbl
}

# %%
spectral_bipartition <- function(data_tbl, res_level) {

  g_chr1 <- igraph::graph_from_data_frame(data_tbl, directed = FALSE)
  # eliminate self loop
  g_chr1 <- igraph::delete_edges(g_chr1, E(g_chr1)[which(igraph::which_loop(g_chr1))])
  chr_mat <- get_adj_mat_fn(g_chr1)
  # whole chromosome laplacian
  lpe_chr1 <- lp_fn(chr_mat)
  # spectral clusters
  smpl_thresh_tbl <- simple_partition_tbl_fn(lpe_chr1, res_level)

  return(smpl_thresh_tbl)
}
# %%
# current_id, "Leaf", n_locs, res_level, depth, perf, parent_id, ids
new_flat_node <- function(id, type, size, res_level, depth, perf, parent_id, ids, null_mu,null_sd,obs,df) {
  structure(list(id = id, type = type, size = size, res_level = res_level, depth = depth, perf = perf, parent_id= parent_id, ids = ids, null_mu = null_mu,null_sd = null_sd, obs = obs,df = df), class = "rc_node")
}
# %%
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
}

# %%
recursive_spectral <- function(MresFile,chrom,ids, res_level = 1, depth = 1, parent_id = NA, threshold,node_registry) {
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
  if (all(data_tbl$bin1_id == data_tbl$bin2_id)){

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
  if (is_ambiguous | is.nan(is_ambiguous) | is.na(is_ambiguous)) {
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
  #perf <- (sum(boot_stat_vec > max(spec_res$stat))/100)
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
res_obj <- waterfall_bhicet(current_MresFile, "chr19", threshold = 0.5)
# %%
saveRDS(res_obj, file = "~/Documents/BHiCeCT2/data/chr19_res_obj.rds")
# %%

library(igraph)
library(hictkR)
library(MASS)
library(dplyr)
library(ggplot2)
# %%
path <- "/home/vipink/Documents/BHiCeCT2/data/HCT116/4DNFIP8RKGDG.mcool"
current_MresFile <- hictkR::MultiResFile(path)
res_obj <- readRDS("~/Documents/BHiCeCT2/data/chr19_res_obj.rds")
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
        get_cl_plot_rectangles(summary_tbl |> slice(idx),current_MresFile)
        })
# %%
ggplot(all_rect_tbl |> arrange(depth)) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = depth),
            color = NA,     # CRITICAL: No borders for 5.5k blocks
            alpha = 0.75) +   # Transparency lets nested TADs show through
  theme_minimal() +             # Removes background noise
  scale_fill_viridis_c(direction = -1, option = "magma")+
  xlim(0,59000000)+
  ylim(0,59000000)+
  coord_fixed()
# %%
# %%
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
  geom_node_point(aes(size = 0.01)) +
  # If you want to see specific "important" labels
  # geom_node_text(aes(label = name), repel = TRUE, size = 1) +
  theme_graph() +
  theme(legend.position = "none")
