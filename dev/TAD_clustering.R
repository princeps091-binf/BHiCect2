# %%
library(devtools)

devtools::load_all()

# %%

library(hictkR)
library(dplyr)

BHiCect_Locus <- function(MresFile, chrom, tad_start, tad_end, start_res, threshold = 0.5) {
  message(sprintf(
    "Targeted Analysis: Isolating matrix window for TAD %s:%d-%d at %d kb",
    chrom, tad_start, tad_end, start_res / 1000
  ))

  # 1. Identify the resolution level index matching the input TAD resolution
  # We search your file's resolution vector (sorted from coarsest to finest)
  all_res <- rev(MresFile$resolutions)
  start_lvl <- which(all_res == start_res)

  if (length(start_lvl) == 0) {
    stop(paste(
      "Error: The specified tad_res (", start_res,
      ") does not exist in your multi-resolution file."
    ))
  }

  # 2. Open the file directly at the native TAD resolution layer
  f <- hictkR::File(MresFile$path, resolution = start_res)

  # 3. Restrict initialization to only the bins contained inside the TAD boundaries
  init_ids <- f$bins |>
    dplyr::filter(chrom == !!chrom & start >= tad_start & start <= tad_end) |>
    dplyr::pull(start)

  if (length(init_ids) == 0) {
    stop("Execution halted: No valid genomic bins were found within the provided coordinates at this resolution.")
  }

  # 4. Initialize the identical environment tracking registry
  node_registry <- new.env(hash = TRUE, parent = emptyenv())
  print(length(init_ids))
  # 5. Invoke the core recursive engine, forcing the root depth to 1 and tagging the external parent
  spec_res <- recursive_spectral(
    MresFile = MresFile,
    chrom = chrom,
    ids = init_ids,
    res_level = start_lvl,
    depth = 1,
    parent_id = "ORIGINAL_INTERVAL",
    threshold = threshold,
    node_registry = node_registry
  )

  # 6. Reconstruct the standard node table layout
  cluster_tibble <- do.call(bind_rows, eapply(node_registry, function(x) node_to_row(x)))

  return(list(
    nodes = cluster_tibble,
    edges = spec_res
  ))
}

# %%

library(readr)

convert_4dn_boundaries_to_tads <- function(boundary_bed_path, min_tad_size = 40000, max_tad_size = 2000000) {
  # 1. Read the 4DN boundary BED file
  # 4DN boundaries usually have columns: chrom, start, end
  boundaries <- readr::read_tsv(boundary_bed_path,
    col_names = c("chrom", "start", "end", "name", "score", "strand"),
    comment = "#"
  )

  # 2. Convert boundary points into adjacent interval blocks
  tad_intervals <- boundaries |>
    dplyr::arrange(chrom, start) |>
    dplyr::group_by(chrom) |>
    dplyr::mutate(
      # The start of the TAD is the end of the previous boundary point
      tad_start = end,
      # The end of the TAD is the start of the next boundary point
      tad_end = lead(start),
      left_border = score,
      right_border = lead(score)
    ) |>
    dplyr::ungroup() |>
    # Remove rows where a trailing neighbor doesn't exist (ends of chromosomes)
    dplyr::filter(!is.na(tad_end)) |>
    # 3. Clean and filter based on biological size thresholds
    dplyr::mutate(tad_size = tad_end - tad_start) |>
    dplyr::filter(tad_size >= min_tad_size & tad_size <= max_tad_size) |>
    dplyr::select(chrom, start = tad_start, end = tad_end, tad_size, left_border, right_border) |>
    dplyr::mutate(insulation = (left_border + right_border) / 2)

  return(tad_intervals)
}

# %%
# This interval based analysis can be centered on genes of interest too!!!
options(scipen = 9999)
tad_tbl <- convert_4dn_boundaries_to_tads("~/Documents/BHiCeCT2/data/GM12878/4DNFIX26B8E9.bed.gz")


ctrl_path <- "/home/vipink/Documents/BHiCeCT2/data/auxin_expt/ctrl/4DNFINIQYFKT.mcool"
ctrl_MresFile <- hictkR::MultiResFile(ctrl_path)

auxin_path <- "/home/vipink/Documents/BHiCeCT2/data/auxin_expt/auxin/4DNFIAVRY6RG.mcool"
auxin_MresFile <- hictkR::MultiResFile(auxin_path)

# %%
tmp_tad <- tad_tbl |>
  arrange(desc(insulation), desc(tad_size)) |>
  slice(4)
TAD_res <- BHiCect_Locus(GM12878_MresFile, tmp_tad$chrom, tmp_tad$start, tmp_tad$end, 5000)

TAD_summary_tbl <- TAD_res$nodes

TAD_summary_tbl <- TAD_summary_tbl |>
  mutate(start = str_split(id, "_") |> purrr::map_int(\(x)as.integer(x[4]))) |>
  mutate(end = str_split(id, "_") |> purrr::map_int(\(x)as.integer(x[5]))) |>
  mutate(chrom = tmp_tad |> dplyr::pull(chrom))

root_id <- TAD_summary_tbl |>
  filter(parent_id == "ORIGINAL_INTERVAL") |>
  pull(id)

plot_exhaustive_locked_split(root_id, TAD_summary_tbl, GM12878_MresFile)

# %%
# MYC example chr8:127_735_434-127_742_951
# Target range -> 500Kb arround:
# MYC "chr8",127235000,128243000 (HCT116 active but GM12878 inactive)
# FOS "chr14",75200000,75350000
# HOX "chr7",27100000,27250000
# CDX2 "chr13", 28450000,28650000
# POU2F2 (Oct-2) "chr19", 48900000, 49150000
# IGF2 - H19 "chr11", 1950000, 2150000
# CCND1 (Cyclin D1) "chr11", 69000000, 70000000 -> HCT116 upregulated (cell cycling)
# Bcl2l11 (BIM) Locus "chr2", 111000000, 112200000 -> GM12878 upregulated (pro-apoptotic)
ctrl_res <- BHiCect_Locus(ctrl_MresFile, "chr8", 127235000, 128243000, 25000)
ctrl_summary_tbl <- ctrl_res$nodes

auxin_res <- BHiCect_Locus(auxin_MresFile, "chr8", 127235000, 128243000, 25000)
auxin_summary_tbl <- auxin_res$nodes

# %%
ctrl_geometry <- compute_cluster_rectangles(ctrl_summary_tbl, ctrl_MresFile)
auxin_geometry <- compute_cluster_rectangles(auxin_summary_tbl, auxin_MresFile)

pctrl <- plot_cluster_heatmap(ctrl_geometry)
pauxin <- plot_cluster_heatmap(auxin_geometry)

patchwork::wrap_plots(c(pctrl, pauxin), nrow = 1)
# %%
current_MresFile <- auxin_MresFile
summary_tbl <- auxin_summary_tbl |>
  mutate(start = str_split(id, "_") |> purrr::map_int(\(x)as.integer(x[4]))) |>
  mutate(end = str_split(id, "_") |> purrr::map_int(\(x)as.integer(x[5]))) |>
  mutate(chrom = "chr8")

root_id <- summary_tbl |>
  filter(parent_id == "ORIGINAL_INTERVAL") |>
  pull(id)

plot_exhaustive_locked_split(root_id, summary_tbl, current_MresFile, "KR")
