# %%
library(devtools)

devtools::load_all()

# %%

library(hictkR)
library(dplyr)

BHiCect_TAD <- function(MresFile, chrom, tad_start, tad_end, tad_res, threshold = 0.5) {
  
  message(sprintf("Targeted Analysis: Isolating matrix window for TAD %s:%d-%d at %d kb", 
                  chrom, tad_start, tad_end, tad_res / 1000))
  
  # 1. Identify the resolution level index matching the input TAD resolution
  # We search your file's resolution vector (sorted from coarsest to finest)
  all_res   <- rev(MresFile$resolutions) 
  start_lvl <- which(all_res == tad_res)
  
  if (length(start_lvl) == 0) {
    stop(paste("Error: The specified tad_res (", tad_res, 
               ") does not exist in your multi-resolution file."))
  }
  
  # 2. Open the file directly at the native TAD resolution layer
  f <- hictkR::File(MresFile$path, resolution = tad_res)
  
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
    MresFile      = MresFile, 
    chrom         = chrom, 
    ids      = init_ids, 
    res_level     = start_lvl, 
    depth         = 1, 
    parent_id     = "EXTERNAL_TAD_CALLER", 
    threshold     = threshold,
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
                                comment = "#")
  
  # 2. Convert boundary points into adjacent interval blocks
  tad_intervals <- boundaries |> 
    dplyr::arrange(chrom, start) |> 
    dplyr::group_by(chrom) |> 
    dplyr::mutate(
      # The start of the TAD is the end of the previous boundary point
      tad_start = end,
      # The end of the TAD is the start of the next boundary point
      tad_end = lead(start)
    ) |> 
    dplyr::ungroup() |> 
    # Remove rows where a trailing neighbor doesn't exist (ends of chromosomes)
    dplyr::filter(!is.na(tad_end)) |> 
    # 3. Clean and filter based on biological size thresholds
    dplyr::mutate(tad_size = tad_end - tad_start) |> 
    dplyr::filter(tad_size >= min_tad_size & tad_size <= max_tad_size) |> 
    dplyr::select(chrom, start = tad_start, end = tad_end, tad_size)
  
  return(tad_intervals)
}

# %%

options(scipen = 9999)
tad_tbl <- convert_4dn_boundaries_to_tads("~/Documents/BHiCeCT2/data/GM12878/4DNFIX26B8E9.bed.gz")


path <- "/home/vipink/Documents/BHiCeCT2/data/GM12878/4DNFI2ZUCIHD.mcool"
current_MresFile <- hictkR::MultiResFile(path)
tmp_tad <- tad_tbl|>arrange(desc(tad_size))|>slice(10)
BHiCect_TAD(current_MresFile,tmp_tad$chrom,tmp_tad$start,tmp_tad$end,50000)
