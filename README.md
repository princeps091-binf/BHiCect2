# BHiCect2

<!-- badges: start -->

<!-- badges: end -->

The goal of BHiCect2 is to cluster HiC data as described in our [manuscript](). <!-- markdownlint-disable-line MD042 -->
Briefly, we decompose intra-chromosomal HiC data into nested clusters of chromosome regions across multiple resolutions
starting from the complete chromosome all the way to DNA-loops at the maximum resolution provided.

## Installation

You can install the current version of BHiCect2 from [GitHub](https://github.com/) with:

```r
# install.packages("devtools")
devtools::install_github("princeps091-binf/BHiCect2")
```

## Example

BHiCect2 offers a core function to cluster the input HiC data.
The expected input data to BHiCect2 is a list of dataframes containing the HiC data in a three columns format for the various resolution provided by the user.

```r
library(BHiCect2)
library(future)
library(furrr)

options(scipen = 9999)
path <- "/home/vipink/Documents/BHiCeCT2/data/HCT116/4DNFIP8RKGDG.mcool"
current_MresFile <- hictkR::MultiResFile(path)
res_obj <- BHiCect(current_MresFile, "chr20", threshold = 0.5)

saveRDS(res_obj, file = "~/Documents/BHiCeCT2/data/chr20_res_obj.rds")

path <- "/home/vipink/Documents/BHiCeCT2/data/HCT116/4DNFIP8RKGDG.mcool"
current_MresFile <- hictkR::MultiResFile(path)
res_obj <- readRDS("~/Documents/BHiCeCT2/data/chr20_res_obj.rds")

summary_tbl <- res_obj$nodes

# 1. Set up the background parallel workers (e.g., using 4 cores)
plan(multisession, workers = 4)

# 2. Run the function (It will instantly execute in parallel)
global_geometry <- compute_cluster_rectangles(summary_tbl, current_MresFile)

# 3. Explicitly close the background connections when done to free up RAM
plan(sequential)

plot_cluster_heatmap(global_geometry)

plot_node_density_boot("D4_R5_29_30000000_44000000", summary_tbl)


plot_node_density_boot("D13_R13_4_56396000_56399000", summary_tbl)


cluster_GRanges <- as_granges_list(summary_tbl, current_MresFile, "chr20")

```
