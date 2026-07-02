# BHiCect2

<!-- badges: start -->

<!-- badges: end -->

The goal of BHiCect2 is to cluster HiC data as described in our [manuscript](). <!-- markdownlint-disable-line MD042 -->
Briefly, we decompose intra-chromosomal HiC data into nested clusters of chromosome regions across multiple resolutions
starting from the complete chromosome all the way to DNA-loops at the maximum resolution provided.


---

# `BHiCect2` — High-Resolution Chromatin Domain & Cluster Detection

`BHiCect2` is an R package tailored for computational biologists and bioinformaticians looking to extract high-confidence chromatin clusters from Hi-C data. By interfacing natively with modern `.mcool` and `.hic` data formats via `hictkR`, `BHiCect2` handles multi-resolution tracking, cluster identification, and bootstrapping-based cluster significance at genomic scale.

## Key Features

* **Native Multi-Res Support:** Direct indexing of `.mcool` files via fast C++ backend interfaces (`hictkR`).
* **Interoperable Data Layouts:** Seamless conversion of hierarchical clusters into Bioconductor-native `GRangesList` objects for downstream enrichment and overlap analysis.
* **Heatmap visualisation** of the multi-resolution architecture uncovered by BHiCect
---

## Quick Start Guide

### 1. BHiCect clustering

Load your multi-resolution Hi-C data and execute the core `BHiCect` clustering algorithm over a target chromosome.

```R
library(BHiCect2)
library(hictkR)

# Prevent scientific notation on genomic coordinates
options(scipen = 9999)

# Load multi-resolution Hi-C data (.mcool format)
mcool_path <- "path/to/your/dataset.mcool"
current_MresFile <- hictkR::MultiResFile(mcool_path)

# Execute tree-based clustering algorithm on target chromosome
res_obj <- BHiCect(current_MresFile, chrom = "chr20", threshold = 0.5)

# (Optional) Persist your computed result object
saveRDS(res_obj, file = "chr20_res_obj.rds")

```

### 2. Heatmap visualisation of multi-resolution clustering

`BHiCect2` clusters are best visualised as an interaction heatmap. You can distribute this workload effortlessly using a `future` execution plan.

```R
library(future)
library(furrr)

# Extract cluster table
summary_tbl <- res_obj$nodes

# 1. Spawn parallel workers based on your infrastructure (e.g., 4 cores)
plan(multisession, workers = 4)
# 2. Compute spatial geometric coordinates in parallel
global_geometry <- compute_cluster_rectangles(summary_tbl, current_MresFile)
# 3. Terminate background processes cleanly to liberate system RAM
plan(sequential)

# Plot global contact heatmap with identified clusters
## Chromosome-wide
plot_cluster_heatmap(global_geometry)
## Interval of interest
plot_cluster_heatmap(global_geometry,xlim=c(3.5e7,4e7))


```
#### Chromosome-wide heatmap
![BHiCect2 Cluster Heatmap](man/figures/global_heatmap_example.png)

#### Region of interest
![BHiCect2 Cluster Heatmap](man/figures/heatmap_example.png)

Since this function outputs a ggplot object it can be integrated with other omics plots using plotgardener or patchwork.

### 3. Diagnostic Visualization & Density Profiling

Evaluate the geometric properties of detected clusters and analyze localized bootstrap densities for micro-domains.

```R
# Bootstrap validation plots for specific cluster targets
plot_node_density_boot("D4_R5_29_30000000_44000000", summary_tbl)
plot_node_density_boot("D13_R13_4_56396000_56399000", summary_tbl)

```
#### Robust cluster
![BHiCect2 Cluster Heatmap](man/figures/good_cluster_example.png)
#### Ambiguous cluster
![BHiCect2 Cluster Heatmap](man/figures/poor_cluster_example.png)
### 4. Downstream Integration (Bioconductor Ecosystem)

Easily transition your clustering results into standard genomic ranges for multi-omics integration.

```R
# Coerce result object into a GenomicRanges::GRangesList object
cluster_GRanges <- as_granges_list(summary_tbl, current_MresFile, chrom = "chr20")

# Ready for downstream packages like GenomicRanges, plyranges, or ChiPseeker
print(cluster_GRanges)

```

---

## Output Formats

* **`res_obj$nodes`**: A clean, parsed hierarchical data frame indexing cluster IDs, parent-child nodes, and raw coordinate metrics.
* **`cluster_GRanges`**: A standard `GRangesList` metadata framework containing distinct genomic coordinates mapping back to cluster groups.
