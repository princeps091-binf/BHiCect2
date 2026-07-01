library(devtools)

devtools::load_all()

options(scipen = 9999)
path <- "/home/vipink/Documents/BHiCeCT2/data/HCT116/4DNFIP8RKGDG.mcool"
current_MresFile <- hictkR::MultiResFile(path)
res_obj <- waterfall_bhicect(current_MresFile, "chr20", threshold = 0.5)
# %%
saveRDS(res_obj, file = "~/Documents/BHiCeCT2/data/chr20_res_obj.rds")
# %%

library(igraph)
library(hictkR)
library(MASS)
library(dplyr)
library(ggplot2)
# %%
path <- "/home/vipink/Documents/BHiCeCT2/data/HCT116/4DNFIP8RKGDG.mcool"
current_MresFile <- hictkR::MultiResFile(path)
res_obj <- readRDS("~/Documents/BHiCeCT2/data/chr20_res_obj.rds")

summary_tbl <- res_obj$nodes

# %%

plot_cluster_heatmap(summary_tbl, current_MresFile)

# %%

plot_node_density_boot("D4_R5_29_30000000_44000000", summary_tbl)

# %%

plot_node_density_boot("D13_R13_4_56396000_56399000", summary_tbl)

# %%

cluster_GRanges <- as_granges_list(summary_tbl, current_MresFile, "chr20")


# %%
summary_tbl <- summary_tbl |>
        mutate(parent_res_level = stringr::str_split(parent_id, "_")[[1]][2]) |>
        mutate(parent_res_level = as.integer(stringr::str_sub(parent_res_level, 2, -1)))
summary_tbl |>
        ggplot(aes(perf, col = type)) +
        geom_density()
# %%
summary_tbl |>
        ggplot(aes(depth, size)) +
        geom_point() +
        scale_y_log10()

# %%

summary_tbl |>
        ggplot(aes(depth, res_level)) +
        geom_point() +
        scale_y_log10()

# %%

summary_tbl |>
        mutate(lbin = log10(size)) |>
        ggplot(aes(lbin, perf)) +
        geom_point() +
        scale_y_log10()
