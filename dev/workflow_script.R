library(devtools)

devtools::load_all()

options(scipen = 9999)
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
summary_tbl <- summary_tbl |>
        mutate(parent_res_level = stringr::str_split(parent_id, "_")[[1]][2]) |>
        mutate(parent_res_level = as.integer(stringr::str_sub(parent_res_level, 2, -1)))
#
