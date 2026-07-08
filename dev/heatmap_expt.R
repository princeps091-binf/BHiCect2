# %%
library(hictkR)
library(dplyr)
library(ggplot2)
library(patchwork)
library(stringr)
library(devtools)

devtools::load_all()

options(scipen = 9999)

# %%
plot_multiresolution_heatmaps <- function(node_id, summary_tbl, MresFile, data_value) {
  
  # 1. Isolate parent and its direct child nodes
  parent_data <- summary_tbl |> dplyr::filter(id == node_id)
  
  current_parents <- node_id
  all_descendants <- list()
  while (length(current_parents) > 0) {
    next_generation <- summary_tbl |> 
      dplyr::filter(parent_id %in% current_parents)
    
    if (nrow(next_generation) > 0) {
      all_descendants[[length(all_descendants) + 1]] <- next_generation
      current_parents <- next_generation$id
    } else {
      current_parents <- c() # Stop when we hit terminal leaves across all branches
    }
  }
  
  if (length(all_descendants) == 0) {
    stop(paste0("Error: Node '", node_id, "' has no descendants; it is a terminal Leaf node."))
  }
  
  # Combine generations into a single comprehensive tracking table
  descendants_tbl <- dplyr::bind_rows(all_descendants)
  descendant_rect_tbl <- compute_cluster_rectangles(descendants_tbl,MresFile)
  # 2. Establish GLOBAL, STRICT viewport coordinate limits based on the parent node
  chrom <- parent_data$chrom
  global_start <- parent_data$start
  global_end   <- parent_data$end
  query_range  <- paste0(chrom, ":", global_start, "-", global_end)
  
  parent_res <- rev(MresFile$resolutions)[parent_data$res_level]
  
  # 3. Plot 1: The Parent Baseline View
  matrix_parent <- hictkR::fetch(hictkR::File(MresFile$path, resolution = parent_res),normalization = data_value, query_range,join=TRUE)
  p_parent <- gg_heatmap_locked(matrix_parent, parent_res, 
                                 title = paste("Parent Base Res:", parent_res / 1000, "kb"),
                                 x_lims = c(global_start, global_end)) +
    # Draw Parent Boundary Bounding Box
    geom_rect(aes(xmin = global_start, xmax = global_end, ymin = global_start, ymax = global_end), 
              color = "red", fill = NA, size = 1) +
    # Overlay all child nodes color-coded by their resolution level
    geom_rect(data = descendant_rect_tbl, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, color = depth), 
              fill = NA, linewidth = 0.8, inherit.aes = FALSE) +
    labs(color = "Child depth") +
    scale_color_viridis_c(option = "viridis", na.value = "transparent") 
  # 4. Generate High-Res Tier Plots with Locked Canvas Viewports
  unique_descendant_res_levels <- unique(descendant_rect_tbl$res_level) |> sort()
  descendant_plots <- list()
  
  for (res_lvl in unique_descendant_res_levels) {
    current_child_res <- rev(MresFile$resolutions)[res_lvl]
    tier_nodes <- descendant_rect_tbl |> dplyr::filter(res_level == res_lvl)
    # Fetch data using the exact same global genomic range
    matrix_tier <- hictkR::fetch(hictkR::File(MresFile$path, resolution = current_child_res), normalization = data_value, query_range,join=TRUE)
    
    # Render with the strict global coordinate limits enforced
    p_tier <- gg_heatmap_locked(matrix_tier, current_child_res, 
                                 title = paste("Child Tier Res:", current_child_res / 1000, "kb"),
                                 x_lims = c(global_start, global_end)) +
      geom_rect(data = tier_nodes, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,color = depth), 
                 fill = NA, linewidth = 0.8, inherit.aes = FALSE) +
    scale_color_viridis_c(option = "viridis", na.value = "transparent") 

    
    descendant_plots[[as.character(current_child_res)]] <- (p_parent | p_tier)
  }
  
  # 5. Composite Assembly via Patchwork
  children_layout <- patchwork::wrap_plots(descendant_plots, ncol = 1)
  
  composite_plot <- children_layout + 
    plot_annotation(title = paste("Exhaustive Multi-Resolution Sub-Tree Decomposition:", node_id),
                    subtitle = paste("Strict Coordinate Lock Area:", query_range))
  
  return(composite_plot)
}

# Core Helper: Forces absolute identity boundaries across coordinate grids
gg_heatmap_locked <- function(df, bin_size, title, x_lims) {
  # 1. Check if data frame is empty to prevent crashes
  if (nrow(df) == 0) {
    df_full <- df
  } else {
    # 2. Mirror the upper triangle to the lower triangle
    df_mirrored <- df |> 
      dplyr::rename(start1 = start2, start2 = start1,
                    end1 = end2, end2 = end1)
    
    # 3. Combine them together, ensuring unique pixel coordinates
    df_full <- dplyr::bind_rows(df, df_mirrored) |> 
      dplyr::distinct(start1, start2, .keep_all = TRUE)
  }
  ggplot(df_full, aes(x = start1, y = start2, fill = count)) +
    geom_tile() +
    scale_fill_viridis_c(option = "magma", trans = "log10", na.value = "transparent") +
    # CRITICAL FIX: Enforces identical canvas boundaries regardless of underlying matrix bin variations
    coord_cartesian(xlim = x_lims, ylim = x_lims, expand = FALSE) +
    theme_minimal() +
    labs(title = title, x = "Genomic Position (bp)", y = "Genomic Position (bp)") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          aspect.ratio = 1) # Force rigid square metrics globally
}

# %%

path <- "/home/vipink/Documents/BHiCeCT2/data/HCT116/4DNFIP8RKGDG.mcool"
HCT116_MresFile <- hictkR::MultiResFile(path)
HCT116_res_obj <- readRDS("~/Documents/BHiCeCT2/data/chr8_res_obj.rds")

# %%
HCT116_summary_tbl <- HCT116_res_obj$nodes
HCT116_summary_tbl <- HCT116_summary_tbl |>
	mutate(start = str_split(id,'_')|> purrr::map_int(\(x)as.integer(x[4])))|>
	mutate(end = str_split(id,'_')|> purrr::map_int(\(x)as.integer(x[5])))|>
	mutate(chrom = 'chr8')
# %%




# %%

HCT116_summary_tbl|>filter(size > 20 & depth < 8)|>arrange(desc(size))|>dplyr::select(id,res_level,size,perf,depth)
tmp_node = "D5_R8_263_56250000_69950000"
plot_exhaustive_locked_split(tmp_node, HCT116_summary_tbl,HCT116_MresFile)
