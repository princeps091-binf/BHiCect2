library(devtools)

devtools::load_all()

options(scipen = 9999)
path <- "/home/vipink/Documents/BHiCeCT2/data/HCT116/4DNFIP8RKGDG.mcool"
current_MresFile <- hictkR::MultiResFile(path)
res_obj <- waterfall_bhicet(current_MresFile, "chr20", threshold = 0.5)
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
res_obj <- readRDS("~/Documents/BHiCeCT2/data/chr10_res_obj.rds")
summary_tbl <- res_obj$nodes
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
##
find_genomic_runs <- function(ids, res_level) {
        if (length(ids) == 0) {
                return(NULL)
        }
        ids <- sort(unique(ids))
        thresh <- res_level
        # Find indices where the jump between positions exceeds the threshold
        breaks <- c(0, which(diff(ids) > thresh), length(ids))

        purrr::map(1:(length(breaks) - 1), ~ {
                start_idx <- breaks[.x] + 1
                end_idx <- breaks[.x + 1]
                # We define the rectangle from start of first bin to END of last bin
                return(list(start = ids[start_idx], end = ids[end_idx] + res_level))
        })
}
# %%
get_cl_plot_rectangles <- function(x, current_MresFile) {
        tmp_contiguous_blocks <- find_genomic_runs(unlist(x |> pull(ids)), rev(current_MresFile$resolutions)[x |> pull(res_level)])
        # Generate all pairwise combinations of runs within this cluster

        rectangle_set_df <- expand.grid(r1 = seq_along(tmp_contiguous_blocks), r2 = seq_along(tmp_contiguous_blocks)) %>%
                purrr::pmap_dfr(function(r1, r2) {
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
        get_cl_plot_rectangles(summary_tbl |> dplyr::slice(idx), current_MresFile)
})

# %%
ggplot(all_rect_tbl |> arrange(depth)) +
        geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = depth),
                color = NA, # CRITICAL: No borders for 5.5k blocks
                alpha = 0.75
        ) + # Transparency lets nested TADs show through
        theme_minimal() + # Removes background noise
        scale_fill_viridis_c(direction = -1, option = "magma") +
        #  xlim(0,7e7)+
        #  ylim(0,7e7)+
        coord_fixed()

# %%
plot_node_density_boot <- function(node_id, summary_tbl) {
        # 1. Define range up to 1.0 only
        x_range <- seq(0, 1, length.out = 200)
        node <- summary_tbl |> filter(id == node_id)
        # 2. Calculate standard density for the valid range
        density_values <- (1 / node$null_sd) * dt((x_range - node$null_mu) / node$null_sd, df = node$df)

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
                        linetype = "dotted", color = "blue"
                ) +
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
