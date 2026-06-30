#' Plot Cluster Rectangles on Hi-C Contact Map Workspace
#'
#' Generates a ggplot2 visualization of hierarchical, non-contiguous cluster blocks.
#' If an interaction data frame or matrix is supplied, it overlays the clusters
#' directly on top of the interaction heatmap.
#'
#' @param summary_tbl A tibble or data.frame containing the cluster summaries.
#'   Must have columns: \code{ids}, \code{res_level}, and \code{depth}.
#' @param current_MresFile A list or object tracking the resolution vector
#'   (e.g., \code{current_MresFile$resolutions}).
#' @param xlim Numeric vector of length 2. The genomic window to visualize.
#'   If \code{NULL}, automatically scales to the minimum and maximum coordinates in \code{summary_tbl}.
#' @param hic_df Optional data.frame containing pre-filtered contact matrix data
#'   with columns \code{binA}, \code{binB}, and \code{count} for background plotting.
#' @param alpha Numeric. Transparency value for cluster rectangles (default: 0.75).
#'
#' @return A \code{ggplot2} object.
#' @import ggplot2
#' @import dplyr
#' @importFrom purrr map map_dfr pmap_dfr
#' @export
plot_cluster_heatmap <- function(summary_tbl,
                                 current_MresFile,
                                 xlim = NULL,
                                 hic_df = NULL,
                                 alpha = 0.75) {
        # ==========================================================================
        # 1. INTERNAL HELPERS (Hidden from user namespace)
        # ==========================================================================

        find_genomic_runs <- function(ids, res_level) {
                if (length(ids) == 0) {
                        return(NULL)
                }
                ids <- sort(unique(ids))
                breaks <- c(0, which(diff(ids) > res_level), length(ids))

                purrr::map(1:(length(breaks) - 1), ~ {
                        start_idx <- breaks[.x] + 1
                        end_idx <- breaks[.x + 1]
                        list(start = ids[start_idx], end = ids[end_idx] + res_level)
                })
        }

        get_cl_plot_rectangles <- function(x, mres) {
                # Dynamically extract resolution level mapping safely
                res_vec <- rev(mres$resolutions)
                current_res <- res_vec[x |> dplyr::pull(res_level)]

                tmp_contiguous_blocks <- find_genomic_runs(unlist(x |> dplyr::pull(ids)), current_res)
                if (is.null(tmp_contiguous_blocks)) {
                        return(dplyr::tibble())
                }

                rectangle_set_df <- expand.grid(
                        r1 = seq_along(tmp_contiguous_blocks),
                        r2 = seq_along(tmp_contiguous_blocks)
                ) %>%
                        purrr::pmap_dfr(function(r1, r2) {
                                run_i <- tmp_contiguous_blocks[[r1]]
                                run_j <- tmp_contiguous_blocks[[r2]]
                                dplyr::tibble(
                                        xmin = run_i$start, xmax = run_i$end,
                                        ymin = run_j$start, ymax = run_j$end
                                )
                        })

                return(rectangle_set_df |> dplyr::mutate(depth = x |> dplyr::pull(depth)))
        }

        # ==========================================================================
        # 2. DATA PROCESSING
        # ==========================================================================

        # Flatten tree summary table to explicit grid coordinates
        all_rect_tbl <- purrr::map_dfr(seq_len(nrow(summary_tbl)), function(idx) {
                get_cl_plot_rectangles(summary_tbl |> dplyr::slice(idx), current_MresFile)
        }) |> dplyr::arrange(depth)

        # Dynamically determine viewport limits if not provided
        if (is.null(xlim)) {
                xlim <- c(min(all_rect_tbl$xmin), max(all_rect_tbl$xmax))
        }
        ylim <- xlim # Maintain strict symmetry for Hi-C matrix diagonals

        # Filter rectangles to visible viewport to speed up rendering
        visible_rects <- all_rect_tbl |>
                dplyr::filter(xmax >= xlim[1] & xmin <= xlim[2] & ymax >= ylim[1] & ymin <= ylim[2])

        # ==========================================================================
        # 3. PLOT CONSTRUCTION
        # ==========================================================================

        p <- ggplot2::ggplot()

        # Layer A: Background Hi-C Heatmap (If supplied by user)
        if (!is.null(hic_df)) {
                # Ensure background matrix columns exist
                if (all(c("binA", "binB", "count") %in% colnames(hic_df))) {
                        p <- p + ggplot2::geom_raster(
                                data = hic_df,
                                ggplot2::aes(x = binA, y = binB, fill = count),
                                show.legend = TRUE
                        ) +
                                # Adjust scale dynamically so it doesn't conflict with cluster depth fill
                                ggplot2::scale_fill_gradient(low = "white", high = "gray20", trans = "log10") +
                                # Allow support for dual-fill scales via ggnewscale package if preferred later
                                ggnewscale::new_scale_fill()
                }
        }

        # Layer B: Non-Contiguous Cluster Overlay
        p <- p +
                ggplot2::geom_rect(
                        data = visible_rects,
                        ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = depth),
                        color = NA,
                        alpha = alpha
                ) +
                ggplot2::scale_fill_viridis_c(direction = -1, option = "magma", name = "Tree Depth") +
                ggplot2::xlim(xlim[1], xlim[2]) +
                ggplot2::ylim(ylim[1], ylim[2]) +
                ggplot2::coord_fixed() +
                ggplot2::theme_minimal() +
                ggplot2::labs(
                        x = "Genomic Coordinate (bp)",
                        y = "Genomic Coordinate (bp)"
                )

        return(p)
}


#' Plot Node Separability Against Bootstrapped Null Distribution
#'
#' Generates a diagnostic density plot for a specific tree node, displaying
#' the empirical or fitted student-t null distribution against the observed score.
#' The function accounts for metric censoring by capturing probability mass
#' exceeding 1.0 as an accumulated visual spike.
#'
#' @param node_id Character. The unique identifier of the node to visualize (e.g., "root.1").
#' @param summary_tbl A data.frame or tibble containing node statistics. Must contain columns:
#'   \code{id}, \code{null_mu}, \code{null_sd}, \code{df}, \code{obs}, and \code{perf}.
#' @param fill_tail Logical. If TRUE, shades the area of the null distribution that is
#'   greater than or equal to the observed score (visualizing the empirical p-value).
#'
#' @return A \code{ggplot2} object.
#' @import ggplot2
#' @importFrom dplyr filter
#' @importFrom tibble tibble
#' @export
plot_node_density_boot <- function(node_id, summary_tbl, fill_tail = TRUE) {
        # ==========================================================================
        # 1. SAFEGUARD AND EXTRACTION
        # ==========================================================================
        node <- summary_tbl |> dplyr::filter(id == node_id)

        if (nrow(node) == 0) {
                stop(paste("Node ID '", node_id, "' not found in the provided summary table."))
        }
        if (nrow(node) > 1) {
                warning(paste("Multiple rows found for Node ID '", node_id, "'. Using the first match."))
                node <- node[1, ]
        }

        # ==========================================================================
        # 2. MATHEMATICAL MODELING (Censored Range)
        # ==========================================================================
        x_range <- seq(0, 1, length.out = 200)

        # Standard scale transformation for Student-t density calculation
        density_values <- (1 / node$null_sd) * stats::dt((x_range - node$null_mu) / node$null_sd, df = node$df)

        # Calculate the "Accumulated Tail" mass (Area from 1 to Infinity)
        t_stat_at_1 <- (1 - node$null_mu) / node$null_sd
        accumulated_mass <- stats::pt(t_stat_at_1, df = node$df, lower.tail = FALSE)

        # Build plotting layout
        plt_df <- tibble::tibble(score = x_range, dens = density_values)
        max_dens <- max(density_values, na.rm = TRUE)

        # ==========================================================================
        # 3. GRAPH PLOTTING
        # ==========================================================================
        p <- ggplot2::ggplot(plt_df, ggplot2::aes(x = score, y = dens))

        # Optional: Shade the p-value tail area for visual diagnostics
        if (fill_tail) {
                obs_capped <- min(node$obs, 1)
                tail_df <- plt_df |> dplyr::filter(score >= obs_capped)

                if (nrow(tail_df) > 0) {
                        p <- p + ggplot2::geom_area(
                                data = tail_df,
                                ggplot2::aes(x = score, y = dens),
                                fill = "firebrick",
                                alpha = 0.15
                        )
                }
        }

        p <- p +
                # Continuous Null Distribution Curve
                ggplot2::geom_line(color = "gray30", linewidth = 0.8) +

                # Visual Spike at 1.0 showing the Censored Boundary Accumulation
                ggplot2::geom_segment(
                        ggplot2::aes(x = 1, xend = 1, y = 0, yend = max_dens),
                        linetype = "dashed", color = "royalblue", linewidth = 0.6
                ) +
                ggplot2::geom_point(
                        ggplot2::aes(x = 1, y = max_dens),
                        color = "royalblue", size = 2.5
                ) +

                # Label the accumulated probability value next to the spike
                ggplot2::annotate(
                        "text",
                        x = 1.02, y = max_dens * 0.9,
                        label = paste0("Mass > 1.0:\n", round(accumulated_mass, 4)),
                        color = "royalblue", hjust = 0, size = 3
                ) +

                # Mark the Observed Score (Capped at 1.0 for viewport intercept alignment)
                ggplot2::geom_vline(
                        xintercept = min(node$obs, 1),
                        color = "firebrick", linewidth = 1, linetype = "solid"
                ) +

                # Annotate the observed value status
                ggplot2::annotate(
                        "text",
                        x = min(node$obs, 1) - 0.02, y = max_dens * 0.5,
                        label = paste("Observed Score =", round(node$obs, 3)),
                        color = "firebrick", angle = 90, vjust = 0, fontface = "bold"
                ) +

                # Titles and Scales
                ggplot2::labs(
                        title = paste("Censored Null Distribution Analysis"),
                        subtitle = paste0("Node: ", node$id, "  |  Pr(Null >= Obs) = ", round(node$perf, 4)),
                        x = "Separability Score",
                        y = "Probability Density"
                ) +
                ggplot2::scale_x_continuous(limits = c(0, 1.15), expand = c(0, 0)) +
                ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.1))) +
                ggplot2::theme_minimal() +
                ggplot2::theme(
                        plot.title = ggplot2::element_text(face = "bold", size = 12),
                        plot.subtitle = ggplot2::element_text(color = "gray40", size = 10),
                        panel.grid.minor = ggplot2::element_blank()
                )

        return(p)
}
