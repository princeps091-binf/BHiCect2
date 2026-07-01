#' Compute Cluster Polygons/Rectangles from Tree Summary via Parallel Processing
#'
#' Evaluates non-contiguous genomic runs from hierarchical partitioning clusters
#' and flattens them into absolute spatial 2D grid coordinates in parallel using furrr.
#'
#' @param summary_tbl A tibble or data.frame containing the cluster summaries.
#'   Must have columns: \code{ids}, \code{res_level}, and \code{depth}.
#' @param current_MresFile A list or object tracking the resolution vector
#'   (e.g., \code{current_MresFile$resolutions}).
#'
#' @return A tidy \code{tibble} containing columns \code{xmin}, \code{xmax}, 
#'   \code{ymin}, \code{ymax}, and \code{depth}.
#' @import dplyr
#' @importFrom purrr map pmap_dfr
#' @importFrom furrr future_map_dfr
#' @export
compute_cluster_rectangles <- function(summary_tbl, current_MresFile) {
        extracted_resolutions <- rev(current_MresFile$resolutions)
        # Internal helper to isolate contiguous blocks along a chromosome arm
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

        # Internal helper to construct pairs of combinations across blocks
        get_cl_plot_rectangles <- function(x, res_vec) {
                current_res <- res_vec[x |> dplyr::pull(res_level)]

                tmp_contiguous_blocks <- find_genomic_runs(unlist(x |> dplyr::pull(ids)), current_res)
                if (is.null(tmp_contiguous_blocks)) {
                        return(dplyr::tibble())
                }

                rectangle_set_df <- expand.grid(
                        r1 = seq_along(tmp_contiguous_blocks),
                        r2 = seq_along(tmp_contiguous_blocks)
                ) |>
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
        # PARALLEL PROCESSING STEP via furrr
        # ==========================================================================
        # We replace purrr::map_dfr with furrr::future_map_dfr.
        # .options = furrr::furrr_options(seed = TRUE) ensures stability if there are internal seeds.
        all_rect_tbl <- furrr::future_map_dfr(
                seq_len(nrow(summary_tbl)), 
                function(idx) {
                        get_cl_plot_rectangles(summary_tbl |> dplyr::slice(idx), extracted_resolutions)
                },
                .options = furrr::furrr_options(seed = TRUE)
        ) |> 
                dplyr::arrange(depth)

        return(all_rect_tbl)
}

#' Plot Cluster Rectangles
#'
#' Generates a ggplot2 visualization of hierarchical, non-contiguous cluster blocks.
#'
#' @param rect_tbl A data frame of calculated rectangle dimensions generated via 
#'   \code{compute_cluster_rectangles()}. Alternatively, a \code{summary_tbl} can be 
#'   passed alongside \code{current_MresFile} to compute boundaries on the fly.
#' @param current_MresFile Optional list tracking the resolution vector. Only required if
#'   passing a raw summary table instead of a pre-calculated rectangle data frame.
#' @param xlim Numeric vector of length 2. The genomic window to visualize.
#' @param hic_df Optional data.frame containing contact matrix coordinates for background plotting.
#' @param alpha Numeric. Transparency value for cluster rectangles (default: 0.75).
#'
#' @return A \code{ggplot2} object.
#' @import ggplot2
#' @import dplyr
#' @export
plot_cluster_heatmap <- function(rect_tbl,
                                 current_MresFile = NULL,
                                 xlim = NULL,
                                 hic_df = NULL,
                                 alpha = 0.75) {
        
        # If user passed a raw summary table, fallback and compute it on the fly
        if (!"xmin" %in% colnames(rect_tbl)) {
                if (is.null(current_MresFile)) {
                        stop("Must provide current_MresFile if rect_tbl is a raw summary table.")
                }
                rect_tbl <- compute_cluster_rectangles(rect_tbl, current_MresFile)
        }

        # Setup standard bounding viewports if missing
        if (is.null(xlim)) {
                xlim <- c(min(rect_tbl$xmin, na.rm = TRUE), max(rect_tbl$xmax, na.rm = TRUE))
        }

        # Filter geometry instantly to minimize canvas compilation overhead
        visible_rects <- rect_tbl |>
                dplyr::filter(xmax >= xlim[1] & xmin <= xlim[2])

        # Initialize Base Canvas Layout
        p <- ggplot2::ggplot()

        # Layer the pre-calculated structural rectangles
        p <- p +
                ggplot2::geom_rect(
                        data = visible_rects,
                        ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = depth),
                        color = NA, alpha = alpha
                ) +
                ggplot2::scale_fill_viridis_c(direction = -1, option = "magma", name = "Tree Depth") +ggplot2::coord_fixed(xlim = xlim, ylim = xlim, expand = FALSE) +
                ggplot2::labs(
                        x = "Genomic Position",
                        y = "Genomic Position",
                ) +
                ggplot2::theme_minimal()

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
        # 2. MATHEMATICAL MODELING (Conditional Framework Selection)
        # ==========================================================================
        x_bar <- node$null_mu
        v     <- (node$null_sd)^2
        
        # Test if Method of Moments is mathematically valid for a Beta distribution
        use_beta <- v < (x_bar * (1 - x_bar))

        if (use_beta) {
                # ------------------------------------------------------------------
                # PRIMARY PATH: Beta Distribution (Bounded)
                # ------------------------------------------------------------------
                model_type <- "Beta (Method of Moments)"
                
                # Analytical derivation of Alpha and Beta parameters
                moment_term <- (x_bar * (1 - x_bar) / v) - 1
                alpha_est   <- x_bar * moment_term
                beta_est    <- (1 - x_bar) * moment_term

                # Generate domain (avoiding absolute 0/1 to protect infinite limits)
                x_range <- seq(1e-4, 1 - 1e-4, length.out = 200)
                density_values <- stats::dbeta(x_range, shape1 = alpha_est, shape2 = beta_est)

                # Upper-tail mass accumulation check (Score >= 0.95)
                accumulated_mass <- stats::pbeta(0.95, shape1 = alpha_est, shape2 = beta_est, lower.tail = FALSE)
                
                # Plotting annotations specific to Beta
                meta_label <- paste0(
                        "Beta Parameters:\n",
                        "\u03b1 = ", round(alpha_est, 2), "\n",
                        "\u03b2 = ", round(beta_est, 2), "\n\n",
                        "Mass \u2265 0.95:\n", round(accumulated_mass, 4)
                )
                marker_x <- 0.95
                
        } else {
                # ------------------------------------------------------------------
                # FALLBACK PATH: Student-t Distribution (Unbounded/Heavy-Noise)
                # ------------------------------------------------------------------
                model_type <- "Student-t Fallback (High-Variance Exception)"
                
                # Unbounded domain mapping matching your original design layout
                x_range <- seq(0, 1, length.out = 200)
                density_values <- (1 / node$null_sd) * stats::dt((x_range - node$null_mu) / node$null_sd, df = node$df)

                # Calculate the original "Accumulated Tail" mass (Area from 1.0 to Infinity)
                t_stat_at_1 <- (1 - node$null_mu) / node$null_sd
                accumulated_mass <- stats::pt(t_stat_at_1, df = node$df, lower.tail = FALSE)
                
                # Plotting annotations specific to Student-t
                meta_label <- paste0(
                        "Student-t Params:\n",
                        "df = ", round(node$df, 1), "\n\n",
                        "Mass \u2265 1.00:\n", round(accumulated_mass, 4)
                )
                marker_x <- 1.00
        }

        # Build combined plotting layout
        plt_df <- tibble::tibble(score = x_range, dens = density_values)
        max_dens <- max(density_values, na.rm = TRUE)

        # ==========================================================================
        # 3. GRAPH PLOTTING
        # ==========================================================================
        p <- ggplot2::ggplot(plt_df, ggplot2::aes(x = score, y = dens))

        # Optional: Shade the empirical p-value tail area (Observed -> Maximum Range)
        if (fill_tail) {
                obs_capped <- min(pmax(node$obs, min(x_range)), max(x_range))
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
                # Continuous Model Null Distribution Curve
                ggplot2::geom_line(color = "royalblue", linewidth = 1.0) +

                # Dynamic Visual Segment indicating high-end threshold accumulation
                ggplot2::geom_segment(
                        ggplot2::aes(x = marker_x, xend = marker_x, y = 0, yend = max_dens * 0.8),
                        linetype = "dashed", color = "gray50", linewidth = 0.6
                ) +
                
                # Render the dynamically prepared parameter string block
                ggplot2::annotate(
                        "text",
                        x = marker_x + 0.02, y = max_dens * 0.8,
                        label = meta_label,
                        color = "gray30", hjust = 0, size = 3
                ) +

                # Mark the Observed Score (Capped safely inside the selected framework viewport)
                ggplot2::geom_vline(
                        xintercept = min(max(node$obs, 0), max(x_range)),
                        color = "firebrick", linewidth = 1, linetype = "solid"
                ) +

                # Annotate the observed value status
                ggplot2::annotate(
                        "text",
                        x = min(max(node$obs, 0), max(x_range)) - 0.02, y = max_dens * 0.5,
                        label = paste("Observed Score =", round(node$obs, 3)),
                        color = "firebrick", angle = 90, vjust = 0, fontface = "bold"
                ) +

                # Titles and Scales adjusted dynamically based on active modeling choice
                ggplot2::labs(
                        title = paste("Null Distribution Analysis:", model_type),
                        subtitle = paste0("Node: ", node$id, "  |  Pr(Null >= Obs) = ", round(node$perf, 4)),
                        x = "Separability Score",
                        y = "Probability Density"
                ) +
                ggplot2::scale_x_continuous(limits = c(0, 1.2), breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
                ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.1))) +
                ggplot2::theme_minimal() +
                ggplot2::theme(
                        plot.title = ggplot2::element_text(face = "bold", size = 11),
                        plot.subtitle = ggplot2::element_text(color = "gray40", size = 9),
                        panel.grid.minor = ggplot2::element_blank()
                )

        return(p)
}
