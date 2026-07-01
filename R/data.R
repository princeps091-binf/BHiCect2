#' Toy Resolution Mapping Result Object
#'
#' A sample output object containing hierarchical partitioning paths, node
#' architectures, and structural data frames for vignette demonstrations.
#'
#' @format A list containing multiple structural slots:
#' \describe{
#'   \item{nodes}{A tibble tracking multi-tier cluster metrics (id, res_level, depth, perf, obs)}
#'   \item{edges}{A tibble reporting all the edges to build the BHiCect clustering tree}
#' }
#' @source Pre-computed test run from the BHiCect2 pipeline.
"res_obj"
