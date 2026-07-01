#' Trial data for trying BHiCect 2.0
#'
#'
#' @format ## `chr_dat_l`
#' A list of tibbles with each elements collecting HiC data for chromosome 22 at a particular resolution for the HMEC cell line from Rao et al. 2014 :
#' \describe{
#'   \item{X1}{HiC bin as integer}
#'   \item{X2}{HiC bin as integer}
#'   \item{X3}{KR normalised HiC signal}
#'   \item{weight}{Box Cox transformed HiC signal}
#'
#'   ...
#' }
"chr_dat_l"

#
#' Produce adjacency matrix from tibble
#'
#' @param g_chr1 Igraph object from which to extract appropriate adjacency matrix
#'
#' @return The square and symmetric Adjacency matrix from the input graph
#' @export
#'
legacy_get_adj_mat_fn <- function(g_chr1) {
  chr_mat <- igraph::as_adjacency_matrix(g_chr1, type = "both", attr = "weight")
  diag(chr_mat) <- 0
  if (any(Matrix::colSums(chr_mat) == 0)) {
    out <- which(Matrix::colSums(chr_mat) == 0)
    chr_mat <- chr_mat[-out, ]
    chr_mat <- chr_mat[, -out]
  }
  return(chr_mat)
}
#
#' Laplacian function
#' @description Function computing the Laplacian
#' @importFrom magrittr %>%
#' @param x Square and symmetric matrix
#'
#' @return list of two elements:
#' * Tibble with two columns, one for each eigen-vector associated with the two smallest eigen-values
#' * Vector of the two smallest eigen-values
#' @export
#'
legacy_lp_fn <- function(x) {
  Dinv <- Matrix::Diagonal(nrow(x), 1 / Matrix::rowSums(x))

  lp_chr1 <- Matrix::Diagonal(nrow(x), 1) - Dinv %*% x
  if (dim(lp_chr1)[1] > 10000) {
    temp <- RSpectra::eigs_sym(lp_chr1, k = 2, sigma = 0, which = "LM", maxitr = 10000)
    tmp_tbl <- tibble::as_tibble(temp[["vectors"]],.name_repair = make.names)
    colnames(tmp_tbl) <- c("fiedler", "zero")
    tmp_tbl <- tmp_tbl %>%
      dplyr::mutate(bins = as.integer(rownames(x)))

    return(list(vectors = tmp_tbl, values = temp[["values"]]))
  } else {
    temp <- eigen(lp_chr1)
    tmp_tbl <- tibble::as_tibble(temp[["vectors"]][, c(length(temp$values) - 1, length(temp$values))],.name_repair = make.names)
    colnames(tmp_tbl) <- c("fiedler", "zero")
    tmp_tbl <- tmp_tbl %>%
      dplyr::mutate(bins = as.integer(rownames(x)))
    return(list(vectors = tmp_tbl, values = temp[["values"]][c(length(temp$values) - 1, length(temp$values))]))
  }
}

#' Sum of square function
#' @description Utility function to compute sum of squared errors
#' @param x numeric vector
#'
#' @return sum of square error
#' @export
#'
legacy_ss <- function(x) {
  sum(scale(x, scale = FALSE)^2)
}

#' Simple partition table function
#' @description Function producing the statistics enabling bi-partition
#'
#' @param lp_res Output from the Laplacian function [lp_fn]
#' @param tmp_res Resolution at which to perform bi-partition
#' @importFrom rlang .data
#' @return Tibble with the partition statistics and ensuing clusters
#' @export
#'
legacy_simple_partition_tbl_fn <- function(lp_res, tmp_res) {
  smpl_thresh_tbl <- lp_res$vectors %>%
    dplyr::mutate(stat = purrr::map_dbl(fiedler, function(x) {
      cl_a <- fiedler[which(fiedler <= x)]
      cl_b <- fiedler[which(fiedler > x)]
      return(legacy_ss(rep(c(mean(cl_a), mean(cl_b)), c(length(cl_a), length(cl_b)))) / legacy_ss(fiedler))
    }))
  smpl_thresh <- smpl_thresh_tbl %>%
    dplyr::slice_max(stat) %>%
    dplyr::select(fiedler) %>%
    unlist()
  smpl_thresh_tbl <- smpl_thresh_tbl %>%
    dplyr::mutate(
      smpl.cl = ifelse(fiedler <= smpl_thresh, 1, 2),
      res = tmp_res
    )
  return(smpl_thresh_tbl)
}

#' Produce bin partitions from clustering step
#' @importFrom rlang !!! .data
#' @param reff_g Parent cluster igraph network object
#' @param smpl_thresh_tbl Tibble with bi-partition statistics
#' @param tmp_res Resolution at which cluster is partitioned
#' @param cl_var Variable from smpl_thresh_tbl to use for partition
#'
#' @return list where each element is the bin content of the produced partitions
#' @export
#'
legacy_partition_fn <- function(reff_g, smpl_thresh_tbl, tmp_res, cl_var) {
  sub_g_list <- list()
  ncl <- as.numeric(unlist(smpl_thresh_tbl %>%
    dplyr::select(!!cl_var) %>%
    dplyr::distinct()))
  for (j in ncl) {
    # create the subnetwork
    tmp_set <- as.integer(smpl_thresh_tbl %>%
      dplyr::filter(!!!rlang::parse_exprs(paste(cl_var, "==", j))) %>%
      dplyr::select(bins) %>%
      unlist())

    sub_g_temp <- igraph::induced_subgraph(reff_g, which(V(reff_g)$name %in% as.character(tmp_set)))

    # save the members of considered cluster
    # create cluster label in considered subnetwork
    temp_name <- paste(tmp_res, length(tmp_set),
      length(E(sub_g_temp)),
      min(tmp_set), max(tmp_set),
      sep = "_"
    )
    sub_g_list[[temp_name]] <- V(sub_g_temp)$name
    rm(sub_g_temp)
  }
  return(sub_g_list)
}

#' Select Best resoluion function
#' @description function selecting the best resolution at which to perform clustering
#' @param i Considered cluster
#' @param tmp_chr_dat_l List of tibble containing the cluster HiC at multiple resolutions
#' @param chr1_tree_cl List with the bin content of the already found clusters
#' @param res_num Named numeric vector with the HiC resolutions present in the data
#' @param tmp_res Current cluster resolution
#' @param tmp_res_set Named numeric vector containing the subset of HiC resolution applicable for the considered cluster
#' @param min_res Character indicating the minimum bin-size available for HiC data
#'
#' @return Character with the best resolution at which to perform bi-partition
#' @export
#'
legacy_select_best_res_fn <- function(i, tmp_chr_dat_l, chr1_tree_cl, res_num, tmp_res, tmp_res_set, min_res) {
  if (tmp_res == min_res) {
    res_select <- 0
  } else {
    # if cluster resolution > 5kb
    # check which resolution best for this cluster starting from the considered cluster resolution
    # vector containing the resolution selection criteria
    # only consider resolutions equal or higher than the considered cluster

    # create vector indicating which resolution adhere to the sparsity criteria
    res_select <- vector("logical", length(tmp_res_set))
    names(res_select) <- tmp_res_set
    # loop through selected resolutions
    for (r in tmp_res_set) {
      # when considered resolution equals the original resolution
      if (r == tmp_res) {
        tmp_dat <- tmp_chr_dat_l[[r]]
        nbin <- length(unique(c(tmp_dat$X1, tmp_dat$X2)))
        res_select[r] <- nrow(tmp_dat) - nbin > nbin * (nbin - 1) / 2 * 0.99
      }
      # when resolution higher than the original resolution
      else {
        # create the higher resolution bins expecting a 0-start counting of bins
        r_bin <- unique(as.integer(sapply(chr1_tree_cl[[i]], function(x) {
          tmp <- seq(as.integer(x), as.integer(x) + res_num[tmp_res], by = res_num[r])
          return(tmp[-length(tmp)])
        })))
        # extract corresponding edgelist from Hi-C data
        tmp_dat <- tmp_chr_dat_l[[r]]
        nbin <- length(r_bin)
        res_select[r] <- nrow(tmp_dat) - nbin > nbin * (nbin - 1) / 2 * 0.99
      }
    }
  }
  # if no higher or equal resolution is appropriate for further clustering keep original resolution
  if (sum(res_select) == 0) {
    tmp_res <- unlist(strsplit(i, split = "_"))[1]
  }
  # set best Hi-C resolution
  if (sum(res_select) > 0) {
    tmp_res <- names(res_select[max(which(res_select))])
  }

  # build the cluster edgelist at best resolution
  return(tmp_res)
}
# Main clustering function

#' BHiCect 2.0
#' @description Main clustering function implementing recursive spectral clustering of HiC data
#' @param res_set Character vector with the HiC resolutions present in the data
#' @param res_num Named numeric vector with the HiC resolutions present in the data
#' @param chr_dat_l List of tibble containing the chromosome-wide HiC at multiple resolutions
#' @param cl_var Variable from smpl_thresh_tbl to use for partition
#' @param nworkers Number of processes to launch for parallel computation
#' @importFrom igraph E V
#' @importFrom magrittr %>%
#' @importFrom future plan multisession sequential
#' @importFrom rlang .data
#' @return List with:
#' * Named list collecting the bin content of each found cluster
#' * Edge-list for the bi-partition tree
#' * Partition statistic table for each cluster
#' @export
#'
legacy_BHiCect <- function(res_set, res_num, chr_dat_l, cl_var = "smpl.cl", nworkers) {
  # initialisation
  # container to save cluster hierarchy as list of lists
  chr1_tree_df <- tibble::tibble(
    from = character(),
    to = character()
  )
  # containers for cluster member list
  chr1_tree_cl <- list()
  # container for cluster statistics
  chr_cl_stat_l <- list()

  min_res <- names(which.min(res_num))
  # compute the best resolution for considered chr given sparsity criteria
  tmp_res <- res_set[max(which(unlist(lapply(chr_dat_l, function(x) {
    nbin <- length(unique(c(x$X1, x$X2)))
    return(nrow(x) - nbin > nbin * (nbin - 1) / 2 * 0.99)
  }))))]
  # if even the coarcest resolution doesn't fulffil the cirteria assign the starting resolution to it
  if (is.na(tmp_res)) {
    tmp_res <- res_set[1]
  }
  g_chr1 <- igraph::graph_from_data_frame(chr_dat_l[[tmp_res]], directed = F)
  # eliminate self loop
  g_chr1 <- igraph::delete_edges(g_chr1, E(g_chr1)[which(igraph::which_loop(g_chr1))])
  chr_mat <- legacy_get_adj_mat_fn(g_chr1)
  # whole chromosome laplacian
  lpe_chr1 <- lp_fn(chr_mat)
  # spectral clusters
  smpl_thresh_tbl <- legacy_simple_partition_tbl_fn(lpe_chr1, tmp_res)
  res_chr1 <- legacy_partition_fn(g_chr1, smpl_thresh_tbl, tmp_res, cl_var)
  print(res_chr1)
  # save cluster membership and expansion
  chr1_tree_cl <- list(chr1_tree_cl, res_chr1)
  chr1_tree_cl <- unlist(chr1_tree_cl, recursive = FALSE)
  print(chr1_tree_cl)
  # save cluster statistics
  chr_cl_stat_l[["Root"]] <- smpl_thresh_tbl
  # temporary list of candidate cluster to further partition
  ok_part <- names(chr1_tree_cl)
  print(ok_part)
  # initiate the tree
  chr1_tree_df <- dplyr::bind_rows(chr1_tree_df, do.call(dplyr::bind_rows, lapply(ok_part, function(x) {
    tibble::tibble(from = "Root", to = x)
  })))

  # given starting best resolution set the further resolutions to consider at later iterations
  tmp_res_set <- res_set[which(res_set == tmp_res):length(res_set)]
  # build a temporary list containing the subset of Hi-C interactions constitutive of the found clusters at all
  # higher resolutions for the considered cluster
  cl_dat_l <- vector("list", length(ok_part))
  for (cl in ok_part) {
    for (r in tmp_res_set) {
      # when considered resolution equals the original resolution
      if (r == tmp_res) {
        tmp_dat <- chr_dat_l[[r]] %>%
          dplyr::filter(.data$X1 %in% as.integer(chr1_tree_cl[[cl]])) %>%
          dplyr::filter(.data$X2 %in% as.integer(chr1_tree_cl[[cl]]))
        cl_dat_l[[cl]][[r]] <- tmp_dat
      }
      # when resolution higher than the original resolution
      else{

        # create the higher resolution bins expecting a 0-start counting of bins
        r_bin<-unique(as.integer(sapply(chr1_tree_cl[[cl]],function(x){

          tmp <- seq(from = as.integer(x), to = (as.integer(x) + res_num[tmp_res]), by = res_num[r])
          return(tmp[-length(tmp)])
        })))

        # extract corresponding edgelist from Hi-C data
        tmp_dat <- chr_dat_l[[r]] %>%
          dplyr::filter(.data$X1 %in% as.integer(r_bin)) %>%
          dplyr::filter(.data$X2 %in% as.integer(r_bin))
        cl_dat_l[[cl]][[r]] <- tmp_dat
      }
    }
  }
  rm(lpe_chr1,res_chr1,r,r_bin,cl,tmp_res,tmp_res_set)


  # recursive looping
  while (length(ok_part) > 0) {
    message(paste("Processing", length(ok_part), "clusters"))
    # Extract for every candidate cluster nested bi-partition
    ## Get:
    ## - Parent partition stats
    ## - Children bin-content
    ## - Children-Parent edge-list
    if (length(ok_part) > 600) {
      future::plan(multisession, workers = nworkers)
    } else {
      future::plan(sequential)
    }
    ok_part_res <- furrr::future_map(ok_part, function(i) {
      tmp_res <- unlist(strsplit(i, split = "_"))[1]
      # Only evaluate higher resolution if the considered cluster was found at coarser resolution than highest resolution
      tmp_res_set <- names(sort(res_num[which(res_num <= res_num[tmp_res])], decreasing = T))
      # Skip cluster composed of loop at highest resolution
      if (length(chr1_tree_cl[[i]]) < 3 & tmp_res == min_res) {
        return(list(edge_tbl = NULL, cl_stat = NULL, cl_content = NULL))
      }

      tmp_chr_dat_l <- cl_dat_l[[i]]
      # Update resolution for current cluster bi-partition
      tmp_res <- legacy_select_best_res_fn(i, tmp_chr_dat_l, chr1_tree_cl, res_num, tmp_res, tmp_res_set, min_res)
      tmp_res_set <- names(sort(res_num[which(res_num <= res_num[tmp_res])], decreasing = T))
      tmp_chr_dat <- tmp_chr_dat_l[[tmp_res]]
      # Skip clusters where the best resolution corresponds to a loop
      if (length(unique(c(tmp_chr_dat$X1, tmp_chr_dat$X2))) < 3) {
        return(list(edge_tbl = NULL, cl_stat = NULL, cl_content = NULL))
      }

      sub_g1 <- igraph::graph_from_data_frame(tmp_chr_dat, directed = F)
      # eliminate self loop
      sub_g1 <- igraph::delete_edges(sub_g1, E(sub_g1)[which(igraph::which_loop(sub_g1))])
      # create the corresponding adjacency matrix
      sub_g1_adj <- legacy_get_adj_mat_fn(sub_g1)
      # Skip clusters where the best resolution corresponds to a loop after eliminating empty bins
      if (nrow(sub_g1_adj) < 1) {
        return(list(edge_tbl = NULL, cl_stat = NULL, cl_content = NULL))
      }

      # Find Fiedler vector
      lpe_sub_g1 <- legacy_lp_fn(sub_g1_adj)
      # Skip clusters where the best resolution corresponds to a loop after eliminating empty bins
      if (nrow(lpe_sub_g1$vectors) < 3) {
        return(list(edge_tbl = NULL, cl_stat = NULL, cl_content = NULL))
      }

      # find actual sub-structures using kmeans on fiedler vector
      subg_thresh_tbl <- legacy_simple_partition_tbl_fn(lpe_sub_g1, tmp_res)
      res_chr1 <- legacy_partition_fn(sub_g1, subg_thresh_tbl, tmp_res, cl_var)
      # Only consider future cluster partition if contain at least 2 bins
      ok_part_temp <- names(res_chr1)[which(unlist(lapply(res_chr1, length)) > 1)]
      tmp_tree_df <- do.call(dplyr::bind_rows, lapply(ok_part_temp, function(x) {
        tibble::tibble(from = i, to = x)
      }))

      return(list(
        edge_tbl = tmp_tree_df,
        cl_stat = subg_thresh_tbl,
        cl_content = res_chr1[ok_part_temp]
      ))
    })
    plan(sequential)

    # Add Parent-Children edges
    tmp_tree_df <- do.call(dplyr::bind_rows, lapply(ok_part_res, "[[", 1))
    chr1_tree_df <- chr1_tree_df %>%
      dplyr::bind_rows(tmp_tree_df)

    # Extract valid children for further clustering
    ok_part_temp <- unique(tmp_tree_df$to)

    # Add Children bin-content
    tmp_tree_cl_content <- unlist(lapply(ok_part_res, "[[", 3), recursive = F)
    chr1_tree_cl <- c(chr1_tree_cl, tmp_tree_cl_content)

    # Add Parent-partition stats
    tmp_cl_stat_l <- lapply(ok_part_res, "[[", 2)
    names(tmp_cl_stat_l) <- ok_part
    chr_cl_stat_l <- c(chr_cl_stat_l, tmp_cl_stat_l)

    # Extract candidate cluster HiC-data for next iteration
    if (length(ok_part_temp) > 600) {
      plan(multisession, workers = nworkers)
    } else {
      plan(sequential)
    }
    tmp_cl_dat_l <- furrr::future_map(ok_part_temp, function(cl) {
      tmp_res <- unlist(strsplit(cl, split = "_"))[1]
      tmp_res_set <- names(sort(res_num[which(res_num <= res_num[tmp_res])], decreasing = T))
      cl_parent <- tmp_tree_df %>%
        dplyr::filter(to == cl) %>%
        # dplyr::select(.data$from) %>%
        dplyr::select("from") %>%
        unlist()
      tmp_chr_dat_l <- cl_dat_l[[cl_parent]]
      cl_dat <- vector("list", length(tmp_res_set))
      names(cl_dat) <- tmp_res_set
      for (r in tmp_res_set) {
        if (r == tmp_res) {
          tmp_dat <- tmp_chr_dat_l[[r]] %>%
            dplyr::filter(.data$X1 %in% as.integer(chr1_tree_cl[[cl]])) %>%
            dplyr::filter(.data$X2 %in% as.integer(chr1_tree_cl[[cl]]))
          cl_dat[[r]] <- tmp_dat
        } else {
          # create the higher resolution bins expecting a 0-start counting of bins
          r_bin <- unique(as.integer(sapply(chr1_tree_cl[[cl]], function(x) {
            tmp <- seq(as.integer(x), as.integer(x) + res_num[tmp_res], by = res_num[r])
            return(tmp[-length(tmp)])
          })))
          # extract corresponding edgelist from Hi-C data
          tmp_dat <- tmp_chr_dat_l[[r]] %>%
            dplyr::filter(.data$X1 %in% as.integer(r_bin)) %>%
            dplyr::filter(.data$X2 %in% as.integer(r_bin))
          cl_dat[[r]] <- tmp_dat
        }
      }
      return(cl_dat)
    })
    plan(sequential)
    names(tmp_cl_dat_l) <- ok_part_temp

    # Update candidate clusters and their HiC-data for next iteration
    ok_part <- ok_part_temp
    cl_dat_l <- tmp_cl_dat_l
  }

  return(list(cl_member = chr1_tree_cl, part_tree = chr1_tree_df, cl_stat = chr_cl_stat_l))
}

#' Produce summary table of cluster genealogy
#'
#' @param chr_spec_res BHiCect result object
#' @param tmp_cl Target cluster name
#' @importFrom magrittr %>%
#' @return Table summarising the nested clusters within the target cluster
#' @export
#'
produce_bpt_cl_lvl_tbl<-function(chr_spec_res,tmp_cl){
  message("Produce tree representation for target cluster")
  chr_bpt<-data.tree::FromDataFrameNetwork(chr_spec_res$part_tree)
  node_ancestor<-chr_bpt$Get(function(x){x$Get('name',traversal='ancestor')})
  node_ancestor<-lapply(node_ancestor,'[',-1)
  node_lvl<-sort(chr_bpt$Get('level'))[-c(1)]
  # extract children for target cluster
  bpt_cl_set<-c(tmp_cl,names(which(purrr::map_lgl(node_ancestor,function(x){
    tmp_cl %in% x
  }))))

  cl_lvl_tbl<-tibble::tibble(node=bpt_cl_set,lvl=node_lvl[bpt_cl_set],bins=chr_spec_res$cl_member[bpt_cl_set]) %>%
    dplyr::mutate(res=stringr::str_split_fixed(.data$node,"_",2)[,1])
  return(cl_lvl_tbl)
}

#' Subset chromosome-wide HiC data to focus on target cluster
#'
#' @param tmp_cl_res_set Subset of resolutions observed in target cluster
#' @param cl_lvl_tbl Summary table of target cluster genealogy
#' @param chr_dat_l Chromosome-wide HiC data for target cluster
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @return List of HiC data with one element per resolution
#' @export
#'
produce_bpt_cl_dat_l<-function(tmp_cl_res_set,cl_lvl_tbl,chr_dat_l){
  message("Subsetting HiC data to target cluster")
  cl_dat_l<-lapply(tmp_cl_res_set,function(x){
    tmp_res_max<-cl_lvl_tbl %>%
      dplyr::filter(.data$res==x)
    bin_range<-range(as.integer(unique(unlist(tmp_res_max$bins))))
    chr_dat_l[[x]]%>%
      dplyr::filter(.data$X1 <= bin_range[2] & .data$X1 >= bin_range[1] & .data$X2 <= bin_range[2] & .data$X2 >= bin_range[1])
  })
  names(cl_dat_l)<-tmp_cl_res_set

  return(cl_dat_l)
}


#' Extract HiC data for full target cluster genealogy
#'
#' @param cl_lvl_tbl Summary table of target cluster genealogy
#' @param cl_dat_l Subset of HiC data for target cluster
#' @param res_num Named vector collecting observed resolution for chromosome-wide HiC data
#' @param tmp_cl_res_set Subset of HiC data resolution observed in target cluster
#' @param nworkers Number of workers for parallel computation
#' @importFrom future plan multisession sequential
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @return Tibble with one row per cluster and a column listing corresponding HiC data at matching resolution and interpolated high-resolution
#' @export
#'
extract_bpt_cl_dat_fn<-function(cl_lvl_tbl,cl_dat_l,res_num,tmp_cl_res_set,nworkers){
  message("Subsetting HiC data for ",nrow(cl_lvl_tbl)," children clusters")

  plan(multisession,workers=nworkers)
  cl_lvl_tbl<-cl_lvl_tbl %>%
           dplyr::mutate(HiC=furrr::future_pmap(list(.data$res,.data$bins),function(res,bins){
             tmp_dat<-cl_dat_l[[res]]%>%
               dplyr::filter(.data$X1 %in% as.integer(bins) & .data$X2 %in% as.integer(bins))

           }))
  plan(sequential)

  hi_res<-res_num[rev(tmp_cl_res_set)[1]]
  message("Interpolating HiC to ", rev(tmp_cl_res_set)[1])
  plan(multisession,workers=nworkers)

  cl_lvl_tbl<-cl_lvl_tbl %>%
    dplyr::mutate(HiC.hires=furrr::future_pmap(list(.data$res,.data$HiC),function(res,HiC){
      if(res_num[res] != hi_res){
        out_tbl<-purrr::map_dfr(1:nrow(HiC),.f = function(i){


          r_bin<-lapply(HiC[i,c(1,2)],function(x){
            tmp<-seq(x,x+res_num[res],by=hi_res)
            return(tmp[-length(tmp)])
          })
          tmp_df<-tidyr::expand_grid(ego=r_bin$X1,alter=r_bin$X2)
          tmp_df<-tmp_df%>%dplyr::mutate(raw=unlist(HiC[i,3]))
          tmp_df<-tmp_df%>%dplyr::mutate(pow=unlist(HiC[i,4]))
          tmp_df<-tmp_df%>%dplyr::mutate(res=res)
          tmp_df<-tmp_df%>%dplyr::mutate(color=unlist(HiC[i,5]))
          tmp_df<-tmp_df %>% dplyr::filter(.data$ego <= .data$alter)
          return(tmp_df)
        })
        return(out_tbl)
      } else{

        return(HiC %>%
          dplyr::rename(ego=.data$X1,alter=.data$X2,raw=.data$X3,pow=.data$weight) %>%
          dplyr::mutate(res=res) %>%
          dplyr::select(.data$ego,.data$alter,.data$raw,.data$pow,.data$res,.data$color))
      }


    }))
  plan(sequential)

  return(cl_lvl_tbl)
}

#' Produce edgelist for heatmap visualisation based on BPT segmentation
#'
#' @param cl_lvl_tbl Summary table for target cluster genealogy
#'
#' @return Edgelist with colormap value for correponding resolutions
#' @export
#'
produce_bpt_color_tbl<-function(cl_lvl_tbl){
  lvl_seq<-sort(unique(cl_lvl_tbl$lvl),decreasing = T)

  tmp_lvl_dat<-cl_lvl_tbl %>%
    dplyr::filter(.data$lvl==lvl_seq[1])
  tmp_seed<-do.call(dplyr::bind_rows,tmp_lvl_dat$HiC.hires)

  message("Joining levels from ",min(lvl_seq)," to ",max(lvl_seq))
  for(i in lvl_seq[-1]){
    tmp_up_lvl_dat_tbl<-cl_lvl_tbl %>%
      dplyr::filter(.data$lvl==i)
    tmp_up_lvl_dat<-do.call(dplyr::bind_rows,tmp_up_lvl_dat_tbl$HiC.hires)

    tmp_seed<- tmp_seed%>%
      dplyr::select(.data$ego,.data$alter,.data$color) %>%
      dplyr::full_join(tmp_up_lvl_dat %>%
                  dplyr::select(.data$ego,.data$alter,.data$color) %>%
                  dplyr::rename(color.b=.data$color)) %>%
      dplyr::mutate(color=ifelse(is.na(.data$color),.data$color.b,.data$color)) %>%
      dplyr::select(-c(.data$color.b))
  }
  return(tmp_seed)
}

#' Wrapper function to produce edgelist for heatmap visualisation with BPT segmentation
#'
#' @param tmp_cl Target cluster
#' @param chr_spec_res BHiCect result object
#' @param tmp_cl_res_set Named vector for HiC resolution observed within target cluster
#' @param res_num Named vector collecting observed resolution for chromosome-wide HiC data
#' @param chr_dat_l Chromosome-wide HiC data for target cluster
#' @param nworkers Number of workers for parallel computation
#'
#' @return 3-column edge-list tibble with edge-weight corresponding to color-map values
#' @export
#'
produce_bpt_heat_tbl<-function(tmp_cl,chr_spec_res,tmp_cl_res_set,res_num,chr_dat_l,nworkers){
  cl_lvl_tbl<-produce_bpt_cl_lvl_tbl(chr_spec_res,tmp_cl)
  cl_dat_l<-produce_bpt_cl_dat_l(tmp_cl_res_set,cl_lvl_tbl,chr_dat_l)
  cl_lvl_tbl<-extract_bpt_cl_dat_fn(cl_lvl_tbl,cl_dat_l,res_num,tmp_cl_res_set,nworkers)
  cl_col_dat<-produce_bpt_color_tbl(cl_lvl_tbl)
  return(cl_col_dat)
}
