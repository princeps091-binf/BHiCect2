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
