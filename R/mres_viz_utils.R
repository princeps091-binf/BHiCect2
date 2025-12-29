
#' Produce cluster specific HiC data-object
#'
#' @param tmp_cl target cluster to visualise
#' @param chr_spec_res BHiCect result object for target cluster chromosome
#' @param tmp_cl_res Resolution at which target cluster was detected
#' @param tmp_cl_res_set Set of resolutions contained within target cluster
#' @param res_num Complete set of resolutions present in HiC data
#' @param chr_dat_l List of chromosome-wide HiC data where each element correspond to a particular resolution
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @return List where each element correspond to cluster HiC data at a particular resolution
#' @export
#'
produce_mres_cl_dat_l<-function(tmp_cl,chr_spec_res,tmp_cl_res,tmp_cl_res_set,res_num,chr_dat_l){
  # Extract corresponding interactions across resolutions
  tmp_cl_bin<-as.integer(chr_spec_res$cl_member[[tmp_cl]])
  ## Subset corresponding bins across resolutions
  tmp_cl_bin_mres_l<-vector('list',length(tmp_cl_res_set))
  names(tmp_cl_bin_mres_l)<-tmp_cl_res_set
  tmp_cl_bin_mres_l[[tmp_cl_res]]<-tmp_cl_bin
  message("Subsetting HiC data to target cluster")
  for(r in tmp_cl_res_set[-1]){
    tmp_cl_bin_mres_l[[r]]<-unlist(lapply(tmp_cl_bin,function(i){
      tmp<-seq(i,i+res_num[tmp_cl_res],by=res_num[r])
      return(tmp[-length(tmp)])
    }))
  }
  tmp_res_set<-names(chr_dat_l)[which(names(chr_dat_l) %in% tmp_cl_res_set)]
  cl_dat_l<-lapply(tmp_res_set,function(r){
    tmp_bin<-tmp_cl_bin_mres_l[[r]]
    return(chr_dat_l[[r]] %>%
             dplyr::filter(.data$X1 %in% tmp_bin & .data$X2 %in% tmp_bin))
  })
  names(cl_dat_l)<-tmp_res_set
  return(cl_dat_l)
}

#' Interpolate HiC data to highest resolution
#'
#' @param tmp_res_set Target cluster set of resolutions
#' @param cl_dat_l List where each element correspond to a tibble with HiC data for the target cluster at particular resolution
#' @param res_num Named vector collecting observed resolution for chromosome-wide HiC data
#' @param nworkers Number of workers used for parallelised computation
#' @importFrom future plan multisession sequential
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @return List of tibbles with corresponding HiC data converted to highest resolution
#' @export
#'
produce_mres_hires_interpolation<-function(tmp_res_set,cl_dat_l,res_num,nworkers){
  cl_hires_dat_l<-lapply(tmp_res_set,function(r){
    plan(multisession,workers=nworkers)
    hi_res<-res_num[rev(tmp_res_set)[1]]
    message("Interpolating ",r," HiC to ", rev(tmp_res_set)[1])
    tmp_dat_tbl<-cl_dat_l[[r]]
    if(res_num[r] != hi_res){
      out_tbl<-furrr::future_map_dfr(1:nrow(tmp_dat_tbl),.f = function(i){


        r_bin<-lapply(tmp_dat_tbl[i,c(1,2)],function(x){
          tmp<-seq(x,x+res_num[r],by=hi_res)
          return(tmp[-length(tmp)])
        })
        tmp_df<-tidyr::expand_grid(ego=r_bin$X1,alter=r_bin$X2)
        tmp_df<-tmp_df%>%dplyr::mutate(raw=unlist(tmp_dat_tbl[i,3]))
        tmp_df<-tmp_df%>%dplyr::mutate(pow=unlist(tmp_dat_tbl[i,4]))
        tmp_df<-tmp_df%>%dplyr::mutate(res=r)
        tmp_df<-tmp_df%>%dplyr::mutate(color=unlist(tmp_dat_tbl[i,5]))
        tmp_df<-tmp_df %>% dplyr::filter(.data$ego <= .data$alter)
        return(tmp_df)
      })
      plan(sequential)
      return(out_tbl)

    } else{
      return(tibble::as_tibble(tmp_dat_tbl) %>%
               dplyr::rename(ego=.data$X1,alter=.data$X2,raw=.data$X3,pow=.data$weight) %>%
               dplyr::mutate(res=r) %>%
               dplyr::select(.data$ego,.data$alter,.data$raw,.data$pow,.data$res,.data$color))

    }
  })
  names(cl_hires_dat_l)<-tmp_res_set
  return(cl_hires_dat_l)
}

#' Convert HiC to colormap
#'
#' @param cl_hires_dat_l List of HiC data at observed resolution for target cluster interpolated to highest resolution
#' @param tmp_res_set Set of resolution at which target cluster exists
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @return tibble with HiC data converted to colormap values
#' @export
produce_mres_color_tbl<-function(cl_hires_dat_l,tmp_res_set){
  lvl_seq<-rev(tmp_res_set)
  tmp_seed<-cl_hires_dat_l[[lvl_seq[[1]]]]
  message("Joining HiC data from",lvl_seq[1]," to ",lvl_seq[length(lvl_seq)])
  for(i in lvl_seq[-1]){
    tmp_up_lvl_dat_tbl<-cl_hires_dat_l[[i]]

    tmp_seed<- tmp_seed%>%
      dplyr::select(.data$ego,.data$alter,.data$color) %>%
      dplyr::full_join(tmp_up_lvl_dat_tbl %>%
                  dplyr::select(.data$ego,.data$alter,.data$color) %>%
                  dplyr::rename(color.b=.data$color)) %>%
      dplyr::mutate(color=ifelse(is.na(.data$color),.data$color.b,.data$color)) %>%
      dplyr::select(-c(.data$color.b))
  }
  return(tmp_seed)
}

#' Wrapper function to produce multi-resolution heatmap table
#'
#' @param tmp_cl target cluster
#' @param chr_spec_res BHiCect clustering result object
#' @param tmp_cl_res target cluster resolution
#' @param tmp_cl_res_set target cluster resolution set
#' @param res_num complete set of resolutions present in HiC data
#' @param chr_dat_l chromosome-wide HiC data corresponding to target cluster
#' @param nworkers number of workers for parallelisation
#'
#' @return 3-column edge-list tibble with edge-weight corresponding to color-map values
#' @export
produce_mres_heat_tbl<-function(tmp_cl,chr_spec_res,tmp_cl_res,tmp_cl_res_set,res_num,chr_dat_l,nworkers){
  cl_dat_l<-produce_mres_cl_dat_l(tmp_cl,chr_spec_res,tmp_cl_res,tmp_cl_res_set,res_num,chr_dat_l)
  hires_cl_dat_l<-produce_mres_hires_interpolation(tmp_cl_res_set,cl_dat_l,res_num,nworkers )
  cl_col_dat<-produce_mres_color_tbl(hires_cl_dat_l,tmp_cl_res_set)

  return(cl_col_dat)
}


