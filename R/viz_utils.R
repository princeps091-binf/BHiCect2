#' get_mres_hic_dat
#'
#' @param dat_folder folder containing the HiC data at multiple resolutions
#' @param chromo chromosome name in UCSC format for the data to visualise
#' @param hub_res_set resolutions at which we wish to visualise the data
#' @description Function to input multiresolution HiC data as a list of dataframes
#' @return list of dataframes where each element is the HiC data at a particular resolution
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
get_mres_hic_dat<-function(dat_folder,chromo,hub_res_set){
  chr_dat_l<-lapply(hub_res_set,function(x)readr::read_delim(file = paste0(dat_folder,x,'/',chromo,'.txt'),delim = '\t',col_names = F))
  names(chr_dat_l)<-hub_res_set
  chr_dat_l<-lapply(chr_dat_l,function(x){
    tmp<-x%>%
      dplyr::filter(!(is.nan(.data$X3)))%>%
      dplyr::filter(.data$X1!=.data$X2)

    tmp_bin<-unique(c(x$X1,x$X2))

    tmp_self<-tibble::tibble(X1=tmp_bin,X2=tmp_bin)
    tmp_self<-tmp_self%>%
      dplyr::left_join(x)

    min_self<-min(tmp_self$X3,na.rm=T)

    tmp_self<-tmp_self %>%
      dplyr::mutate(X3=ifelse(is.na(.data$X3),min_self,.data$X3))

    return(tmp %>%
             dplyr::bind_rows(tmp_self))
  })
  chr_dat_l<-purrr::map(chr_dat_l,function(x){
    if(nrow(x)>1e5){
      preprocessParams.r <- replicate(50,caret::BoxCoxTrans(sample(x$X3,size = 1e5,replace = F),na.rm = T)$lambda)
      preprocessParams<-round(mean(preprocessParams.r),digits = 1)
    } else{
      preprocessParams <- caret::BoxCoxTrans(x$X3,na.rm = T)$lambda

    }
    x <- tibble::as_tibble(x) %>%
      dplyr::mutate(weight=((.data$X3 ^ preprocessParams) - 1) / preprocessParams)

    m.w<-min(x$weight,na.rm = T)
    x <- x %>%
      dplyr::mutate(weight= .data$weight+(1-m.w))

    return(x)
  })
  return(chr_dat_l)
}
# Produce color-scale separating each resolution into separate color-channel

#' make_mres_color_map
#'
#' @param hub_res_set set of resolutions present in considered cluster visualised
#' @param res_set set of resolutions present in HiC data
#' @param chr_dat_l list object containing a tibble for every HiC resolution data
#'
#' @return the input list of dataframes with an additional color map column
#' @export
#'
make_mres_color_map<-function(hub_res_set,res_set,chr_dat_l){
  out_chr_dat_l<-lapply(seq_along(hub_res_set),function(x){
    tmp_dat<-chr_dat_l[[hub_res_set[x]]]
    res_idx<-which(res_set == hub_res_set[x])
    toMin<-(res_idx-1)*100 +1
    toMax<-(res_idx-1)*100 +99
    tmp_dat$color<-toMin+(tmp_dat$weight-min(tmp_dat$weight))/(max(tmp_dat$weight)-min(tmp_dat$weight))*(toMax-toMin)
    return(tmp_dat)
  })
  names(out_chr_dat_l)<-hub_res_set
  return(out_chr_dat_l)
}


#' Import R objects
#'
#' @param file string containing the path to the R object to import
#'
#' @return R objects
#' @export
#'
get_obj_in_fn<-function(file){
  out_tbl<-get(base::load(file))
  tmp_obj<-names(mget(base::load(file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(out_tbl)
}


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
