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

