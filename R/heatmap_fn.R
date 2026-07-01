#' Build heatmap matrix
#'
#' This function creates a heatmap matrix from cluster data, which can be used for visualizing relationships or interactions between different elements in a dataset.
#'
#' @param cl_mat A data frame containing an edge-list of clusters. This means that each row should represent a connection between two elements (called `ego` and `alter`), along with their associated values.
#' @param res A numeric value representing the resolution of the heatmap. This dictates the size of the bins used for grouping the data.
#' @param val_col A string indicating the name of the column in `cl_mat` that contains the values to be used for coloring the heatmap.
#'
#' @importFrom rlang .data !!
#'
#' @return A sparse matrix where each entry corresponds to a cluster relationship, designed for heatmap visualization. The matrix will have cluster values filled in based on the provided edge-list.
#' @export
#'
full_f_mat<-function(cl_mat,res,val_col){

  range_5kb<-range(unique(c(cl_mat$ego,cl_mat$alter)))
  bin_5kb<-seq(range_5kb[1],range_5kb[2],by=res)
  #add the bins not present in original Hi-C dataset
  #miss_bin<-bin_5kb[which(!(bin_5kb %in% unique(c(mat_df$X1,mat_df$X2))))]

  id_conv<-seq_along(bin_5kb)
  names(id_conv)<-bin_5kb

  cl_mat$ego_id<-id_conv[as.character(cl_mat$ego)]
  cl_mat$alter_id<-id_conv[as.character(cl_mat$alter)]

  chr_mat<-Matrix::sparseMatrix(i=cl_mat$ego_id,
                                cl_mat$alter_id,
                                x=cl_mat[[val_col]],
                                symmetric = T,
                                dimnames = list(names(id_conv),names(id_conv)))
  return(chr_mat)
}

#' Utility function to produce color maps
#'
#' @param res_set Set of resolution over which to map colors
#' @param palette Palette to use for color map
#'
#' @return list with color map hexdecimal values and corresponding breaks
#' @export
#'
custom_color_map_fn<-function(res_set,palette){
  p_breaks<-c(unlist(lapply(seq_along(res_set),function(x){
    seq((x-1)*100,(x-1)*100 +100,length.out = 101)[-101]
  })),(length(res_set)-1)*100 +100)

  p_color<-rev(RColorBrewer::brewer.pal(n=length(res_set),name = palette))
  p_col<-unlist(lapply(seq_along(res_set),function(x){
    grDevices::colorRampPalette(c("black",p_color[x]))(100)
  }))
  return(list(p.col=p_col,p.breaks=p_breaks))
}

#' Function to plot chromosome-wide heatmaps
#'
#' @param chr_dat_l List with each element containing as tibble for Hi-C data of corresponding resolution
#' @param chromo Character indicating the chromosome symbol
#' @param chr_res Character indicating the resolution at which to plot the chromsome heatmap
#' @param res_num Named integer vector mapping Resolution name to their corresponding bin size in basepair unit
#' @param color_scale Column name for the HiC values to use in the heatmap
#' @param mode Transformation to apply to the color-map of the heatmap, choose between either linear of quantile
#'
#' @return Figure containing the chromosome-wide heatmap at the indicated resolution
#' @export
#'
chromosome_heatmap<-function(chr_dat_l,chromo,chr_res,res_num,color_scale,mode='linear'){

  tmp_chr_dat<-chr_dat_l[[chr_res]]

  tmp_chr_tbl<-tmp_chr_dat %>%
    dplyr::select(.data$X1,.data$X2,!!color_scale)
  colnames(tmp_chr_tbl)<-c('ego','alter','value')

  tmp_mat<-full_f_mat(tmp_chr_tbl,res_num[chr_res],val_col = "value")
  mat_bin<-sort(as.numeric(colnames(tmp_mat)))
  tick_label<-paste0(mat_bin[range(which(mat_bin%%1e5 == 0))]/1e6,'Mb')
  tick_pos<-range(which(mat_bin%%1e5 == 0))/length(mat_bin)

  def.par <- graphics::par(no.readonly = TRUE)
  graphics::par(bg=NA)
  graphics::layout(matrix(1:2,ncol=2), width = c(2,0.75),height = c(2,2))
  if(!(mode %in% c('linear','quantile'))){
    rlang::abort("You can choose between either 'linear' or 'quantile' as values for the mode parameter, displaying linear mode")

  }
  if(mode == 'linear'){
    graphics::image(as.matrix(tmp_mat),
                    col=viridisLite::viridis(100),
                    axes = FALSE,useRaster=T,ylab=chromo)

  }
  if(mode == 'quantile'){
    chr_col_specs<-list(p.col=viridisLite::viridis(100),p.break=c(0,stats::quantile(tmp_chr_dat[[color_scale]],seq(0,1,length.out=100))))

    graphics::image(as.matrix(tmp_mat),
                    col=chr_col_specs$p.col,
                    breaks=chr_col_specs$p.break,
                    axes = FALSE,useRaster=T,ylab=chromo)

  }

  graphics::axis(2, at=tick_pos, labels=tick_label, lwd=1, pos=-1e-2)
  col_scale<- matrix(rep(rev(viridisLite::viridis(100)),50),ncol=50)

  legend_image <- grDevices::as.raster(col_scale)
  plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '')
  graphics::rasterImage(legend_image, 0, 0, 1,1)
  graphics::par(def.par)
}

#' Base R cluster heatmap
#'
#' @param cl_col_dat2 Tibble collecting all heatmap entries to plot with corresponding coordinates
#' @param col_specs List collecting the necessary components to build color-map
#' @param res_num Named vector with the resolutions and their labels present in the HiC data
#' @param chromo String indicating the chromosome containing the considered cluster
#'
#' @return Heatmap plot produced by the plotgardener backend
#' @export
#'
base_heatmap<-function(cl_col_dat2,col_specs,res_num,chromo){
  def.par <- graphics::par(no.readonly = TRUE)
  min_res<-min(res_num)
  tmp_mat<-full_f_mat(cl_col_dat2,min_res,val_col = "color")

  mat_bin<-sort(as.numeric(colnames(tmp_mat)))
  tick_label<-paste0(mat_bin[range(which(mat_bin%%1e5 == 0))]/1e6,'Mb')
  tick_pos<-range(which(mat_bin%%1e5 == 0))/length(mat_bin)
  graphics::par(bg=NA)
  graphics::layout(matrix(1:2,ncol=2), width = c(2,0.75),height = c(2,2))
  graphics::image(as.matrix(tmp_mat),
        col=col_specs$p.col,
        breaks=col_specs$p.breaks,
        axes = FALSE,useRaster=T,ylab=chromo)
  graphics::axis(2, at=tick_pos, labels=tick_label, lwd=1, pos=-1e-2)

  col_scale<- matrix(rep(rev(col_specs$p.col),50),ncol=50)

  legend_image <- grDevices::as.raster(col_scale)
  plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '')
  graphics::text(x=1.5, y = seq(1/6,1,by=1/6)-1/12, srt = -90, labels = names(res_num))
  graphics::rasterImage(legend_image, 0, 0, 1,1)
  graphics::par(def.par)
  }