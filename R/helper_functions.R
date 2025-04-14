#' LiDAR Metrics Helper Function
#'
#' @return Lists to be used in pixel_metrics() and voxel_metrics()
#' @export

#Quantile 99th function

q99_f<-function(z){
  list(p_99 = quantile(z,probs = 0.99))
}

#' LiDAR Metrics Vegetation Cover Helper Function
#'
#' @return Lists to be used in pixel_metrics() and voxel_metrics()
#' @export
#'

cov_f<-function(z,cov_grid){
  maxZ<-max(z,na.rm = T)
  area<-ifelse(maxZ>1.3,cov_grid^2,0)
  list(Cov = area)
}

#' LiDAR Metrics Roughness Helper Function
#'
#' @return Lists to be used in pixel_metrics() and voxel_metrics()
#' @export

#Roughness metrics function

roughness_metrics_f<-function(z,shannon_cut){

  #Coefficient of variation
  CV <- sd(z, na.rm = TRUE) / mean(z, na.rm = TRUE)

  #Roughness according to Herrero-Huerta et al. (2020)
  roughness <- IQR(z, na.rm = TRUE)^(median(z, na.rm = TRUE))

  #Shannon according to Davison et al. (2023)
  z_cut <- cut(z, shannon_cut)
  z_cut <- na.omit(z_cut)
  prop <- as.numeric(table(z_cut)) / length(z_cut)
  prop <- prop[prop > 0]
  shannon <- (-1) * sum(prop * log(prop))

  list(
    height_cv = CV,
    canopy_roughness = roughness,
    canopy_shannon = shannon
  )
}

#' LiDAR Metrics Voxel Helper Function
#'
#' @return Lists to be used in pixel_metrics() and voxel_metrics()
#' @export
#'
#Note voxel area function

vox_f<-function(vox_res){
  list(vol = vox_res^3)
}
