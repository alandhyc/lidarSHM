#' LiDAR Metrics Roughness Helper Function
#'
#' @return Lists to be used in pixel_metrics() and voxel_metrics()
#' @export

#Roughness metrics function

roughness_metrics_f<-function(Z,shannon_cut){

  #Coefficient of variation
  CV <- sd(Z, na.rm = TRUE) / mean(Z, na.rm = TRUE)

  #Shannon according to Davison et al. (2023)
  z_cut <- cut(Z, shannon_cut)
  z_cut <- na.omit(z_cut)
  prop <- as.numeric(table(z_cut)) / length(z_cut)
  prop <- prop[prop > 0]
  shannon <- (-1) * sum(prop * log(prop))

  list(
    height_cv = CV,
    canopy_shannon = shannon
  )
}
