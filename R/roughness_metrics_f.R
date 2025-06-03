#' LiDAR Metrics Roughness Helper Function
#'
#' @return Lists to be used in pixel_metrics() and voxel_metrics()
#' @export

#Roughness metrics function

roughness_metrics_f<-function(Z,shannon_cut){

  #Coefficient of variation
  CV <- sd(Z, na.rm = TRUE) / mean(Z, na.rm = TRUE)

  #Robust coefficient of variation

  RCV <- IQR(Z,na.rm = T) / median(Z, na.rm = T)

  #Root mean square of LiDAR

  tmean<-mean(Z,na.rm = T)
  RMS<-sqrt(mean((Z-tmean)^2, na.rm = T))

  #Shannon according to Davison et al. (2023)
  z_cut <- cut(Z, shannon_cut)
  z_cut <- na.omit(z_cut)
  prop <- as.numeric(table(z_cut)) / length(z_cut)
  prop <- prop[prop > 0]
  shannon <- (-1) * sum(prop * log(prop))

  list(
    lidar_height_cv = CV,
    lidar_rcv = RCV,
    lidar_rms = RMS,
    lidar_canopy_shannon = shannon
  )
}
