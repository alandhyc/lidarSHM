#' LiDAR Metrics Helper Functions
#'
#' @return Lists to be used in pixel_metrics() and voxel_metrics()
#' @export

#Quantile 99th function

q99_f<-function(z){
  list(p_99 = quantile(z,probs = 0.99))
}

#Vegetation cover function

cov_f<-function(z,cov_grid){
  maxZ<-max(z,na.rm = T)
  area<-ifelse(maxZ>1.3,cov_grid^2,0)
  list(Cov = area)
}

#Roughness metrics function

roughness_metrics_f<-function(z,shannon_cut){

  #Coefficient of variation
  CV <- sd(Z, na.rm = TRUE) / mean(Z, na.rm = TRUE)

  #Roughness according to Herrero-Huerta et al. (2020)
  roughness <- IQR(Z, na.rm = TRUE)^(median(Z, na.rm = TRUE))

  #Shannon according to Davison et al. (2023)
  z_cut <- cut(Z, shannon_cut)
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

#Note voxel area function

vox_f<-function(vox_res){
  list(vol = vox_res^3)
}
