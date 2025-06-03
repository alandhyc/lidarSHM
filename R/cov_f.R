#' LiDAR Metrics Vegetation Cover Helper Function
#'
#' @return Lists to be used in pixel_metrics() and voxel_metrics()
#' @export

cov_f<-function(Z,h_cutoff){
  Z_canopy<-length(which(Z>h_cutoff))
  Z_total<-length(which(is.na(Z)==F))
  list(Cov = Z_canopy/Z_total)
}
