#' LiDAR Metrics Vegetation Cover Helper Function
#'
#' @return Lists to be used in pixel_metrics() and voxel_metrics()
#' @export

cov_f<-function(Z,cov_grid){
  maxZ<-max(Z,na.rm = T)
  area<-ifelse(maxZ>1.3,cov_grid^2,0)
  list(Cov = area)
}
