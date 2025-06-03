#' LiDAR Metrics Helper Function
#'
#' @return Lists to be used in pixel_metrics() and voxel_metrics()
#' @export

#Quantile 99th function

q99_f<-function(z){
  list(p_99 = quantile(z,probs = 0.999))
}
