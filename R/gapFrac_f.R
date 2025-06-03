#' LiDAR Metrics Gap Fraction Helper Function
#'
#' @return Lists to be used in pixel_metrics() and voxel_metrics()
#' @export

gapFrac_f<-function(Z,gap_thres){
  Z_gap<-length(which(Z<gap_thres))
  Z_total<-length(which(is.na(Z)==F))
  list(gapFrac = Z_gap/Z_total)
}
