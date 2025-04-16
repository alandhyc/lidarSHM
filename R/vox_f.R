#' LiDAR Metrics Voxel Helper Function
#'
#' @return Lists to be used in pixel_metrics() and voxel_metrics()
#' @export
#'
#Note voxel area function

vox_f<-function(vox_res){
  list(vol = vox_res^3)
}
