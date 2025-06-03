#' Function to calculate normalised rumple index (rumple_norm) based on point cloud
#'
#' @param chm A canopy height model in SpatRaster format
#' @param agg Aggregation factor when calculating the metrics. For instance, if the original resolution is 1m and we want to calculate CHM metrics at 50m grid, agg=50.
#' @param metrics Metrics to calculate.
#' @return A SpatRaster of metrics, each band being one metric, metric name in band name
#' @import terra
#' @import tools
#' @import geometry
#'
CHM_rumple_norm_point<-function(chm_vals, chm_r){

  if (length(chm_vals) < 3) return(NA)

  # Normalize by CHMmax (99.9th percentile)

  chm_max <- quantile(chm_vals, probs = 0.999, na.rm = TRUE)

  if (chm_max <= 0 || is.na(chm_max)) return(NA)

  chm_norm <- chm_r / chm_max



  # Convert to XYZ point cloud

  chm_df <- terra::as.data.frame(chm_norm, xy = TRUE, na.rm = TRUE)

  colnames(chm_df) <- c("x", "y", "z")



  if (nrow(chm_df) < 3) return(NA)



  # 3D surface area

  tin <- geometry::convhulln(chm_df, options = "FA")

  surface_area <- tin$area



  # 2D ground area from XY projection

  ground_area <- nrow(chm_df) * prod(terra::res(chm))


  # Rumple = surface / ground

  rumple_norm <- surface_area / ground_area

  return(rumple_norm)

}
