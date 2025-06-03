#' Function to calculate normalised rumple index (rumple_norm)
#'
#' @param chm A canopy height model in SpatRaster format
#' @param agg Aggregation factor when calculating the metrics. For instance, if the original resolution is 1m and we want to calculate CHM metrics at 50m grid, agg=50.
#' @param metrics Metrics to calculate.
#' @return A SpatRaster of metrics, each band being one metric, metric name in band name
#' @import terra
#' @import tools
#' @import geometry

CHM_rumple_norm<-function(chm_vals, chm_r){

  if (length(chm_vals) < 3) return(NA)

  # Normalize by CHMmax (99.9th percentile)

  chm_max <- quantile(chm_vals, probs = 0.999, na.rm = TRUE)

  if (chm_max <= 0 || is.na(chm_max)) return(NA)

  chm_norm <- chm_r / chm_max

  #Calculate a slope raster

  chm_norm_slope<-terra::terrain(chm_norm,v = "slope")*pi/180

  #Use 1/cos(x) to create an area raster assuming the resolution is 1m (it's going to cancel out anyway)

  area_r<-1/cos(chm_norm_slope)

  #Sum the area

  area<-terra::global(area_r,fun = "sum",na.rm = T)
  area<-area[1,1]

  #Number of non-NA pixels

  non_na<-terra::global(area_r,fun = "notNA",na.rm = T)
  non_na<-non_na[1,1]

  rumple_norm<-area/non_na

  return(rumple_norm)

}

