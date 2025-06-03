#' CHM Metrics Function
#'
#' This function is to calculate various CHM-based canopy metrics in remote sensing.
#'
#' We referenced Zhang et al. (2024) when choosing metrics, but also added in some other useful metrics such as RCV.
#' The normalised rumple index is often the slowest step, so we offer a fast alternative where surface area is calculated from the CHM directly rather than converting the CHM to a triangulated point cloud.
#'
#' @param chm A canopy height model in SpatRaster format
#' @param agg Aggregation factor when calculating the metrics. For instance, if the original resolution is 1m and we want to calculate CHM metrics at 50m grid, agg=50.
#' @param gap_thres This is the cutoff height when calculating gap fraction.
#' @param metrics Metrics to calculate, defaults to "all" to calculate all metrics. Alternatively metrics could be specified by names: c("chm_maxH","chm_sdH","chm_height_cv","chm_rcv","chm_rms","chm_skewH","chm_kurH","chm_gapFrac","chm_grndFrac","chm_rumple_norm")
#' @param rumple_fast Logical. Whether to use a fast approximation of CHM surface area or triangulating point clouds.
#' @return A SpatRaster of metrics, each band being one metric, metric name in band name
#' @import terra
#' @import moments
#' @import raster
#' @import spex
#' @import sf
#' @import fasterize
#' @import tools
#' @import geometry
#' @export

chm<-terra::rast("F:/Global_wind_project/Wind_vs_aridity/GCA_chms_all_4326/20090708_chm_lspikefree_842.tif")

agg<-50
gap_thres = 2
metrics = "all"

CHM_metrics<-function(chm,agg,gap_thres = 2,metrics = "all",rumple_fast = T){

  # First choose the metrics we want to calculate

  chosen_metrics<-c("chm_maxH",
                    "chm_sdH","chm_height_cv","chm_rcv","chm_rms","chm_skewH","chm_kurH",
                    "chm_gapFrac","chm_grndFrac",
                    "chm_rumple_norm")

  if(metrics != "all"){
    chosen_metrics<-chosen_metrics[which(chosen_metrics %in% metrics)]
  }

  if(length(chosen_metrics)==0){
    cat("Error: no valid metric names supplied")
    return(NULL)
  }


  #Aggregate the raster to create a grid at the desired resolution

  chm_agg<-terra::aggregate(chm,
                            fact = agg,
                            fun = "mean",
                            na.rm = T)

  chm_agg_copy<-chm_agg

  #Turn this into an sf object

  chm_agg[,]<-1:terra::ncell(chm_agg)
  chm_agg<-terra::as.polygons(chm_agg)
  chm_agg<-sf::st_as_sf(chm_agg)

  #Now loop through each polygon (pixel) one by one

  metrics_sf<-lapply(1:nrow(chm_agg),function(row){

    #Crop out the tile and extract values

    cell<-chm_agg[row,]
    tchm<-terra::crop(chm,terra::ext(cell))
    chm_values<-terra::values(tchm,mat = F)
    chm_values<-chm_values[is.na(chm_values)==F]


    #1. Height metrics

    #Calculate chm_maxH if requested

    if("chm_maxH" %in% chosen_metrics){

      if (length(chm_values) == 0){

        cell$chm_maxH<-NA

      } else {

        cell$chm_maxH<-quantile(chm_values,
                                probs = 0.999,
                                na.rm = TRUE)
        }


    }



    #2. Height distribution metrics

    #Calculate sd if requested

    if("chm_sdH" %in% chosen_metrics){

      if (length(chm_values) == 0){

        cell$chm_sdH<-NA

      } else {

        cell$chm_sdH<-sd(chm_values,
                         na.rm = TRUE)

      }

    }

    if("chm_height_cv" %in% chosen_metrics){

      if (length(chm_values) == 0){

        cell$chm_height_cv<-NA

      } else {

        cell$chm_height_cv<-sd(chm_values,na.rm = TRUE)/mean(chm_values, na.rm = T)

      }

    }

    #Calculate RCV if requested

    if("chm_rcv" %in% chosen_metrics){

      if (length(chm_values) == 0){

        cell$chm_rcv<-NA

      } else{

        cell$chm_rcv<-IQR(chm_values,na.rm = TRUE)/median(chm_values, na.rm = T)

      }
    }

    #Calculate RMS if requested

    if("chm_rms" %in% chosen_metrics){

      if (length(chm_values) == 0){

        cell$chm_rms<-NA

      } else {

        tmean<-mean(chm_values,na.rm = T)
        trms<-sqrt(mean((chm_values-tmean)^2))
        cell$chm_rms<-trms

        }
    }

    #Calculate skew if requested

    if("chm_skewH" %in% chosen_metrics){

      if (length(chm_values) == 0){

        cell$chm_skewH<-NA

      } else {

        cell$chm_skewH<-moments::skewness(chm_values)

      }
    }

    #Calculate kurtosis if requested

    if("chm_kurH" %in% chosen_metrics){

      if (length(chm_values) == 0){

        cell$chm_kurH<-NA

      } else{

        cell$chm_kurH<-moments::kurtosis(chm_values)

      }
    }


    #3. Canopy openness metrics

    #Calculate gap fraction if requested

    if("chm_gapFrac" %in% chosen_metrics){

      if (length(chm_values) == 0){

        cell$chm_gapFrac<-NA

      } else {

        cell$chm_gapFrac<-
          sum(chm_values < gap_thres, na.rm = TRUE) / length(chm_values)

        }
    }

    #Calculate ground fraction if requested

    if("chm_grndFrac" %in% chosen_metrics){

      if (length(chm_values) == 0){

        cell$chm_grndFrac<-NA

      } else {

        cell$chm_grndFrac<-
          sum(chm_values < 0.5, na.rm = TRUE) / length(chm_values)

        }
    }


    #4. Vegetation complexity metrics

    #Calculate rumple norm if requested

    if("chm_rumple_norm" %in% chosen_metrics){

      if(rumple_fast==T){

        cell$chm_rumple_norm<-CHM_rumple_norm(
          chm_vals = chm_values,
          chm_r = tchm
        )

      } else {

        cell$chm_rumple_norm<-CHM_rumple_norm_point(
          chm_vals = chm_values,
          chm_r = tchm
        )

      }

    }

    #Now return the cell

    return(cell)

    }) #End of lapply loop to add metrics to sf


  #Bind the list of sf back to a single sf

  metrics_sf<-do.call(rbind,metrics_sf)

  #Translate polygons back into raster using fasterize()


  chosen_metrics<-c("chm_maxH",
                    "chm_sdH","chm_height_cv","chm_rcv","chm_rms","chm_skewH","chm_kurH",
                    "chm_gapFrac","chm_grndFrac",
                    "chm_rumple_norm")

  final_r<-lapply(chosen_metrics,function(field){

    fast_r<-terra::rasterize(
      terra::vect(metrics_sf),
      chm_agg_copy,
      field = field,
      fun = "mean"
    )

    return(fast_r)
  })

  final_r<-terra::rast(final_r)


  #Rename raster bands

  names(final_r)<-chosen_metrics

  #Return the metrics raster

  return(final_r)

  #Cleanup

  rm(final_r,metrics_sf,chm_agg_copy,chm_agg,chm_agg_copy,chm)
  terra::tmpFiles(current=TRUE, orphan=FALSE, old=FALSE, remove=T)

}
