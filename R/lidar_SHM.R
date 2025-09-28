#' SHM creation function
#'
#' A function to create shrub height models (SHM) from LiDAR point clouds. The function takes a LAS or LAScatalog object and builds an SHM based on the parameters specified.
#'
#' Default parameter values are from Chan et al. (2026) and represent values optimised for a airborne LiDAR scan in Scotland. Ericaceous shrubs in Scotland are short (mean height 0.36m), so some adjustments might need to be made if the study system contains taller shrubs.
#'
#'
#' @param las A LAS or LAScatalog object created by the `lidR` package
#' @param DTM_fp Character. Filepath to a standard DTM in the format "D:/path/DTM.tif". This is used to filter out trees before stretching. The quality of this DTM is not critical as the normalisation is later reversed, but a raster of reasonable resolution (e.g. 2m) is recommended.
#' @param DTM_res Numeric. Resolution of the new DTM produced for the second ground classification step. Defaults to 2.3.
#' @param max_shrub_ht Numeric. Points above this height is regarded as trees and filtered out after first normalisation. Defaults to 1.025.
#' @param stretch_factor Numeric. Determines how much the point cloud is being stretched out before ground classfication. Defaults to 22.
#' @param mcc_s Numeric. Scale parameter fed into the `mcc()` function during ground classification. Defaults to 0.5.
#' @param mcc_t Numeric. Curvature threshold fed into the `mcc()` function during ground classification. Defaults to 0.525.
#' @param shm_res Numeric. Resolution of the output SHM. Defaults to 0.2m.
#' @return A SpatRaster. The shrub height model (SHM).
#' @import raster
#' @import lidR
#' @import RMCC
#' @export

#LAScatalog function

lidar_SHM<-function(las,
                    DTM_fp,
                    DTM_res = 2.3,
                    max_shrub_ht = 1.025,
                    stretch_factor = 22,
                    mcc_s = 0.5,
                    mcc_t = 0.525,
                    shm_res = 0.2){

    require(raster)

    if(is(las,"LAS")){

      #Normalise height based on a coarse DTM that essentially gives average height of shrubs+ground

      DTM<-raster::raster(DTM_fp)

      las<-lidR::normalize_height(las,
                                  lidR::tin(),
                                  dtm = DTM,
                                  add_lasattribute = T)

      #Remove trees based on a threshold called "max_shrub_ht"

      las<-lidR::filter_poi(las,Z<max_shrub_ht)
      las<-lidR::filter_poi(las,Z>(-1)) #Added to remove artifacts

      #Add an attribute that indicates how much we want to stretch

      las<-lidR::add_lasattribute(las,
                                  x = las$Z*(stretch_factor-1),
                                  name = "Z_stretch",
                                  desc = "stretch_coefficient")


      #Remove the original normalisation

      las<-lidR::unnormalize_height(las)

      #Stretch Z to get a pseudoforest

      las$Z<-las$Z+las$Z_stretch

      #Use MCC to identify ground points in the pseudoforest

      las<-lidR::classify_ground(las,lidR::mcc(s = mcc_s,t = mcc_t))

      #Normalise based on the MCC ground classification

      new_dtm<-lidR::rasterize_terrain(las,
                                       res = DTM_res,
                                       algorithm = lidR::tin())

      las<-lidR::normalize_height(las,lidR::tin(),dtm = new_dtm)

      chm<-lidR::rasterize_canopy(las,
                                  res = shm_res,
                                  algorithm = lidR::p2r(subcircle = 0.1))

      chm<-chm/stretch_factor

      return(chm)

    } # End of is(las,"LAS")


    if(is(las,"LAScatalog")){

      options<-list(
        need_output_file=T,
        need_buffer = T
      )

      res<-lidR::catalog_map(las,
                             lidar_SHM,
                             DTM_fp = DTM_fp,
                             DTM_res = DTM_res,
                             max_shrub_ht=max_shrub_ht,
                             stretch_factor=stretch_factor,
                             mcc_s=mcc_s,
                             mcc_t=mcc_t,
                             shm_res=shm_res,
                             .options = options)

    } #End of is(las,"LAScatalog")

  }

