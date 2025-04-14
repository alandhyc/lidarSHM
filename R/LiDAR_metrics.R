#' LiDAR Metrics Function
#'
#' This function is used to calculate various LiDAR metrics used in ecology. The function takes LAS or LAScatalog objects created with the lidR package as inputs and returns a multi-layered terra SpatRaster, with each band corresponding to a metric. The choice of metrics, band names, and definitions are mostly based on Shokirov et al. (2023), with some minor changes being made for metrics that were not clearly defined in the paper. Notably, the roughness_L1, roughness_L2, roughness_L3 have not been included in the current edition.
#'Shokirov, S. et al. (2023) ‘Habitat highs and lows: Using terrestrial and UAV LiDAR for modelling avian species richness and abundance in a restored woodland’, Remote Sensing of Environment. Elsevier, 285, p. 113326. doi: 10.1016/J.RSE.2022.113326.
#'
#' @param las A LAS or LAScatalog object
#' @param res Resolution in meters
#' @param mcc_s For MCC ground classification, see support document for details
#' @param mcc_t For MCC ground classification, see support document for details
#' @param cov_grid Resolution of grid when tallying pixels when calculating canopy cover
#' @param zmax Approximate maximum height of trees
#' @param shannon_cut Bins when calculating shannon diversity index for height distribution
#' @param vox_res Size of voxels when calculating vegetation volume
#' @param L1_range Height boundaries when defining first height layer (same for L2/L3)
#' @return A SpatRaster of metrics, each band being one metric, metric name in band name
#' @import lidR
#' @export

LiDAR_metrics<-function(las,
                        res=5,
                        h_cutoff=1.3,
                        mcc_s=1.5,
                        mcc_t=0.3
                        # cov_grid=0.25,
                        # zmax = 35,
                        # shannon_cut=c(-1,2,5,10,15,35),
                        # vox_res=0.5,
                        # L1_range=c(0,1),
                        # L2_range=c(1,10),
                        # L3_range=c(10,35)
                        ){

  if(is(las,"LAS")){

    las<-classify_ground(las,mcc(s=mcc_s,t=mcc_t))
    las<-normalize_height(las,knnidw())

    #Most metrics carries a filter of Z>1.3

    las_filtered<-filter_poi(las,Z>h_cutoff)

    #If las_filtered is empty, create a stack of empty rasters




    if(nrow(las_filtered@data)==0){

      #Empty point cloud, create a stack from las

      empty_raster<-pixel_metrics(las,~list(zmax = mean(Z,na.rm = T)),res = res)
      terra::values(empty_raster)<-0

      std_metrics<-terra::rast(replicate(12,empty_raster))
      names(std_metrics)<-c("maxH","meanH","stdH","skewH","kurH","p_05","p_10","p_25","p_50","p_75","p_90","p_95")

      q99<-empty_raster
      names(q99)<-"p_99"

      VCI_combined<-terra::rast(replicate(5,empty_raster))
      names(VCI_combined)<-paste0("VCI_",c(2,5,10,15,20))


    } else {

      std_metrics<-pixel_metrics(las_filtered,.stdmetrics_z,res = res)

      #Select the relevant metrics

      std_metrics<-std_metrics[[c("zmax","zmean","zsd","zskew","zkurt","zq5","zq10","zq25","zq50","zq75","zq90","zq95")]]

      names(std_metrics)<-c("maxH","meanH","stdH","skewH","kurH","p_05","p_10","p_25","p_50","p_75","p_90","p_95")

      #Also add the 99th quantile

      # q99_exp<-substitute(lidRmetrics::q99_f(z,zmax),list(zmax = zmax))

      q99<-pixel_metrics(las_filtered,
                         res = res,
                         func = q99_f(z = Z))

      std_metrics<-c(std_metrics,q99)

      #VCI

      # vci2exp<-substitute(lidR::VCI(z,zmax,by),list(zmax = zmax, by = 2))
      # vci5exp<-substitute(lidR::VCI(z,zmax,by),list(zmax = zmax, by = 5))
      # vci10exp<-substitute(lidR::VCI(z,zmax,by),list(zmax = zmax, by = 10))
      # vci15exp<-substitute(lidR::VCI(z,zmax,by),list(zmax = zmax, by = 15))
      # vci20exp<-substitute(lidR::VCI(z,zmax,by),list(zmax = zmax, by = 20))

      vci2<-pixel_metrics(las_filtered,VCI(z=Z,by = 2,zmax = 35),res = res) #temp
      vci5<-pixel_metrics(las_filtered,VCI(z=Z,by = 5,zmax = 35),res = res) #temp
      vci10<-pixel_metrics(las_filtered,VCI(z=Z,by = 10,zmax = 35),res = res) #temp
      vci15<-pixel_metrics(las_filtered,VCI(z=Z,by = 15,zmax = 35),res = res) #temp
      vci20<-pixel_metrics(las_filtered,VCI(z=Z,by = 20,zmax = 35),res = res) #temp

      VCI_combined<-c(vci2,vci5,vci10,vci15,vci20)
      names(VCI_combined)<-paste0("VCI_",c(2,5,10,15,20))

    } # End of check empty else{}


    #Canopy cover (cov)
    #For each pixel, the area of 0.25m (or cov_grid) pixels being canopy (max_z>1.3) is divided by the total area of the pixel

    # cov_exp<-substitute(cov_f(z,cov_grid),list(cov_grid = cov_grid))

    cov<-pixel_metrics(
      las,
      res = cov_grid,
      func = cov_f(z = Z,cov_grid = 0.25)
    )

    cov<-terra::resample(cov,vci2,method = "sum")
    cov<-cov/(res^2)
    names(cov)<-"Cov"

    #Roughness metrics
    #Takes the filtered point cloud (>1.3m) as input


    if(nrow(las_filtered@data)==0){

      rough<-terra::rast(replicate(3,empty_raster))
      names(rough)<-c("height_cv","canopy_roughness","canopy_shannon")

    } else {

      # rough_exp<-substitute(roughness_metrics_f(z,shannon_cut),list(shannon_cut = shannon_cut))

      rough<-pixel_metrics(
        las,
        func = roughness_metrics_f(z = Z, shannon_cut = c(-1,2,5,10,15,35)), #temp
        res = res)

      rough<-terra::resample(rough,q99)
    }


    #Vegetation volume
    #Number of 0.5 m3 voxels divided by 8

    las_nonground<-filter_poi(las_filtered,Classification==1)

    # vox_expr<-substitute(vox_f(vox_res),list(vox_res = vox_res))

    las_vox<-voxel_metrics(las_nonground,
                           func = vox_f(vox_res = 0.5), #temp
                           res = vox_res)

    las_vox<-LAS(las_vox)

    Tvolume<-pixel_metrics(las_vox,~list(Tvolume = sum(vol)),res = res)

    Tvolume<-terra::resample(Tvolume,q99)


    #Layer metrics
    #Let's first create subsets of the point clouds

    L1<-filter_poi(las_nonground,Z<=L1_range[2] & Z>=L1_range[1])
    L2<-filter_poi(las_nonground,Z<=L2_range[2] & Z>L2_range[1])
    L3<-filter_poi(las_nonground,Z<=L3_range[2] & Z>L3_range[1])

    empty_raster<-Tvolume
    terra::values(empty_raster)<-0

    #vlayer
    #volume in each layer

    #L1

    if(nrow(L1@data)!=0){
      vlayer_L1<-voxel_metrics(L1,
                               vox_f(vox_res = 0.5), #temp
                               res = vox_res)
      vlayer_L1<-LAS(vlayer_L1)
      vlayer_L1<-pixel_metrics(vlayer_L1,~list(vlayer_L1=sum(vol)),res = res)

    } else {
      vlayer_L1<-empty_raster
      names(vlayer_L1)<-"vlayer_L1"
    }

    #L2

    if(nrow(L2@data)!=0){
      vlayer_L2<-voxel_metrics(L2,
                               vox_f(vox_res = 0.5), #temp
                               res = vox_res)
      vlayer_L2<-LAS(vlayer_L2)
      vlayer_L2<-pixel_metrics(vlayer_L2,~list(vlayer_L2=sum(vol)),res = res)

    } else {
      vlayer_L2<-empty_raster
      names(vlayer_L2)<-"vlayer_L2"
    }

    #L3

    if(nrow(L3@data)!=0){
      vlayer_L3<-voxel_metrics(L3,
                               vox_f(vox_res = 0.5), #temp
                               res = vox_res)
      vlayer_L3<-LAS(vlayer_L3)
      vlayer_L3<-pixel_metrics(vlayer_L3,~list(vlayer_L3=sum(vol)),res = res)
    } else {
      vlayer_L3<-empty_raster
      names(vlayer_L3)<-"vlayer_L3"
    }

    #Match ext
    vlayer_L1<-terra::resample(vlayer_L1,q99)
    vlayer_L2<-terra::resample(vlayer_L2,q99)
    vlayer_L3<-terra::resample(vlayer_L3,q99)

    #mean, sd, roughness, and vci

    #Define empty raster
    terra::values(empty_raster)<-NA
    empty_raster<-c(empty_raster,empty_raster)
    names(empty_raster)<-c("meanH","sdH")

    #Function

    #L1

    if(nrow(L1@data)!=0){
      mean_sd_L1<-pixel_metrics(L1,
                                res = res,
                                func = ~list(meanH = mean(z,na.rm = T),
                                             sdH = sd(z,na.rm = T)))
    } else {
      mean_sd_L1<-empty_raster
    }

    #L2

    if(nrow(L2@data)!=0){
      mean_sd_L2<-pixel_metrics(L2,
                                res = res,
                                func = ~list(meanH = mean(z,na.rm = T),
                                             sdH = sd(z,na.rm = T)))
    } else {
      mean_sd_L2<-empty_raster
    }

    #L3

    if(nrow(L3@data)!=0){
      mean_sd_L3<-pixel_metrics(L3,
                                res = res,
                                func = ~list(meanH = mean(z,na.rm = T),
                                             sdH = sd(z,na.rm = T)))
    } else {
      mean_sd_L3<-empty_raster
    }

    mean_sd_L1<-resample(mean_sd_L1,q99)
    mean_sd_L2<-resample(mean_sd_L2,q99)
    mean_sd_L3<-resample(mean_sd_L3,q99)

    names(mean_sd_L1)<-paste0(names(mean_sd_L1),"_L1")
    names(mean_sd_L2)<-paste0(names(mean_sd_L2),"_L2")
    names(mean_sd_L3)<-paste0(names(mean_sd_L3),"_L3")



    #VCI by layer
    #The paper did not specify bin size
    #Here we do 5 bins for L1 and L2, then use the bin size of L2 for L3 (because L3 theoretically has no upper bound)

    #Define empty raster

    empty_raster<-empty_raster[[1]]

    #Define bin size

    # bin_size_L1<-(L1_range[2]-L1_range[1])/5
    # bin_size_L2<-(L2_range[2]-L2_range[1])/5

    # vci_L1_exp<-substitute(lidR::VCI(z,zmax,by),list(zmax = L1_range[2], by = bin_size_L1))
    # vci_L2_exp<-substitute(lidR::VCI(z,zmax,by),list(zmax = L2_range[2], by = bin_size_L2))
    # vci_L3_exp<-substitute(lidR::VCI(z,zmax,by),list(zmax = L3_range[2], by = bin_size_L2))

    #L1

    if(nrow(L1@data)!=0){
      vci_L1<-pixel_metrics(L1,func = VCI(z = Z,zmax = 1,by = 0.2),res = res) #temp
    } else {
      vci_L1<-empty_raster
    }

    #L2

    if(nrow(L2@data)!=0){
      vci_L2<-pixel_metrics(L2,func = VCI(z = Z,zmax = 10,by = 1.8),res = res) #temp
    } else {
      vci_L2<-empty_raster
    }

    #L3

    if(nrow(L3@data)!=0){
      vci_L3<-pixel_metrics(L3,func = VCI(z = Z,zmax = 35,by = 1.8),res = res) #temp
    } else {
      vci_L3<-empty_raster
    }

    vci_L1<-terra::resample(vci_L1,q99)
    vci_L2<-terra::resample(vci_L2,q99)
    vci_L3<-terra::resample(vci_L3,q99)

    VCI_L123<-c(vci_L1,vci_L2,vci_L3)
    names(VCI_L123)<-paste0("vci_L",c(1,2,3))

    #Put everything back together

    all_metrics<-c(std_metrics,VCI_combined,cov,rough,Tvolume,vlayer_L1,vlayer_L2,vlayer_L3,mean_sd_L1[[1]],mean_sd_L2[[1]],mean_sd_L3[[1]],mean_sd_L1[[2]],mean_sd_L2[[2]],mean_sd_L3[[2]],VCI_L123)

    return(all_metrics)

  } # End of is(las,"LAS")

  if(is(las,"LAScatalog")){

    options<-list(
      need_output_file=T
    )

    mtrc_result<-catalog_map(las,
                             LiDAR_metrics,
                             res = res,
                             h_cutoff = h_cutoff,
                             mcc_s = mcc_s,
                             mcc_t = mcc_t,
                             # cov_grid = cov_grid,
                             # zmax = zmax,
                             # shannon_cut = shannon_cut,
                             # vox_res = vox_res,
                             # L1_range = L1_range,
                             # L2_range = L2_range,
                             # L3_range = L3_range,
                             .options = options)

  } #End of is(las,"LAScatalog")

}
