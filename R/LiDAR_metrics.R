#' LiDAR Metrics Function
#'
#' This function is used to calculate various LiDAR metrics used in ecology. The function takes LAS or LAScatalog objects created with the lidR package as inputs and returns a multi-layered terra SpatRaster, with each band corresponding to a metric.
#'
#' @param las A LAS or LAScatalog object
#' @param ground_classified True if the LAS object supplied is already ground classified
#' @param h_cutoff Ground and shrub points are removed when calculating metrics such as SD, skew and kurtosis. This represents the cutoff height, 1.3m by default (see Shokirov et al. (2023)).
#' @param gap_thres Threshold for calculating gap fraction. Default at 2m.
#' @param res Resolution in meters
#' @param mcc_s For MCC ground classification, see support document for details
#' @param mcc_t For MCC ground classification, see support document for details
#' @param cov_grid Resolution of grid when tallying pixels when calculating canopy cover
#' @param zmax Approximate maximum height of trees
#' @param shannon_cut Bins when calculating shannon diversity index for height distribution
#' @param vox_res Size of voxels when calculating vegetation volume
#' @param L1_range Height boundaries when defining first height layer (same for L2/L3)
#' @param metrics Metrics to calculate. Defaults to "all", but could be any subset of c("lidar_maxH", "lidar_meanH", "lidar_stdH", "lidar_skewH", "lidar_kurH", "lidar_p_05", "lidar_p_10", "lidar_p_25", "lidar_p_50", "lidar_p_75", "lidar_p_90", "lidar_p_95", "lidar_p_999", "lidar_vci_2m", "lidar_vci_5m", "lidar_vci_10m", "lidar_vci_15m", "lidar_vci_20m", "lidar_Cov","lidar_gapFrac", "lidar_grndFrac", "lidar_height_cv","lidar_rcv","lidar_rms", "lidar_canopy_shannon", "lidar_Tvolume", "lidar_vlayer_L1", "lidar_vlayer_L2", "lidar_vlayer_L3", "lidar_meanH_L1", "lidar_sdH_L1", "lidar_meanH_L2", "lidar_sdH_L2", "lidar_meanH_L3", "lidar_sdH_L3", "lidar_vci_L1", "lidar_vci_L2", "lidar_vci_L3")
#' @return A SpatRaster of metrics, each band being one metric, metric name in band name
#' @import terra
#' @import moments
#' @import sf
#' @import tools
#' @import geometry
#' @import lidR
#' @import RMCC
#' @export

# las<-lidR::readLAS("F:/Cairngorms_project/LiDAR/Temp/280000_805000.laz")
#
# ground_classified = F
# res=5
# h_cutoff=1.3
# gap_thres = 2
# mcc_s=1.5
# mcc_t=0.3
# zmax = 35
# shannon_cut=c(-1,2,5,10,15,35)
# vox_res=0.5
# L1_range=c(0,1)
# L2_range=c(1,10)
# L3_range=c(10,35)
# metrics = "all"

LiDAR_metrics<-function(las,
                        ground_classified = F,
                        res=5,
                        h_cutoff=1.3,
                        gap_thres = 2,
                        mcc_s=1.5,
                        mcc_t=0.3,
                        zmax = 35,
                        shannon_cut=c(-1,2,5,10,15,35),
                        vox_res=0.5,
                        L1_range=c(0,1),
                        L2_range=c(1,10),
                        L3_range=c(10,35),
                        metrics = "all"
                        ){

  # 1. Metric names ####
  #Define the subset of metrics we want

  metric_names<-c("lidar_maxH","lidar_meanH","lidar_stdH","lidar_skewH","lidar_kurH","lidar_p_05","lidar_p_10","lidar_p_25","lidar_p_50","lidar_p_75","lidar_p_90","lidar_p_95","lidar_p_999",
                  paste0("lidar_vci_",c(2,5,10,15,20),"m"),
                  "lidar_Cov","lidar_gapFrac","lidar_grndFrac",
                  "lidar_height_cv","lidar_rcv","lidar_rms","lidar_canopy_shannon",
                  "lidar_Tvolume","lidar_ePAI",
                  "lidar_vlayer_L1","lidar_vlayer_L2","lidar_vlayer_L3",
                  "lidar_meanH_L1","lidar_sdH_L1","lidar_meanH_L2","lidar_sdH_L2","lidar_meanH_L3","lidar_sdH_L3",
                  paste0("lidar_vci_L",c(1,2,3)))

  additional_metrics<-c()

  if(metrics == "all"){
    metric_names<-metric_names
  } else {
    metric_names<-metrics
  }

  #Now calculate some metrics

  if(is(las,"LAS")){

    # 2. Ground classification and normalisation ####

    r_list<-list() #Create a list to store results

    if(ground_classified==F){
      las<-lidR::classify_ground(las,lidR::mcc(s=mcc_s,t=mcc_t))
    }
    las<-lidR::normalize_height(las,lidR::knnidw())

    #Most metrics carries a filter of Z>1.3

    las_filtered<-lidR::filter_poi(las,Z>h_cutoff)

    # 3. Empty las handling for std.metrics ####

    #If las_filtered is empty, create a stack of empty rasters

    empty_raster<-lidR::pixel_metrics(las,~list(zmax = mean(Z,na.rm = T)),res = res)
    terra::values(empty_raster)<-0

    #Empty point cloud, create a stack from las

    if(nrow(las_filtered@data)==0){

      std_names<-c("lidar_maxH","lidar_meanH","lidar_stdH","lidar_skewH","lidar_kurH","lidar_p_05","lidar_p_10","lidar_p_25","lidar_p_50","lidar_p_75","lidar_p_90","lidar_p_95")

      if(sum(std_names %in% metric_names)>=1){
        std_metrics<-terra::rast(replicate(12,empty_raster))
        names(std_metrics)<-std_names
        wanted<-std_names[which(std_names) %in% metric_names]
        std_metrics<-std_metrics[[wanted]]
        r_list$std_metrics<-std_metrics
      }

      if("lidar_p_999" %in% metrics){
        q999<-empty_raster
        names(q999)<-"lidar_p_999"
        r_list$q999<-q999
      }

      if(sum(paste0("lidar_vci_",c(2,5,10,15,20),"m") %in% metric_names) >=1){
        lidar_VCI_combined<-terra::rast(replicate(5,empty_raster))
        names(lidar_VCI_combined)<-paste0("lidar_vci_",c(2,5,10,15,20),"m")
        wanted<-paste0("lidar_vci_",c(2,5,10,15,20),"m")[which(paste0("lidar_vci_",c(2,5,10,15,20),"m") %in% metric_names)]
        lidar_VCI_combined<-lidar_VCI_combined[[wanted]]
        r_list$lidar_VCI_combined<-lidar_VCI_combined
      }



    } else {

      # 4. Standard metrics ####

      #Point cloud not empty, calculation needed
      std_names<-c("lidar_maxH","lidar_meanH","lidar_stdH","lidar_skewH","lidar_kurH","lidar_p_05","lidar_p_10","lidar_p_25","lidar_p_50","lidar_p_75","lidar_p_90","lidar_p_95")

      if(sum(std_names %in% metric_names)>=1){
        #Calculate metrics
        std_metrics<-lidR::pixel_metrics(las_filtered,lidR::.stdmetrics_z,res = res)
        std_metrics<-std_metrics[[c("zmax","zmean","zsd","zskew","zkurt","zq5","zq10","zq25","zq50","zq75","zq90","zq95")]]

        names(std_metrics)<-std_names

        #Select the relevant metrics
        wanted<-std_names[which(std_names %in% metric_names)]
        std_metrics<-std_metrics[[wanted]]
        std_metrics<-terra::resample(std_metrics,empty_raster)
        r_list$std_metrics<-std_metrics
      }

      #5. 99.9th quantile ####

      #Also add the 999th quantile
      if("lidar_p_999" %in% metric_names){
        q999<-lidR::pixel_metrics(las_filtered,
                                  res = res,
                                  func = ~q99_f(z = Z))
        q999<-terra::resample(q999,empty_raster)
        names(q999)<-"lidar_p_999"
        r_list$q999<-q999
      }

      #6. VCI ####

      las_filtered_vci<-las_filtered
      colnames(las_filtered_vci@data)[colnames(las@data)=="Z"]<-"z" #VCI() takes z not Z

      if("lidar_vci_2m" %in% metric_names){
        vci2exp<-substitute(~lidR::VCI(z,zmax,by),list(zmax = zmax, by = 2))
        vci2<-lidR::pixel_metrics(las_filtered_vci,eval(vci2exp),res = res)
        vci2<-terra::resample(vci2,empty_raster)
        names(vci2)<-"lidar_vci_2m"
        r_list$lidar_VCI_2<-vci2
      }

      if("lidar_vci_5m" %in% metric_names){
        vci5exp<-substitute(~lidR::VCI(z,zmax,by),list(zmax = zmax, by = 5))
        vci5<-lidR::pixel_metrics(las_filtered_vci,eval(vci5exp),res = res)
        vci5<-terra::resample(vci5,empty_raster)
        names(vci5)<-"lidar_vci_5m"
        r_list$lidar_VCI_5<-vci5
      }

      if("lidar_vci_10m" %in% metric_names){
        vci10exp<-substitute(~lidR::VCI(z,zmax,by),list(zmax = zmax, by = 10))
        vci10<-lidR::pixel_metrics(las_filtered_vci,eval(vci10exp),res = res)
        vci10<-terra::resample(vci10,empty_raster)
        names(vci10)<-"lidar_vci_10m"
        r_list$lidar_VCI_10<-vci10
      }

      if("lidar_vci_15m" %in% metric_names){
        vci15exp<-substitute(~lidR::VCI(z,zmax,by),list(zmax = zmax, by = 15))
        vci15<-lidR::pixel_metrics(las_filtered_vci,eval(vci15exp),res = res)
        vci15<-terra::resample(vci15,empty_raster)
        names(vci15)<-"lidar_vci_15m"
        r_list$lidar_VCI_15<-vci15
      }

      if("lidar_vci_20m" %in% metric_names){
        vci20exp<-substitute(~lidR::VCI(z,zmax,by),list(zmax = zmax, by = 20))
        vci20<-lidR::pixel_metrics(las_filtered_vci,eval(vci20exp),res = res)
        vci20<-terra::resample(vci20,empty_raster)
        names(vci20)<-"lidar_vci_20m"
        r_list$lidar_VCI_20<-vci20
      }

    } # End of check empty else{}


    #7. Canopy cover (cov) ####
    #For each pixel, we divide the number of first returns >1.3m with total number of first returns

    las_first<-lidR::filter_first(las)

    if("lidar_Cov" %in% metric_names){

      cov_exp<-substitute(~cov_f(Z,h_cutoff),list(h_cutoff = h_cutoff))

      cov<-lidR::pixel_metrics(
        las_first,
        res = res,
        func = eval(cov_exp)
      )

      names(cov)<-"lidar_Cov"

      cov<-terra::resample(cov,empty_raster)

      r_list$lidar_Cov<-cov
    }


    # 8. Gap fraction ####

    if("lidar_gapFrac" %in% metric_names){

      gapFrac_exp<-substitute(~gapFrac_f(Z,gap_thres),list(gap_thres = gap_thres))

      gapFrac<-lidR::pixel_metrics(
        las_first,
        res = res,
        func = eval(gapFrac_exp)
      )

      gapFrac<-terra::resample(gapFrac,empty_raster)

      names(gapFrac)<-"lidar_gapFrac"

      r_list$lidar_gapFrac<-gapFrac
    }


    # 9. Ground fraction ####

    if("lidar_grndFrac" %in% metric_names){

      grndFrac<-lidR::pixel_metrics(
        las_first,
        ~list(grndFrac=length(which(Z<0.5))/length(is.na(Z)==F)),
        res = res
        )

      grndFrac<-terra::resample(grndFrac,empty_raster)
      names(grndFrac)<-"lidar_grndFrac"

      r_list$lidar_grndFrac<-grndFrac
    }


    # 10. Roughness metrics ####
    #Takes the filtered point cloud (>1.3m) as input

    rough_names<-c("lidar_height_cv","lidar_rcv","lidar_rms","lidar_canopy_shannon")

    if(sum(rough_names %in% metric_names) >= 1){

      if(nrow(las_filtered@data)==0){

        #Empty point cloud

        rough<-terra::rast(replicate(4,empty_raster))
        names(rough)<-rough_names

        rough<-rough[[rough_names[which(rough_names %in% metric_names)]]]

        r_list$rough<-rough

      } else {

        #Something in the point cloud

        rough_exp<-substitute(~roughness_metrics_f(Z,shannon_cut),list(shannon_cut = shannon_cut))

        rough<-lidR::pixel_metrics(
          las,
          func = eval(rough_exp),
          res = res)

        rough<-terra::resample(rough,empty_raster)
        names(rough)<-rough_names

        rough<-rough[[rough_names[which(rough_names %in% metric_names)]]]

        r_list$rough<-rough
      }
    }

    # 11. Vegetation volume ####

    #Number of 0.5 m3 voxels divided by 8

    las_nonground<-lidR::filter_poi(las,Classification==1)

    if("lidar_Tvolume" %in% metric_names){
      if(nrow(las_nonground@data)==0){

        #Empty point cloud

        Tvolume<-empty_raster
        names(Tvolume)<-"lidar_Tvolume"

        r_list$Tvolume<-Tvolume

      } else{

        vox_expr<-substitute(~vox_f(vox_res),list(vox_res = vox_res))

        las_vox<-lidR::voxel_metrics(las_nonground,
                                     func = eval(vox_expr),
                                     res = vox_res)

        las_vox<-lidR::LAS(las_vox)

        Tvolume<-lidR::pixel_metrics(las_vox,
                                     ~list(Tvolume = sum(vol)),
                                     res = res)

        Tvolume<-terra::resample(Tvolume,empty_raster)
        names(Tvolume)<-"lidar_Tvolume"

        r_list$Tvolume<-Tvolume

      } #End of empty check

    } #End of metrics name check

    # 12. ePAI ####

    if("lidar_ePAI" %in% metric_names){
      if(nrow(las@data)==0){

        #Empty point cloud

        ePAI<-empty_raster
        names(ePAI)<-"lidar_ePAI"

        r_list$lidar_ePAI<-ePAI

      } else{

        ePAI_expr<-substitute(~ePAI_f(Z,ReturnNumber,NumberOfReturns,ScanAngleRank,h_cutoff),list(h_cutoff = h_cutoff))

        ePAI<-lidR::pixel_metrics(las,
                                  func = eval(ePAI_expr),
                                  res = res)


        ePAI<-terra::resample(ePAI,empty_raster)
        names(ePAI)<-"lidar_ePAI"

        r_list$lidar_ePAI<-ePAI

      } #End of empty check
    } #End of name check

    # 13. Layer metrics ####
    #Let's first create subsets of the point clouds

    L1<-lidR::filter_poi(las_nonground,Z<=L1_range[2] & Z>=L1_range[1])
    L2<-lidR::filter_poi(las_nonground,Z<=L2_range[2] & Z>L2_range[1])
    L3<-lidR::filter_poi(las_nonground,Z<=L3_range[2] & Z>L3_range[1])

    ## 13.1 Volume by layer (vlayer) ####
    #volume in each layer

    #L1

    if("lidar_vlayer_L1" %in% metric_names){

      if(nrow(L1@data)!=0){

        vlayer_L1<-lidR::voxel_metrics(L1,
                                       eval(vox_expr),
                                       res = vox_res)
        vlayer_L1<-lidR::LAS(vlayer_L1)
        vlayer_L1<-lidR::pixel_metrics(vlayer_L1,
                                       ~list(vlayer_L1=sum(vol)),
                                       res = res)
        vlayer_L1<-terra::resample(vlayer_L1,empty_raster)
        names(vlayer_L1)<-"lidar_vlayer_L1"
        r_list$vlayer_L1<-vlayer_L1

      } else {

        vlayer_L1<-empty_raster
        names(vlayer_L1)<-"lidar_vlayer_L1"
        r_list$vlayer_L1<-vlayer_L1

      } #End of empty check

    } #End of metric names check


    #L2

    if("lidar_vlayer_L2" %in% metric_names){

      if(nrow(L2@data)!=0){
        vlayer_L2<-lidR::voxel_metrics(L2,
                                       eval(vox_expr),
                                       res = vox_res)
        vlayer_L2<-lidR::LAS(vlayer_L2)
        vlayer_L2<-lidR::pixel_metrics(vlayer_L2,
                                       ~list(vlayer_L2=sum(vol)),
                                       res = res)
        vlayer_L2<-terra::resample(vlayer_L2,empty_raster)
        names(vlayer_L2)<-"lidar_vlayer_L2"
        r_list$vlayer_L2<-vlayer_L2

      } else {
        vlayer_L2<-empty_raster
        names(vlayer_L2)<-"lidar_vlayer_L2"
        r_list$vlayer_L2<-vlayer_L2
      } #End of empty check
    } #End of metric name check


    #L3

    if("lidar_vlayer_L3" %in% metric_names){

      if(nrow(L3@data)!=0){
        vlayer_L3<-lidR::voxel_metrics(L3,
                                       eval(vox_expr),
                                       res = vox_res)
        vlayer_L3<-lidR::LAS(vlayer_L3)
        vlayer_L3<-lidR::pixel_metrics(vlayer_L3,
                                       ~list(vlayer_L3=sum(vol)),
                                       res = res)
        vlayer_L3<-terra::resample(vlayer_L3,empty_raster)
        names(vlayer_L3)<-"lidar_vlayer_L3"
        r_list$vlayer_L3<-vlayer_L3

      } else {

        vlayer_L3<-empty_raster
        names(vlayer_L3)<-"lidar_vlayer_L3"
        r_list$vlayer_L3<-vlayer_L3

      } #End of empty check

    } #End of metric name check



    ## 13.2 mean, sd, roughness, and vci####

    meansd_names<-c("lidar_meanH_L1","lidar_sdH_L1","lidar_meanH_L2","lidar_sdH_L2","lidar_meanH_L3","lidar_sdH_L3")

    if(sum(meansd_names %in% metric_names) >= 1){

      #Define empty raster

      empty_raster2<-c(empty_raster,empty_raster)

      #Function

      #L1

      if(nrow(L1@data)!=0){

        mean_sd_L1<-lidR::pixel_metrics(L1,
                                        res = res,
                                        func = ~list(lidar_meanH = mean(Z,na.rm = T),
                                                     lidar_sdH = sd(Z,na.rm = T)))

      } else {

        mean_sd_L1<-empty_raster2

      }

      #L2

      if(nrow(L2@data)!=0){

        mean_sd_L2<-lidR::pixel_metrics(L2,
                                        res = res,
                                        func = ~list(lidar_meanH = mean(Z,na.rm = T),
                                                     lidar_sdH = sd(Z,na.rm = T)))

      } else {

        mean_sd_L2<-empty_raster2

      }

      #L3

      if(nrow(L3@data)!=0){
        mean_sd_L3<-lidR::pixel_metrics(L3,
                                        res = res,
                                        func = ~list(lidar_meanH = mean(Z,na.rm = T),
                                               lidar_sdH = sd(Z,na.rm = T)))

      } else {

        mean_sd_L3<-empty_raster2

      }

      mean_sd_L1<-terra::resample(mean_sd_L1,empty_raster)
      mean_sd_L2<-terra::resample(mean_sd_L2,empty_raster)
      mean_sd_L3<-terra::resample(mean_sd_L3,empty_raster)

      names(mean_sd_L1)<-c("lidar_meanH_L1","lidar_sdH_L1")
      names(mean_sd_L2)<-c("lidar_meanH_L2","lidar_sdH_L2")
      names(mean_sd_L3)<-c("lidar_meanH_L3","lidar_sdH_L3")

      mean_sd_L123<-c(mean_sd_L1,mean_sd_L2,mean_sd_L3)
      mean_sd_L123<-mean_sd_L123[[meansd_names[which(meansd_names %in% metric_names)]]]

      r_list$mean_sd_L123<-mean_sd_L123

    } #End of metric name check




    ## 13.3 VCI ####
    #The paper did not specify bin size
    #Here we do 5 bins for L1 and L2, then use the bin size of L2 for L3 (because L3 theoretically has no upper bound)

    #Define empty raster

    empty_raster<-empty_raster[[1]]

    #Define bin size

    bin_size_L1<-(L1_range[2]-L1_range[1])/5
    bin_size_L2<-(L2_range[2]-L2_range[1])/5

    vci_L1_exp<-substitute(~lidR::VCI(z,zmax,by),list(zmax = L1_range[2], by = bin_size_L1))
    vci_L2_exp<-substitute(~lidR::VCI(z,zmax,by),list(zmax = L2_range[2], by = bin_size_L2))
    vci_L3_exp<-substitute(~lidR::VCI(z,zmax,by),list(zmax = L3_range[2], by = bin_size_L2))

    #Change Z to z for VCI()

    L1_VCI<-L1
    colnames(L1_VCI@data)[colnames(L1_VCI@data)=="Z"]<-"z"
    L2_VCI<-L2
    colnames(L2_VCI@data)[colnames(L2_VCI@data)=="Z"]<-"z"
    L3_VCI<-L3
    colnames(L3_VCI@data)[colnames(L3_VCI@data)=="Z"]<-"z"

    #L1
    if("lidar_vci_L1" %in% metric_names){

      if(nrow(L1@data)!=0){
        vci_L1<-lidR::pixel_metrics(L1_VCI,func = eval(vci_L1_exp),res = res)
        vci_L1<-terra::resample(vci_L1,empty_raster)
      } else {
        vci_L1<-empty_raster
      }

      names(vci_L1)<-"lidar_vci_L1"
      r_list$lidar_vci_L1<-vci_L1

    }

    #L2
    if("lidar_vci_L2" %in% metric_names){

      if(nrow(L2@data)!=0){
        vci_L2<-lidR::pixel_metrics(L2_VCI,func = eval(vci_L2_exp),res = res)
        vci_L2<-terra::resample(vci_L2,empty_raster)
      } else {
        vci_L2<-empty_raster
      }

      names(vci_L2)<-"lidar_vci_L2"
      r_list$lidar_vci_L2<-vci_L2

    }


    #L3
    if("lidar_vci_L3" %in% metric_names){

      if(nrow(L3@data)!=0){
        vci_L3<-lidR::pixel_metrics(L3_VCI,func = eval(vci_L3_exp),res = res)
        vci_L3<-terra::resample(vci_L3,empty_raster)
      } else {
        vci_L3<-empty_raster
      }

      names(vci_L3)<-"lidar_vci_L3"
      r_list$lidar_vci_L3<-vci_L3

    }


    # 14. Final raster####

    final_names<-lapply(r_list,function(r){
      return(terra::names(terra::rast(r)))
      })
    final_names<-do.call(c,final_names)
    all_metrics<-terra::rast(r_list)
    names(all_metrics)<-final_names

    return(all_metrics)

    rm(las,las_filtered,empty_raster,las_first,las_filtered_vci,las_nonground,L1,L2,L3,L1_VCI,L2_VCI,L3_VCI)
    gc()

  } # End of is(las,"LAS")

  if(is(las,"LAScatalog")){

    options<-list(
      need_output_file=T
    )

    mtrc_result<-catalog_map(las,
                             LiDAR_metrics,
                             ground_classified = ground_classified,
                             res = res,
                             h_cutoff = h_cutoff,
                             gap_thres = gap_thres,
                             mcc_s = mcc_s,
                             mcc_t = mcc_t,
                             zmax = zmax,
                             shannon_cut = shannon_cut,
                             vox_res = vox_res,
                             L1_range = L1_range,
                             L2_range = L2_range,
                             L3_range = L3_range,
                             metrics = metrics,
                             .options = options)

  } #End of is(las,"LAScatalog")

}
