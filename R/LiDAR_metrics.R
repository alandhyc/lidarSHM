#' LiDAR Metrics Function
#'
#' This function is used to calculate various LiDAR metrics used in ecology. The function takes LAS or LAScatalog objects created with the lidR package as inputs and returns a multi-layered terra SpatRaster, with each band corresponding to a metric. The choice of metrics, band names, and definitions are mostly based on Shokirov et al. (2023), with some minor changes being made for metrics that were not clearly defined in the paper. Notably, the roughness_L1, roughness_L2, roughness_L3 have not been included in the current edition.
#'
#' Shokirov, S. et al. (2023) ‘Habitat highs and lows: Using terrestrial and UAV LiDAR for modelling avian species richness and abundance in a restored woodland’, Remote Sensing of Environment. Elsevier, 285, p. 113326. doi: 10.1016/J.RSE.2022.113326.
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
#' @param metrics Metrics to calculate. Defaults to "shokirov", but could be any subset of c("maxH", "meanH", "stdH", "skewH", "kurH", "p_05", "p_10", "p_25", "p_50", "p_75", "p_90", "p_95", "p_99", "VCI_2", "VCI_5", "VCI_10", "VCI_15", "VCI_20", "Cov", "height_cv", "canopy_roughness", "canopy_shannon", "Tvolume", "vlayer_L1", "vlayer_L2", "vlayer_L3", "meanH_L1", "sdH_L1", "meanH_L2", "sdH_L2", "meanH_L3", "sdH_L3", "vci_L1", "vci_L2", "vci_L3")
#' @return A SpatRaster of metrics, each band being one metric, metric name in band name
#' @import lidR
#' @export

LiDAR_metrics<-function(las,
                        res=5,
                        h_cutoff=1.3,
                        mcc_s=1.5,
                        mcc_t=0.3,
                        cov_grid=0.25,
                        zmax = 35,
                        shannon_cut=c(-1,2,5,10,15,35),
                        vox_res=0.5,
                        L1_range=c(0,1),
                        L2_range=c(1,10),
                        L3_range=c(10,35),
                        metrics = "shokirov"
                        ){

  #Define the subset of metrics we want

  metric_names<-c("maxH","meanH","stdH","skewH","kurH","p_05","p_10","p_25","p_50","p_75","p_90","p_95","p_99",
                  paste0("VCI_",c(2,5,10,15,20)),
                  "Cov",
                  "height_cv","canopy_shannon",
                  "Tvolume",
                  "vlayer_L1","vlayer_L2","vlayer_L3",
                  "meanH_L1","sdH_L1","meanH_L2","sdH_L2","meanH_L3","sdH_L3",
                  paste0("vci_L",c(1,2,3)))

  additional_metrics<-c()

  if(metrics == "shokirov"){
    metric_names<-metric_names
  } else {
    metric_names<-metrics
  }

  #Now calculate some metrics

  if(is(las,"LAS")){

    r_list<-list() #Create a list to store results

    las<-classify_ground(las,mcc(s=mcc_s,t=mcc_t))
    las<-normalize_height(las,knnidw())

    #Most metrics carries a filter of Z>1.3

    las_filtered<-filter_poi(las,Z>h_cutoff)


    #If las_filtered is empty, create a stack of empty rasters

    empty_raster<-pixel_metrics(las,~list(zmax = mean(Z,na.rm = T)),res = res)
    terra::values(empty_raster)<-0

    #Empty point cloud, create a stack from las

    if(nrow(las_filtered@data)==0){

      std_names<-c("maxH","meanH","stdH","skewH","kurH","p_05","p_10","p_25","p_50","p_75","p_90","p_95")

      if(sum(std_names %in% metric_names)>=1){
        std_metrics<-terra::rast(replicate(12,empty_raster))
        names(std_metrics)<-std_names
        wanted<-std_names[which(std_names) %in% metric_names]
        std_metrics<-std_metrics[[wanted]]
        r_list$std_metrics<-std_metrics
      }

      if("p_99" %in% metrics){
        q99<-empty_raster
        names(q99)<-"p_99"
        r_list$q99<-q99
      }

      if(sum(paste0("VCI_",c(2,5,10,15,20)) %in% metric_names) >=1){
        VCI_combined<-terra::rast(replicate(5,empty_raster))
        names(VCI_combined)<-paste0("VCI_",c(2,5,10,15,20))
        wanted<-paste0("VCI_",c(2,5,10,15,20))[which(paste0("VCI_",c(2,5,10,15,20)) %in% metric_names)]
        VCI_combined<-VCI_combined[[wanted]]
        r_list$VCI_combined<-VCI_combined
      }



    } else {

      #Point cloud not empty, calculation needed
      std_names<-c("maxH","meanH","stdH","skewH","kurH","p_05","p_10","p_25","p_50","p_75","p_90","p_95")

      if(sum(std_names %in% metric_names)>=1){
        #Calculate metrics
        std_metrics<-pixel_metrics(las_filtered,.stdmetrics_z,res = res)
        std_metrics<-std_metrics[[c("zmax","zmean","zsd","zskew","zkurt","zq5","zq10","zq25","zq50","zq75","zq90","zq95")]]

        names(std_metrics)<-std_names

        #Select the relevant metrics
        wanted<-std_names[which(std_names %in% metric_names)]
        std_metrics<-std_metrics[[wanted]]
        std_metrics<-terra::resample(std_metrics,empty_raster)
        r_list$std_metrics<-std_metrics
      }

      #Also add the 99th quantile
      if("p_99" %in% metric_names){
        q99<-pixel_metrics(las_filtered,
                           res = res,
                           func = ~q99_f(z = Z))
        q99<-terra::resample(q99,empty_raster)
        r_list$q99<-q99
      }

      #VCI

      las_filtered_vci<-las_filtered
      colnames(las_filtered_vci@data)[colnames(las@data)=="Z"]<-"z" #VCI() takes z not Z

      if("VCI_2" %in% metric_names){
        vci2exp<-substitute(~lidR::VCI(z,zmax,by),list(zmax = zmax, by = 2))
        vci2<-pixel_metrics(las_filtered_vci,eval(vci2exp),res = res)
        vci2<-terra::resample(vci2,empty_raster)
        r_list$VCI_2<-vci2
      }

      if("VCI_5" %in% metric_names){
        vci5exp<-substitute(~lidR::VCI(z,zmax,by),list(zmax = zmax, by = 5))
        vci5<-pixel_metrics(las_filtered_vci,eval(vci5exp),res = res)
        vci5<-terra::resample(vci5,empty_raster)
        r_list$VCI_5<-vci5
      }

      if("VCI_10" %in% metric_names){
        vci10exp<-substitute(~lidR::VCI(z,zmax,by),list(zmax = zmax, by = 10))
        vci10<-pixel_metrics(las_filtered_vci,eval(vci10exp),res = res)
        vci10<-terra::resample(vci10,empty_raster)
        r_list$VCI_10<-vci10
      }

      if("VCI_15" %in% metric_names){
        vci15exp<-substitute(~lidR::VCI(z,zmax,by),list(zmax = zmax, by = 15))
        vci15<-pixel_metrics(las_filtered_vci,eval(vci15exp),res = res)
        vci15<-terra::resample(vci15,empty_raster)
        r_list$VCI_15<-vci15
      }

      if("VCI_20" %in% metric_names){
        vci20exp<-substitute(~lidR::VCI(z,zmax,by),list(zmax = zmax, by = 20))
        vci20<-pixel_metrics(las_filtered_vci,eval(vci20exp),res = res)
        vci20<-terra::resample(vci20,empty_raster)
        r_list$VCI_20<-vci20
      }

    } # End of check empty else{}


    #Canopy cover (cov)
    #For each pixel, the area of 0.25m (or cov_grid) pixels being canopy (max_z>1.3) is divided by the total area of the pixel

    if("Cov" %in% metric_names){

      cov_exp<-substitute(~cov_f(Z,cov_grid),list(cov_grid = cov_grid))

      cov<-pixel_metrics(
        las,
        res = cov_grid,
        func = eval(cov_exp)
      )

      cov<-terra::resample(cov,vci2,method = "sum")
      cov<-cov/(res^2)
      names(cov)<-"Cov"

      cov<-terra::resample(cov,empty_raster)

      r_list$Cov<-cov
    }


    #Roughness metrics
    #Takes the filtered point cloud (>1.3m) as input

    rough_names<-c("height_cv","canopy_shannon")

    if(sum(rough_names %in% metric_names) >= 1){

      if(nrow(las_filtered@data)==0){

        #Empty point cloud

        rough<-terra::rast(replicate(3,empty_raster))
        names(rough)<-rough_names

        rough<-rough[[rough_names[which(rough_names %in% metric_names)]]]

        r_list$rough<-rough

      } else {

        #Something in the point cloud

        rough_exp<-substitute(~roughness_metrics_f(Z,shannon_cut),list(shannon_cut = shannon_cut))

        rough<-pixel_metrics(
          las,
          func = eval(rough_exp),
          res = res)

        rough<-terra::resample(rough,empty_raster)

        rough<-rough[[rough_names[which(rough_names %in% metric_names)]]]

        r_list$rough<-rough
      }
    }

    #Vegetation volume
    #Number of 0.5 m3 voxels divided by 8

    las_nonground<-filter_poi(las_filtered,Classification==1)

    if("Tvolume" %in% metric_names){
      if(nrow(las_nonground@data)==0){

        #Empty point cloud

        Tvolume<-empty_raster
        names(Tvolume)<-"Tvolume"

        r_list$Tvolume<-Tvolume

      } else{

        vox_expr<-substitute(~vox_f(vox_res),list(vox_res = vox_res))

        las_vox<-voxel_metrics(las_nonground,
                               func = eval(vox_expr),
                               res = vox_res)

        las_vox<-LAS(las_vox)

        Tvolume<-pixel_metrics(las_vox,~list(Tvolume = sum(vol)),res = res)

        Tvolume<-terra::resample(Tvolume,empty_raster)

        r_list$Tvolume<-Tvolume

      } #End of empty check

    } #End of metrics name check


    #Layer metrics
    #Let's first create subsets of the point clouds

    L1<-filter_poi(las_nonground,Z<=L1_range[2] & Z>=L1_range[1])
    L2<-filter_poi(las_nonground,Z<=L2_range[2] & Z>L2_range[1])
    L3<-filter_poi(las_nonground,Z<=L3_range[2] & Z>L3_range[1])

    #vlayer
    #volume in each layer

    #L1

    if("vlayer_L1" %in% metric_names){

      if(nrow(L1@data)!=0){

        vlayer_L1<-voxel_metrics(L1,
                                 eval(vox_expr),
                                 res = vox_res)
        vlayer_L1<-LAS(vlayer_L1)
        vlayer_L1<-pixel_metrics(vlayer_L1,~list(vlayer_L1=sum(vol)),res = res)
        vlayer_L1<-terra::resample(vlayer_L1,empty_raster)
        r_list$vlayer_L1<-vlayer_L1

      } else {

        vlayer_L1<-empty_raster
        names(vlayer_L1)<-"vlayer_L1"
        r_list$vlayer_L1<-vlayer_L1

      } #End of empty check

    } #End of metric names check


    #L2

    if("vlayer_L2" %in% metric_names){

      if(nrow(L2@data)!=0){
        vlayer_L2<-voxel_metrics(L2,
                                 eval(vox_expr),
                                 res = vox_res)
        vlayer_L2<-LAS(vlayer_L2)
        vlayer_L2<-pixel_metrics(vlayer_L2,~list(vlayer_L2=sum(vol)),res = res)
        vlayer_L2<-terra::resample(vlayer_L2,empty_raster)
        r_list$vlayer_L2<-vlayer_L2

      } else {
        vlayer_L2<-empty_raster
        names(vlayer_L2)<-"vlayer_L2"
        r_list$vlayer_L2<-vlayer_L2
      } #End of empty check
    } #End of metric name check


    #L3

    if("vlayer_L3" %in% metric_names){

      if(nrow(L3@data)!=0){
        vlayer_L3<-voxel_metrics(L3,
                                 eval(vox_expr),
                                 res = vox_res)
        vlayer_L3<-LAS(vlayer_L3)
        vlayer_L3<-pixel_metrics(vlayer_L3,~list(vlayer_L3=sum(vol)),res = res)
        vlayer_L3<-terra::resample(vlayer_L3,empty_raster)
        r_list$vlayer_L3<-vlayer_L3

      } else {

        vlayer_L3<-empty_raster
        names(vlayer_L3)<-"vlayer_L3"
        r_list$vlayer_L3<-vlayer_L3

      } #End of empty check

    } #End of metric name check



    #mean, sd, roughness, and vci

    meansd_names<-c("meanH_L1","sdH_L1","meanH_L2","sdH_L2","meanH_L3","sdH_L3")

    if(sum(meansd_names %in% metric_names) >= 1){

      #Define empty raster

      empty_raster2<-c(empty_raster,empty_raster)
      names(empty_raster2)<-c("meanH","sdH")

      #Function

      #L1

      if(nrow(L1@data)!=0){
        mean_sd_L1<-pixel_metrics(L1,
                                  res = res,
                                  func = ~list(meanH = mean(Z,na.rm = T),
                                               sdH = sd(Z,na.rm = T)))
      } else {
        mean_sd_L1<-empty_raster2
      }

      #L2

      if(nrow(L2@data)!=0){
        mean_sd_L2<-pixel_metrics(L2,
                                  res = res,
                                  func = ~list(meanH = mean(Z,na.rm = T),
                                               sdH = sd(Z,na.rm = T)))
      } else {
        mean_sd_L2<-empty_raster2
      }

      #L3

      if(nrow(L3@data)!=0){
        mean_sd_L3<-pixel_metrics(L3,
                                  res = res,
                                  func = ~list(meanH = mean(Z,na.rm = T),
                                               sdH = sd(Z,na.rm = T)))
      } else {
        mean_sd_L3<-empty_raster2
      }

      mean_sd_L1<-terra::resample(mean_sd_L1,empty_raster)
      mean_sd_L2<-terra::resample(mean_sd_L2,empty_raster)
      mean_sd_L3<-terra::resample(mean_sd_L3,empty_raster)

      names(mean_sd_L1)<-paste0(names(mean_sd_L1),"_L1")
      names(mean_sd_L2)<-paste0(names(mean_sd_L2),"_L2")
      names(mean_sd_L3)<-paste0(names(mean_sd_L3),"_L3")

      mean_sd_L123<-c(mean_sd_L1,mean_sd_L2,mean_sd_L3)
      mean_sd_L123<-mean_sd_L123[[meansd_names[which(meansd_names %in% metric_names)]]]

      r_list$mean_sd_L123<-mean_sd_L123

    } #End of metric name check




    #VCI by layer
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
    if("vci_L1" %in% metric_names){

      if(nrow(L1@data)!=0){
        vci_L1<-pixel_metrics(L1_VCI,func = eval(vci_L1_exp),res = res)
        vci_L1<-terra::resample(vci_L1,empty_raster)
      } else {
        vci_L1<-empty_raster
      }

      names(vci_L1)<-"vci_L1"
      r_list$vci_L1<-vci_L1

    }

    #L2
    if("vci_L2" %in% metric_names){

      if(nrow(L2@data)!=0){
        vci_L2<-pixel_metrics(L2_VCI,func = eval(vci_L2_exp),res = res)
        vci_L2<-terra::resample(vci_L2,empty_raster)
      } else {
        vci_L2<-empty_raster
      }

      names(vci_L2)<-"vci_L2"
      r_list$vci_L2<-vci_L2

    }


    #L3
    if("vci_L3" %in% metric_names){

      if(nrow(L3@data)!=0){
        vci_L3<-pixel_metrics(L3_VCI,func = eval(vci_L3_exp),res = res)
        vci_L3<-terra::resample(vci_L3,empty_raster)
      } else {
        vci_L3<-empty_raster
      }

      names(vci_L3)<-"vci_L3"
      r_list$vci_L3<-vci_L3

    }


    #Put everything back together

    final_names<-lapply(r_list,function(r){
      return(terra::names(terra::rast(r)))
      })
    final_names<-do.call(c,final_names)
    all_metrics<-terra::rast(r_list)
    names(all_metrics)<-final_names

    return(all_metrics)

    rm(las,las_filtered,empty_raster,las_filtered_vci,las_nonground,L1,L2,L3,L1_VCI,L2_VCI,L3_VCI)
    gc()

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
                             cov_grid = cov_grid,
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
