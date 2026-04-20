#' LiDAR Metrics Gap Fraction Helper Function
#'
#' @return Lists to be used in pixel_metrics() and voxel_metrics()
#' @export

gapFrac_f<-function(Z,ReturnNumber,NumberOfReturns,ScanAngleRank,gap_thres = 2){
  
  #First get the variables related to gap fraction

  weighted_GR<-data.frame(Z = Z,nrtn = NumberOfReturns)
  weighted_GR<-weighted_GR[which(Z<gap_thres),]

  weighted_GR$nrtn = 1/weighted_GR$nrtn

  weighted_GR<-sum(weighted_GR$nrtn) #Sum of weighted ground returns

  ttl_rtn<-length(which(ReturnNumber==1)) #Total number of first returns

  #Then get average scanangle
  #Scan angles ranged from around -30 to 30
  #took the absolute to prevent negative/positive angles from cancelling out each other
  #Slope could be a problem but here we don't have info to tackle it

  avg_angle<-mean(abs(ScanAngleRank))
  avg_angle<-avg_angle*pi/180 #Degrees to radians

  #Now get angle corrected gap fraction

  gap_frac<-(weighted_GR/ttl_rtn)^(cos(avg_angle))

  list(gapFrac = gap_frac)

}