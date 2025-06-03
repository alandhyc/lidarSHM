#' LiDAR Metrics ePAI Helper Function
#'
#' @return Lists to be used in pixel_metrics() and voxel_metrics()
#' @export

ePAI_f<-function(Z,ReturnNumber,NumberOfReturns,ScanAngleRank,h_cutoff){

  #First get the variables related to gap fraction

  weighted_NR<-data.frame(Z = Z,nrtn = NumberOfReturns)
  weighted_NR<-weighted_NR[which(Z>h_cutoff),]

  weighted_NR$nrtn = 1/weighted_NR$nrtn

  weighted_NR<-sum(weighted_NR$nrtn) #Sum of weighted canopy returns

  ttl_rtn<-length(which(ReturnNumber==1)) #Total number of first returns

  #Then get average scanangle
  #Scan angles ranged from around -30 to 30
  #took the absolute to prevent negative/positive angles from cancelling out each other
  #Slope could be a problem but here we don't have info to tackle it

  avg_angle<-mean(abs(ScanAngleRank))
  avg_angle<-avg_angle*pi/180 #Degrees to radians

  #now extinction coefficient k
  #The values come from the text in Campbell (1990), note that slightly different numbers were given in the abstract of the paper (1.774 instead of 1.744), but I presume the values in the text is the correct one

  x<-2
  k<-(x^2 + (tan(avg_angle)^2))^0.5 / (x + 1.744*(x+1.182)^(-0.733))

  #Now combine everything together

  PAI<-(-1) * log(1 - weighted_NR/ttl_rtn) * cos(avg_angle) / k

  list(lidar_ePAI = PAI)

}
