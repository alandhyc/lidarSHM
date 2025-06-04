lidarSHM
================
Aland Chan
2025-06-04

## lidarSHM

This package is designed for calculating various lidar metrics from
lidar point clouds and canopy height models.

The package has three main functions:

1.  `lidar_SHM` (in development) calculates a shrub height model (SHM)
    from lidar point cloud data.

2.  `CHM_metrics` calculates a number of useful structural metrics from
    a canopy height model or shrub height model in one single step. The
    function takes a CHM SpatRaster and a vector of metric names;
    returns a SpatRaster with multiple layers (names of layers
    corresponding to metric names). Definitions of metrics are based on
    Zhang et al. (2024) with several minor modifications and additions.

3.  `liDAR_metrics` calculates a number of useful structural metrics
    from point clouds in one single step. The function takes either a
    LAS or LAScatalog object created with the `lidR::` package as
    inputs. If a `LAS` object is supplied, the function returns a
    multi-layered SpatRaster, with layers corresponding to metrics users
    requested. Alternatively, if a `LAScatalog` object is supplied, it
    writes the SpatRasters as `.tif` files in folder specified by
    `opt_output_files()`. By default, the function performs mcc ground
    classification, but users could also supply a ground-classified
    point cloud if they want to skip this time-consuming step or prefer
    another ground classification algorithm. For the function to operate
    properly, <las@data> should have the columns named `Z` (height),
    `ReturnNumber`, `NumberOfReturns`, and `ScanAngleRank` (angle in
    degrees). The last two are specifically used for ePAI calculations.
    The definitions and names of metrics are mostly based on Shokirov et
    al. (2023).

## List of metrics

Some of the groupings and definitions are derived from Zhang et
al. (2024) and Shokirvo et al. (2023).

| Type | Derived.from | Metric | Code | Description | Function |
|:---|:---|:---|:---|:---|:---|
| Shrub metrics | LiDAR | Shrub height map (SHM) | SHM | Canopy height model of shrub layer (height \< 1.3m). We use our own ground classification method to ensure that the shrubs are preserved | lidar_SHM |
| Height metrics | LiDAR | LiDAR mean height | lidar_meanH | Mean height of canopy (points \> 1.3 m). | LiDAR_metrics |
| Height metrics | CHM | CHM max height | chm_maxH | Maximum height of canopy (99.9th percentile) | CHM_metrics |
| Height metrics | LiDAR | LiDAR max height | lidar_maxH | Maximum height of canopy (99.9th percentile; points \> 1.3 m). | LiDAR_metrics |
| Height metrics | LiDAR | LiDAR canopy height percentiles | lidar_p_05, lidar_p_10, lidar_p_25, lidar_p_50, lidar_p_75, lidar_p_90, lidar_p_95, lidar_p_999 | Canopy height percentiles (points \> 1.3 m). Canopy height percentiles are the height below which a specified percentage of total point clouds were located. For example, p_05 = 2 m means that 5% of vegetation points are found below 2 m. lidar_p_999 represents the 99.9th quantile, which is equivalent to max height but with extreme values filtered out. | LiDAR_metrics |
| Height distribution metrics | CHM | Standard deviation of CHM | chm_sdH | Standard deviation of canopy height from CHM | CHM_metrics |
| Height distribution metrics | LiDAR | Standard deviation of LiDAR | lidar_stdH | Standard deviation of canopy height (points \> 1.3 m) | LiDAR_metrics |
| Height distribution metrics | CHM | Coefficient of variation of CHM | chm_height_cv | Coefficient of variation of canopy height from CHM | CHM_metrics |
| Height distribution metrics | LiDAR | Coefficient of variation of LiDAR | lidar_height_cv | Coefficient of variation of height, (points \> 1.3 m). Indicates the canopy height variation. | LiDAR_metrics |
| Heigh distribution metrics | CHM | Robust coefficient of variation of CHM | chm_rcv | The coefficient of variation is often not very stable when you have extreme values, the robust variant of the coefficient of variation use the interquartile range and median instead. | CHM_metrics |
| Heigh distribution metrics | CHM | Robust coefficient of variation of lidar points | lidar_rcv | The coefficient of variation is often not very stable when you have extreme values, the robust variant of the coefficient of variation use the interquartile range divided by the median instead. | CHM_metrics |
| Height distribution metrics | CHM | Root mean square of CHM | chm_rms | Similar to standard deviation as a measure of height distribution (spread out or close to mean). Formula: sqrt(mean((zi - zm)^2)) | CHM_metrics |
| Height distribution metrics | LiDAR | Root mean square of LiDAR | lidar_rms | Similar to standard deviation as a measure of height distribution (spread out or close to mean). Formula: sqrt(mean((zi - zm)^2)) | LiDAR_metrics |
| Height distribution metrics | CHM | Skew of canopy heights from CHM | chm_skewH | Skewness of canopy height (points \> 1.3 m). Negative skewness means that the distribution is dominated by higher points (upper canopy is dominant) but a few extreme lower points. Positive skewness means that the distribution dominated by lower points (lower canopy is dominant) but a few extreme higher points. | CHM_metrics |
| Height distribution metrics | LiDAR | Skew of canopy heights from LiDAR | lidar_skewH | Skewness of canopy height (points \> 1.3 m). Negative skewness means that the distribution is dominated by higher points (upper canopy is dominant) but a few extreme lower points. Positive skewness means that the distribution dominated by lower points (lower canopy is dominant) but a few extreme higher points. | LiDAR_metrics |
| Height distribution metrics | LiDAR | Kurtosis of LiDAR point cloud | lidar_kurH | Kurtosis of canopy height (points \> 1.3 m). Negative kurtosis means the distribution of points centered around the mean (mid-canopy is dominant). Positive kurtosis means the point distribution is heavy on tails and less around the mean (lower and upper canopy is dominant). | LiDAR_metrics |
| Canopy openness metrics | LiDAR | Canopy cover | lidar_Cov | Fraction of canopy cover, (Fraction of first returns \> 1.3 m) | LiDAR_metrics |
| Canopy openness metrics | CHM | Gap fraction from CHM | chm_gapFrac | Proportion of pixels \<2m height | CHM_metrics |
| Canopy openness metrics | LiDAR | Gap fraction from LiDAR points | lidar_gapFrac | Proportion of lidar first returns \<2m in height | LiDAR_metrics |
| Canopy openness metrics | CHM | Gap fraction from CHM | chm_grndFrac | Proportion of pixels \<0.5m in height | CHM_metrics |
| Canopy openness metrics | LiDAR | Gap fraction from LiDAR points | lidar_grndFrac | Proportion of lidar first returns classified as ground | LiDAR_metrics |
| Vegetation complexity metrics | CHM | Normalised rumple index | chm_rumple_norm | Ratio of canopy surface area to projected ground area derived from a height-normalized CHM | CHM_metrics |
| Vegetation complexity metrics | LiDAR | Effective plant area index (m2/m2) | lidar_ePAI | Estimating surface area of plant material. In LiDAR contexts this is related to the extinction coefficient in Beer Lambert law. See Kwak et al. (2010) and Zhu et al. (2020). | LiDAR_metrics |
| Vegetation complexity metrics | LiDAR | Vertical complexity index | lidar_vci_2m, lidar_vci_5m, lidar_vci_10m, lidar_vci_15m, lidar_vci_20m | Vertical complexity indexes (VCI) at 2 m, 5 m, 10 m, 15 m, 20 m height bins, (points \> 1.3 m). | LiDAR_metrics |
|  |  |  |  | VCI = (∑i=1HB\[(pi ln (pi))\])/ ln (HB) |  |
|  |  |  |  | Where VCI in a vertical complexity index, HB is the total number of height bins, and pi is the proportional abundance of LiDAR returns in height bin i. |  |
|  |  |  |  | A VCI value close to one indicates that most height bins have an equal amount of vegetation. VCI value decreases if the distribution of canopy in the height bin becomes more uneven (van Ewijk et al., 2011). |  |
| Vegetation complexity metrics | LiDAR | Normalised Shannon diversity index | lidar_canopy_shannon | Normalized Shannon diversity index of canopy (Pretzsch, 2009), (points \> 1.3 m). Indicates canopy height diversity. | LiDAR_metrics |
| Vegetation complexity metrics | LiDAR | Total vegetation volume (m3) | lidar_Tvolume | Total vegetation volume (m3) – number of 0.5 m3 voxels divided by 8 (ground points excluded). | LiDAR_metrics |
| Vegetation complexity metrics | LiDAR | Vegetation volume at 0-1m height | lidar_vlayer_L1 | Vegetation volume (m3) in 1st layer (points 0-1 m¸ ground points excluded). | LiDAR_metrics |
| Vegetation complexity metrics | LiDAR | Vegetation volume at 1-10m height | lidar_vlayer_L2 | Vegetation volume (m3) in 2st layer (points 1 m–10 m). | LiDAR_metrics |
| Vegetation complexity metrics | LiDAR | Vegetation volume at \>10m height | lidar_vlayer_L3 | Vegetation volume (m3) in 3st layer (points 10 m and above). | LiDAR_metrics |
| Vegetation complexity metrics | LiDAR | Mean height at 0-1m, 1-10m, and \>10m strata | lidar_meanH_L1, lidar_meanH_L2, lidar_meanH_L3 | Mean height of 1st, 2nd, 3rd layer. | LiDAR_metrics |
| Vegetation complexity metrics | LiDAR | SD of height at 0-1m, 1-10m, and \>10m strata | lidar_sdH_L1, lidar_sdH_L2, lidar_sdH_L3 | Standard deviation of vegetation height in 1st, 2nd, 3rd layer. | LiDAR_metrics |
| Vegetation complexity metrics | LiDAR | VCI at 0-1m, 1-10m, and \>10m strata | lidar_vci_L1, lidar_vci_L2, lidar_vci_L3 | Vertical complexity indexes of 1st, 2nd, 3rd layer (van Ewijk et al., 2011). Vertical distribution of vegetation across different layers. | LiDAR_metrics |

## Installation

``` r
library(devtools)

devtools::install_github("alandhyc/lidarSHM",force = T)
```

## Examples

### CHM_metrics

Calculating CHM metrics from a CHM.

``` r
library(terra)

#

chm_r<-terra::rast("path/to/chm")

#Calculating all metrics
#Aggregate 50x50 pixels into one (resotuion 1m for CHM, 50m for metrics)
#Defining "gap" as things under 2m

chm_metrics<-CHM_metrics(chm = chm_r,
                         agg = 50,
                         gap_thres = 2,
                         metrics = "all")

#Calculating specific metrics
#Just do max height and standard deviation

chm_metrics_subset<-CHM_metrics(chm = chm_r,
                                agg = 50,
                                gap_thres = 2,
                                metrics = c("chm_maxH","chm_sdH"))
```

### LiDAR_metrics

Applying the function across a `LAS` object

``` r
library(lidR)

las<-readLAS("/path/to/las/file.las")

#With all options
#Perform ground classification
#Resolution at 50m
#Calculate most metrics based on points >1.3m ("canopy" definition)
#Anything <2m treated as a gap
#ground classification with mcc using mcc_s as s and mcc_t as t
#Maximum height of tree approximately 35m
#height classes -1-2, 2-5, 5-10, 10-15, 15-35 when calculating shannon index
#voxel resolution 0.5 m3 when calculating canopy volumn 
#Stratifying point cloud into three layers (0-1, 1-10, 10-35) when calculated layer-specific metrics in Shakirov et al. (2023)
#Calculate all metrics

output_metrics<-LiDAR_metrics(las,
                              ground_classified = F,
                              res=50,
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
                              )

#Calculating specific metrics
#Point cloud already ground classified

output_metrics<-LiDAR_metrics(las,
                              ground_classified = T,
                              res=50,
                              zmax = 35,
                              metrics = c("lidar_maxH","lidar_ePAI")
                              )
```

Alternatively, we could apply the function to a `LAScatalog` object with
parallel computing and output written to a folder.

``` r
library(lidR)
library(lidRmetrics)
library(future)

#Set up parallel computation
ncores<-5
plan(multisession, workers = ncores)
set_lidr_threads(ncores)

#Load las catalog

ctg<-readLAScatalog("/path/to/las/files")

output_dir<-"/path/to/output/directory"


opt_output_files(ctg)<-paste0(output_dir,"/{ORIGINALFILENAME}")
opt_chunk_buffer(ctg)<-20 #Set a buffer

LiDAR_metrics(las,
              ground_classified = F,
              res=50,
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
              )
```

## Acknowledgements

Evie Huhtala for providing part of the code in the `CHM_metrics()`
function

## References

Zhang, B., Fischer, F.J., Prober, S.M., Yeoh, P.B., Gosper, C.R.,
Zdunic, K. and Jucker, T. (2024), Robust retrieval of forest canopy
structural attributes using multi-platform airborne LiDAR. Remote Sens
Ecol Conserv, 10: 725-742. <https://doi.org/10.1002/rse2.398>

Shokirov, S. et al. (2023) ‘Habitat highs and lows: Using terrestrial
and UAV LiDAR for modelling avian species richness and abundance in a
restored woodland’, Remote Sensing of Environment. Elsevier, 285,
p. 113326. doi: 10.1016/J.RSE.2022.113326.

## Contact

For questions or suggestions, open an issue or contact `hyc43@cam.ac.uk`
