lidRmetrics
================
Aland Chan
2025-04-15

## lidRmetrics

This function is used to calculate various LiDAR metrics used in
ecology. The function takes LAS or LAScatalog objects created with the
lidR package as inputs and returns a multi-layered terra SpatRaster,
with each band corresponding to a metric. The choice of metrics, band
names, and definitions are mostly based on Shokirov et al. (2023), with
some minor changes being made for metrics that were not clearly defined
in the paper. Notably, the roughness_L1, roughness_L2, roughness_L3 have
not been included in the current edition. Additionally, canopy_roughness
currently gives the wrong result and needs fixing in the near future.

Shokirov, S. et al. (2023) ‘Habitat highs and lows: Using terrestrial
and UAV LiDAR for modelling avian species richness and abundance in a
restored woodland’, Remote Sensing of Environment. Elsevier, 285,
p. 113326. doi: 10.1016/J.RSE.2022.113326.

## Installation

``` r
library(devtools)

devtools::install_github("alandhyc/lidRmetrics",force = T)
```

## Example

Applying the function across a las catalog, writing `.tif` files in
`output_dir`

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

lidRmetrics::LiDAR_metrics(las = ctg,
                           res=5,
                           h_cutoff=1.3,
                           mcc_s=1.5,
                           mcc_t=0.3)
```

Processing a .las file

``` r
mylas<-readLAS("/path/to/lasfile.las")

metrics<-lidRmetrics::LiDAR_metrics(las = mylas,
                                    res=5,
                                    h_cutoff=1.3,
                                    mcc_s=1.5,
                                    mcc_t=0.3)
```

## Contact

For questions or suggestions, open an issue or contact `hyc43@cam.ac.uk`
