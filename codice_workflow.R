######################################################################################################

library(rgee)
library(rgeeExtra)
library(raster)
library(sf)
library(scales)  # Scale functions for visualization
library(RColorBrewer) 
library(elsa) # https://github.com/babaknaimi/elsa

###############################################################################
# setwd("C:/Users/giuli/OneDrive/Desktop/ISPRA_HOTSPOT")

setwd("/home/alf/Scrivania/lav_gee")
source("aux_gee.R")

ee_Initialize(drive = TRUE)
ee_user_info()

##############################################################################

start_date = '2017-06-01';
end_date   = '2021-08-31';


################################################################################
# create useful palette 

cmap1 = c('blue', 'cyan', 'green', 'yellow', 'red');
cmap2 = c('#F2F2F2','#EFC2B3','#ECB176','#E9BD3A','#E6E600','#63C600','#00A600'); 

###############################################################################
# create a region of interest

geometry_piana <- ee$Geometry$Rectangle(
  coords = c(10.8, 43.70 , 11.4, 43.99),
  proj = "EPSG:4326",
  geodesic = FALSE
)

####################################################################################
# load features from gee

Firenze_comune = ee$FeatureCollection("users/giuliaguerri/Firenze_boundaries");
Roma_comune = ee$FeatureCollection("users/giuliaguerri/Roma_boundaries");


# load features from local

Firenze_comune  <- st_read("comuni/Firenze_comune_3035.shp") %>% st_transform("EPSG:4326") %>% sf_as_ee()

Map$centerObject(geometry_piana)



####################################################################################
# Retrieve dem data

db_elev= ee$Image('USGS/SRTMGL1_003');
elevation = db_elev$select('elevation')$clip(geometry_piana)

# ee_print(elevation)


##################################################################################################################################
# Masking sea land

landMask <- ee$Image("CGIAR/SRTM90_V4")$mask()

########################################################################
## Landsat-8 Collection 1 Tier 1

landsat_SR_date <- ee$ImageCollection('LANDSAT/LC08/C01/T1_SR')$
  filterDate(start_date, end_date)$
  filterBounds(geometry_piana)$
  filter(ee$Filter$calendarRange(6,8,'month'))$
  filter(ee$Filter$lte("CLOUD_COVER", 50))$
  filterMetadata("CLOUD_COVER", 'less_than',20)$
  map(cloudmask_sr)
  
 spectral_indices <- landsat_SR_date %>%
 ee$ImageCollection$median() %>% 
 ee_Image_spectralIndex(c("NDVI", "SAVI"))

 names(spectral_indices)

##################################################################################
# Asset consumo suolo

SC = ee$Image('users/andreastrollo/sc2020/SC_LAEA_2020_v2')

# ee_print(SC)

viz <- list(
  max = 230,
  min = 1,
  palette = c("#000000","#FB0E0E")
)

Map$addLayer(SC$clip(geometry_piana), viz, 'SC' );

##############################################################################
# SC reclass

SC_reclass = ee$Image(0)$where(SC$lt(200)$And(SC$neq(2)), 1)$updateMask(SC$gt(0))

# ee_print(SC_reclass )

viz <- list(
  max = 1,
  min = 0,
  palette = c("#000000","#FB0E0E")
)

Map$addLayer(SC_reclass, viz, 'SC_reclass' )

##############################################################################
# Focal

focal = SC_reclass$reduceNeighborhood(
  reducer = ee$Reducer$mean(),
  kernel = ee$Kernel$circle(300, 'meters'))

# ee_print(focal)

focal_reclass = ee$Image(0)$where(focal$gt(0.5), 3)$where(focal$lte(0.5), 2)$where(focal$lte(0.1), 1)$updateMask(SC$gt(0));

# ee_print(focal_reclass )


############################################################################################
# Stima  LST con protocollo ermida
############################################################################################

lsmod <- 'users/sofiaermida/landsat_smw_lst:modules/Landsat_LST.js'
Sys.sleep(2)
message("Ci proviamo...!!")
Sys.sleep(2)
mod <- module(lsmod)
Sys.sleep(2)
LST_ermidacol <- mod$collection("L8", start_date, end_date, Firenze_comune$geometry()$bounds(), TRUE)
LST_ermida=LST_ermidacol$select("LST")$median()

## ee_print(LST_ermida)

##########################################################################################
# masking 

LST_hotspot_EW = LST_ermida$updateMask(focal_reclass$gt(0)) # masking by focal
LST_hotspot_elev_EW= LST_hotspot_EW$updateMask(elevation$lt(300)); # masking by DEM


###########################################################################################


Map$addLayer(LST_ermida$clip(Firenze_comune$geometry()$bounds()),visParamsLST,'LST Ermida Firenze') 

###############################################################################################################################

LST_ermida_ras<- ee_as_raster(
    image = LST_hotspot_elev_EW,
    region = Firenze_comune$geometry()$bounds(),
    scale =100,
    via = "drive"
  )

elevation_FI<- ee_as_raster(
  image = elevation,
  region = Firenze_comune$geometry()$bounds(),
  scale =100,
  via = "drive"
)


writeRaster(LST_ermida_ras,"LST_ermida_ras.tif", overwrite=TRUE)
writeRaster(elevation_FI,"elevation_FI.tif", overwrite=TRUE)


#######################################################################################################

LST_ermida_ras=raster("LST_ermida_ras.tif")
LST_ermida_ras_3035=projectRaster(LST_ermida_ras,crs="+init=epsg:3035")
LST_ermida_ras_3035=raster("LST_ermida_ras_3035.tif")
Firenze_comune  <- st_read("comuni/Firenze_comune_3035.shp")

gstar <- elsa::lisa(LST_ermida_ras_3035, d1 = 0, d2 = 22, statistic = "localG*")

hotspot<- reclassify(gstar, c(-Inf,-2.58,-3,
                               -2.58,-1.96,-2,
                               -1.96,-1.64,-1,
                               -1.64,1.64,0,
                                1.64,1.96,1, # 95 97.5
                                1.96,2.58,2, # 97.5 99
                                2.58,Inf,3)) # 99

plot(hotspot,col=rev(rainbow(7)))

plot(raster("hotspot_arcgis_FI.tif"),col=rev(rainbow(7)))

# labels=c("Cold spot: 99% confidence", 
#           "Cold spot: 95% confidence", 
#           "Cold spot: 90% confidence", 
#           "Not significant",
#           "Hot spot: 90% confidence", 
#           "Hot spot: 95% confidence", 
#           "Hot spot: 99% confidence")


#######################################################################################
# Reference

# [1] Ermida, S.L., Soares, P., Mantas, V., Göttsche, F.-M., Trigo, I.F., 2020. Google Earth Engine open-source code for Land Surface Temperature estimation from the Landsat series. 
#      Remote Sensing, 12 (9), 1471; https://doi.org/10.3390/rs12091471
#      https://github.com/sofiaermida/Landsat_SMW_LST
      
      
# [2] Naimi B, Hamm NA, Groen TA, Skidmore AK, Toxopeus AG, Alibakhshi S (2019). “ELSA: An
#      Entropy-based Local indicator of Spatial Association.” _Spatial Statistics_, *29*, 66-88. doi: 10.1016/j.spasta.2018.10.001 (URL:
#       https://doi.org/10.1016/j.spasta.2018.10.001).
#       https://cran.r-project.org/web/packages/elsa/vignettes/elsa.html


# [3] C Aybar, Q Wu, L Bautista, R Yali and A Barja (2020) rgee: An R package for
#   interacting with Google Earth Engine Journal of Open Source Software URL
#   https://github.com/r-earthengine/rgeeExtra/. 
