######################################################################################################
# Cloud and shadow mask


cloudmask_toa <- function(image) {
  qa <- image$select("BQA")
  mask <- qa$bitwiseAnd(bitwShiftL(1, 4))$eq(0)
  image$updateMask(mask)
}

#  Cloudmask for SR data

cloudmask_sr <- function(image) {
  cloudShadowBitMask <- bitwShiftL(1, 3)
  cloudsBitMask <- bitwShiftL(1, 5)
  qa <- image$select('pixel_qa')
  mask <- qa$bitwiseAnd(cloudShadowBitMask)$eq(0)$
    And(qa$bitwiseAnd(cloudsBitMask)$eq(0))
  image$updateMask(mask)
}
######################################################################################################

######################################################################################################

albedo_func = function(img) {
 
    alb = img$expression(
      expression ="((0.356*blue)+(0.130*red)+(0.373*nir)+(0.085*swir)+(0.072*swir2)- 0.018)/ 1.016",
    opt_map = list(
      'blue'=img$select('B2')$multiply(0.0001),
      'red'=img$select('B4')$multiply(0.0001),
      'nir'= img$select('B5')$multiply(0.0001),
      'swir'= img$select('B6')$multiply(0.0001),
      'swir2'= img$select('B7')$multiply(0.0001)
    ))$rename('albedo')$toFloat()$addBands(img)
    
return(alb$copyProperties(img , list('system:time_start')))

}

######################################################################################################


###   Normalized Difference Vegetation Index (NDVI)

addNDVI <- function(image) {
  ndvi <- image$normalizedDifference(c('B5', 'B4'))$rename('NDVI');
  return(image$addBands(ndvi))
}

##############################################################################################

fvc_func_sobrino=function(img,max=0.5,min=0.2) {
  
  fvc = img$expression(
    expression ="((NDVI-NDVI_min)/(NDVI_max-NDVI_min))**2",
    opt_map = list(
      'NDVI'=img$select('NDVI'),
      NDVI_max = ee$Number(max),## NDVI. max();
      NDVI_min = ee$Number(min) ## NDVI. min();
    ))$rename('FVC')$toFloat()$addBands(img)
  
  return(fvc$copyProperties(img , list('system:time_start')))
  
}
######################################################################################################


EM_func=function(img){ 
EM = img$expression(
  'a*FVC+b*(1-FVC)+((1-b)*(1-FVC)*0.55*a)',
  opt_map = list('FVC' =img$select("FVC"),
                 'a'=ee$Number(0.982),
                 'b'=ee$Number(0.928)
  ))$rename('emissivity')$toFloat()$addBands(img); 
return(EM$copyProperties(img , list('system:time_start')))
}

 ######################################################################################################
 

TOA_func_sr=function(img){ 
  toa = img$expression(
    expression ="((B10*0.1)*coefrad_mult+coefrad_add)-correct_factor",
    opt_map = list(
      'B10' = img$select('B10'),
      'EM' = img$select('emissivity'),
      'coefrad_mult' = ee$Number(0.0003342),     ## Rescaling factor, Band 10
      'coefrad_add' = ee$Number(0.1),            ## Rescaling factor, Band 10
      'correct_factor' = ee$Number(0.29)        ## Correction factor for band 10
      
    ))$rename('TOA')$toFloat()$addBands(img);
  
  return(toa$copyProperties(img , list('system:time_start')))
  
}

  ######################################################################################################

TB_func_sr=function(img){ 
TB = img$expression(
  'k2/log((k1/TOA)+1)',
  opt_map = list('TOA'= img$select('TOA'),
                 'k1'= ee$Number(774.8853), ## Thermal constant, Band 10
                 'k2'= ee$Number(132.0789) ## Thermal constant, Band 10
))$rename('TB')$toFloat()$addBands(img); 
return(TB$copyProperties(img , list('system:time_start')))

}
    
######################################################################################################



LST_func_sr=function(img){   
  LST = img$expression(
    '(TB/(1+(10.8*(TB/14388))*log(EM)))', opt_map = 
      list(
        'TB' = img$select('TB'),
        'EM'= img$select('emissivity')))$rename('LST')$toFloat()$addBands(img)
  
  return(LST$copyProperties(img , list('system:time_start')))
  
}


hotspot_maps=function (x, dist = NULL, p = 0.05) 
{
  if (isFALSE(grepl("0.0001|0.001|0.01|0.05|0.1", p))) {
    stop("P is not one of 0.1, 0.05, 0.01, 0.001 or 0.0001")
  }
  g <- elsa::lisa(x, d1 = 0, d2 = dist, statistic = "localG*")
  data=getValues(g)
  pv <-pnorm(data, mean(data,na.rm =T),sd(data,na.rm =T))
  pv <- as.data.frame(pv)
  FDR <- round(p.adjust(pv[, 1], "fdr"), 7)
  hotspot=ifelse(FDR > 0 & FDR < p, 1, ifelse(FDR < 0 & FDR < p, -1, NA))
  hotspot_ras=g*0+hotspot
  return(hotspot_ras)
}

getis_points=function (x, dist = NULL, p = 0.05) 
{
  if (isFALSE(grepl("0.0001|0.001|0.01|0.05|0.1", p))) {
    stop("P is not one of 0.1, 0.05, 0.01, 0.001 or 0.0001")
  }
  g <- elsa::lisa(x, d1 = 0, d2 = dist, statistic = "G*")
  data=getValues(g)
  pv <- as.data.frame(pnorm(data))
  FDR <- round(p.adjust(pv[, 1], "fdr"), 7)
  g <- round(g, 2)
  shp <- rasterToPoints(g, na.rm = T,spatial=T)
  names(shp) <- "Z.scores"
  shp$FDR<- FDR
  shp$hotspot <- ifelse(shp$Z.scores > 0 & shp$FDR < p, 1, 
                        ifelse(shp$Z.scores < 0 & shp$FDR < p, -1, NA))
  return(shp)
}


# breaks <- c(min(sac.tracts$localg), -2.58, -1.96, -1.65, 1.65, 1.96, 2.58, max(sac.tracts$localg))
#
######################################################################################################
# palette

visParamsLST = list(
  min = 280,
  max = 292,
  bands = "LST",
  palette = c('F2F2F2','EFC2B3','ECB176','E9BD3A','E6E600','63C600','00A600'
  )
  
);



visParamsNDVI = list(
  min = 0.0,
  max = 1,
  bands = "NDVI",
  palette = c(
    'FFFFFF', 'CE7E45', 'DF923D', 'F1B555', 'FCD163', '99B718', '74A901',
    '66A000', '529400', '3E8601', '207401', '056201', '004C00', '023B01',
    '012E01', '011D01', '011301'
  )
)

