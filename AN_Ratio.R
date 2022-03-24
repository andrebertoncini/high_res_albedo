#This code is used to retrive the MODIS Albedo-to-Nadir (AN) Reflectance ratios,
#which is based on Shuai et al. (2011) "An algorithm for the retrieval of 30-m 
#snow-free albedo from Landsat surface reflectance and MODIS BRDF" - Remote Sensing
#of Environment.

#IMPORTANT: If you choose to use the RTLSR_Kernel.R script, you will need to
#comment out all the lines that involve the snow kernel and the f_iso parameter.

#Author: Andre Bertoncini

#Institution: Centre for Hydrology - University of Saskatchewan

library(XML)
library(pracma)
library(raster)


#Set target date (same as Sentinel-2 image)

Dm <- "08"
M <- "08"
Y <- "2018"

img_date <- paste0(Y, M, Dm)
modis_date <- paste0(Y, ".", M, ".", Dm)
Dm_j <- (as.numeric(Dm) - 32) + (floor(275*(as.numeric(M)/9))) + (2*floor(3/(as.numeric(M)+1)))
Dm_j <- ifelse(Y == "2016" | Y == "2020", Dm_j + 1, Dm_j) #for leap years: 2016 and 2020 
date_j <- paste0(Y, Dm_j)


#Set folders for Sentinel-2 metadata

list_folders <- list.dirs("/path to Sentinel-2 folder")

folder <- grep(pattern = paste0(img_date), list_folders)


#Generate 16-day MODIS surface reflectance timeseries

setwd("/path to MODIS images")

AOI_sent <- extent(399960, 509760, 5690220, 5800020)

modis_list <- list()

for (i in 0:15) {
  
  date_list <- as.character((as.numeric(date_j) - 8) + i)  
  
  modis_list <- append(modis_list, list.files(getwd(), pattern = paste0("MOD09GA.*.", date_list, ".*.tif$"), recursive = T))
  
}

modis_list <- unlist(modis_list)


for (w in c(1,2,3,4,6,7)) {
  
  band <- w
  
  
  setwd("/path to MODIS images")  
  
  ref_i <- stack()
  
  for (i in 0:15) {
    
    modis_ref <- raster(modis_list[i*12 + (band + 5)])
    
    modis_ref <- crop(modis_ref, AOI_sent)
    
    modis_ref[modis_ref == -28672] <- NA
    modis_ref[modis_ref > 16000] <- NA
    modis_ref[modis_ref < -100] <- NA
    modis_ref <- modis_ref*0.0001
    
    plot(modis_ref)
    
    ref_i <- stack(ref_i, modis_ref)
  }  
}




pb <- txtProgressBar(min = 1, max = 6, style = 3)

for (w in c(1,2,3,4,6,7)) {
  
  band = w
  
  setwd("/path to output AN ratios")
  
  k1 <- 1.247
  k2 <- 1.186
  k3 <- 5.157
  alpha <- 0.3


  #Generate geometry information from Sentinel-2 XML data

  doc = xmlParse(paste0(list_folders[[folder[7]]], "/MTD_TL.xml"))
  #doc = xmlParse(paste0(list_folders[[127]], "/MTD_TL.xml"))

  sun_zenith <- xmlToDataFrame(getNodeSet(doc, "//Tile_Angles/Sun_Angles_Grid/Zenith/Values_List"))

  sun_zenith <- strsplit(t(sun_zenith), "[ ]")

  sun_zenith <- as.numeric(unlist(sun_zenith))

  sun_zenith_m <- t(matrix(sun_zenith, ncol = 23, nrow = 23))

  theta_s_i <- raster(sun_zenith_m)

  extent(theta_s_i) <- modis_ref

  theta_s_i <- resample(theta_s_i, modis_ref)*(pi/180)

  plot(theta_s_i)


  sun_azimuth <- xmlToDataFrame(getNodeSet(doc, "//Tile_Angles/Sun_Angles_Grid/Azimuth/Values_List"))

  sun_azimuth <- strsplit(t(sun_azimuth), "[ ]")

  sun_azimuth <- as.numeric(unlist(sun_azimuth))

  sun_azimuth_m <- t(matrix(sun_azimuth, ncol = 23, nrow = 23))

  phi_s <- raster(sun_azimuth_m)

  extent(phi_s) <- modis_ref

  phi_s <- resample(phi_s, modis_ref)*(pi/180)

  plot(phi_s)


  view_zenith <- xmlToDataFrame(getNodeSet(doc, "//Tile_Angles/Viewing_Incidence_Angles_Grids/Zenith/Values_List"))

  view_zenith <- strsplit(t(view_zenith), "[ ]")

  view_zenith <- as.numeric(unlist(view_zenith))

  mean_v_vector <- vector()

  for (i in 1:length(view_zenith)) {

    mean_v_vector[i] <- mean(view_zenith[c(i, i+529, i+(2*529), i+(3*529), i+(4*529), i+(5*529), i+(6*529),
                                           i+(7*529), i+(8*529), i+(9*529), i+(10*529), i+(11*529), i+(12*529),
                                           i+(13*529), i+(14*529), i+(15*529), i+(16*529), i+(17*529), i+(18+529),
                                           i+(19*529), i+(20*529), i+(21*529), i+(22*529))], na.rm = T)
  }

  view_zenith_m <- t(matrix(mean_v_vector, ncol = 23, nrow = 23))

  theta_v <- raster(view_zenith_m)

  extent(theta_v) <- modis_ref

  theta_v <- resample(theta_v, modis_ref)*(pi/180)

  plot(theta_v)


  view_azimuth <- xmlToDataFrame(getNodeSet(doc, "//Tile_Angles/Viewing_Incidence_Angles_Grids/Azimuth/Values_List"))

  view_azimuth <- strsplit(t(view_azimuth), "[ ]")

  view_azimuth <- as.numeric(unlist(view_azimuth))

  mean_p_vector <- vector()

  for (i in 1:length(view_azimuth)) {

    mean_p_vector[i] <- mean(view_azimuth[c(i, i+529, i+(2*529), i+(3*529), i+(4*529), i+(5*529), i+(6*529),
                                            i+(7*529), i+(8*529), i+(9*529), i+(10*529), i+(11*529), i+(12*529),
                                            i+(13*529), i+(14*529), i+(15*529), i+(16*529), i+(17*529), i+(18+529),
                                            i+(19*529), i+(20*529), i+(21*529), i+(22*529))], na.rm = T)
  }

  view_azimuth_m <- t(matrix(mean_p_vector, ncol = 23, nrow = 23))

  phi_v <- raster(view_azimuth_m)

  extent(phi_v) <- modis_ref

  phi_v <- resample(phi_v, modis_ref)*(pi/180)


  phi = abs(phi_s - phi_v)

  phi_f <- phi

  phi_f[phi < pi] <- NA

  phi[phi >= pi] <- NA

  phi_f <- abs(phi_f - (2*pi))

  phi_final <- cover(phi, phi_f)

  phi <- phi_final

  plot(phi)
  
  
  #Import f_iso, f_vol, f_geo, and f_snow parameters and evaluation metrics
  
  f_iso <- raster(paste0("/path to ouput MODIS BRDF kernel parameters and evaluation/f_iso_", img_date, "_B", band, ".tif"))
  
  f_vol <- raster(paste0("/path to ouput MODIS BRDF kernel parameters and evaluation/f_vol_", img_date, "_B", band, ".tif"))
  
  f_geo <- raster(paste0("/path to ouput MODIS BRDF kernel parameters and evaluation/f_geo_", img_date, "_B", band, ".tif"))
  
  f_snow <- raster(paste0("/path to ouput MODIS BRDF kernel parameters and evaluation/f_snow_", img_date, "_B", band, ".tif"))
  
  rmse <- raster(paste0("/path to ouput MODIS BRDF kernel parameters and evaluation/rmse_", img_date, "_B", band, ".tif"))
  
  wod_wdr <- raster(paste0("/path to ouput MODIS BRDF kernel parameters and evaluation/wod_wdr_", img_date, "_B", band, ".tif"))
  
  wod_wsa <- raster(paste0("/path to ouput MODIS BRDF kernel parameters and evaluation/wod_wsa_", img_date, "_B", band, ".tif"))
  
  
  #Remove pixels that do not meet MODIS thresholds following Shuai et al. (2008) "Quality assessment of BRDF/albedo 
  #retrievals in MODIS operational system" - Geophysical Research Letters
  
  f_iso[rmse > 0.08] <- NA
  f_iso[wod_wdr > 1.65] <- NA
  f_iso[wod_wsa > 2.5] <- NA
  
  f_vol[rmse > 0.08] <- NA
  f_vol[wod_wdr > 1.65] <- NA
  f_vol[wod_wsa > 2.5] <- NA
  
  f_geo[rmse > 0.08] <- NA
  f_geo[wod_wdr > 1.65] <- NA
  f_geo[wod_wsa > 2.5] <- NA
  
  f_snow[rmse > 0.08] <- NA
  f_snow[wod_wdr > 1.65] <- NA
  f_snow[wod_wsa > 2.5] <- NA
  
  
  #Calculation of BRDF reflectance (R_omega) at Sentinel-2 observation and illumination geometry
  #(theta_s, theta_v, phi).
  
  theta_s <- theta_s_i


  cos_E = ((sin(theta_s))*(sin(theta_v))*(cos(phi))) + ((cos(theta_s))*(cos(theta_v)))

  E = (acos(cos_E))*(180/pi)

  P_E = (11.1*(exp(-0.087*(180 - E)))) + (1.1*(exp(-0.014*(180 - E))))

  R_0 = (k1 + (k2*(cos(theta_s) + cos(theta_v))) + (k3*(cos(theta_s)*cos(theta_v))) + P_E)/(4*(cos(theta_s) + cos(theta_v)))

  theta_st = atan(tan(theta_s))

  theta_vt = atan(tan(theta_v))

  efe = sqrt(((tan(theta_st))^2) + ((tan(theta_vt))^2) - (2*(tan(theta_st)*tan(theta_vt)*cos(phi))))

  cos_k = 2*((sqrt((efe^2) + ((tan(theta_st)*tan(theta_vt)*sin(phi))^2)))/((1/cos(theta_st))+(1/cos(theta_vt))))

  cos_k[cos_k < -1] <- -1
  cos_k[cos_k > 1] <- 1

  k = (acos(cos_k))

  over = (1/pi)*(k - (sin(k)*cos_k))*((1/cos(theta_st)) + (1/cos(theta_vt)))

  cos_Et = (sin(theta_st)*sin(theta_vt)*cos(phi)) + (cos(theta_st)*cos(theta_vt))
  
  
  k_snow_omega = (R_0*(1-(alpha*(cos(E*(pi/180)))*(exp(-cos(E*(pi/180))))))) + ((0.4076*alpha) - 1.1081)
  
  k_vol_omega = (((((pi/2) - (E*(pi/180)))*cos_E) + sin(E*(pi/180)))/(cos(theta_s) + cos(theta_v))) - (pi/4)

  k_geo_omega = over - ((1/cos(theta_st)) + (1/cos(theta_vt)) - ((1/2)*(1 + cos_Et)*((1/cos(theta_vt))*(1/cos(theta_st)))))
  
  
  R_omega = f_iso + f_vol*k_vol_omega + f_geo*k_geo_omega + f_snow*k_snow_omega
  
  R_omega[R_omega < 0] <- NA

  plot(R_omega)


  #Calculation of BRDF reflectance (R_l_theta_s) at Sentinel-2 illumination angle (theta_s).

  rm(theta_v)
  rm(phi)


  k_vol_theta_s <- function(theta_v, phi) {

    cos_E = ((sin(theta_s))*(sin(theta_v))*(cos(phi))) + ((cos(theta_s))*(cos(theta_v)))

    E = (acos(cos_E))*(180/pi)


    k_vol_omega = ((((((pi/2) - (E*(pi/180)))*cos_E) + sin(E*(pi/180)))/(cos(theta_s) + cos(theta_v))) - (pi/4))*(sin(theta_v))*(cos(theta_v))

  }



  k_geo_theta_s <- function(theta_v, phi) {

    theta_st = atan(tan(theta_s))

    theta_vt = atan(tan(theta_v))

    efe = sqrt(((tan(theta_st))^2) + ((tan(theta_vt))^2) - (2*(tan(theta_st)*tan(theta_vt)*cos(phi))))

    cos_k = 2*((sqrt((efe^2) + ((tan(theta_st)*tan(theta_vt)*sin(phi))^2)))/((1/cos(theta_st))+(1/cos(theta_vt))))

    cos_k[cos_k < -1] <- -1
    cos_k[cos_k > 1] <- 1

    k = (acos(cos_k))

    over = (1/pi)*(k - (sin(k)*cos_k))*((1/cos(theta_st)) + (1/cos(theta_vt)))

    cos_Et = (sin(theta_st)*sin(theta_vt)*cos(phi)) + (cos(theta_st)*cos(theta_vt))


    k_geo_omega = (over - ((1/cos(theta_st)) + (1/cos(theta_vt)) - ((1/2)*(1 + cos_Et)*((1/cos(theta_vt))*(1/cos(theta_st))))))*(sin(theta_v))*(cos(theta_v))

  }
  
  
  
  k_snow_theta_s <- function(theta_v, phi) {
    
    cos_E = ((sin(theta_s))*(sin(theta_v))*(cos(phi))) + ((cos(theta_s))*(cos(theta_v)))
    
    E = (acos(cos_E))*(180/pi)
    
    P_E = (11.1*(exp(-0.087*(180 - E)))) + (1.1*(exp(-0.014*(180 - E))))
    
    R_0 = (k1 + (k2*(cos(theta_s) + cos(theta_v))) + (k3*(cos(theta_s)*cos(theta_v))) + P_E)/(4*(cos(theta_s) + cos(theta_v)))
    
    
    k_snow_omega = ((R_0*(1-(alpha*(cos(E*(pi/180)))*(exp(-cos(E*(pi/180))))))) + ((0.4076*alpha) - 1.1081))*(sin(theta_v))*(cos(theta_v))
    
  }
  
  
  theta_s_vec <- as.vector(theta_s_i)
  
  k_vol_l_theta_s <- vector()
  
  k_geo_l_theta_s <- vector()
  
  k_snow_l_theta_s <- vector()
  
  
  for (i in 1:length(theta_s_vec)) {
    
    theta_s <- theta_s_vec[i]
    
    k_vol_l_theta_s[i] <- (integral2(k_vol_theta_s, 0, pi/2, 0, 2*pi, reltol = 1e-2)$Q)*(1/pi)
    
  }
  
  
  for (i in 1:length(theta_s_vec)) {
    
    theta_s <- theta_s_vec[i]
    
    k_geo_l_theta_s[i] <- (integral2(k_geo_theta_s, 0, pi/2, 0, 2*pi, reltol = 1e-2)$Q)*(1/pi)
    
  }
  
  
  for (i in 1:length(theta_s_vec)) {
    
    theta_s <- theta_s_vec[i]
    
    k_snow_l_theta_s[i] <- (integral2(k_snow_theta_s, 0, pi/2, 0, 2*pi, reltol = 1e-2)$Q)*(1/pi)
    
  }
  
  
  k_vol_l_theta_raster <- raster(matrix(k_vol_l_theta_s, ncol = 237, nrow = 237))
  
  extent(k_vol_l_theta_raster) <- modis_ref
  
  
  k_geo_l_theta_raster <- raster(matrix(k_geo_l_theta_s, ncol = 237, nrow = 237))
  
  extent(k_geo_l_theta_raster) <- modis_ref
  
  
  k_snow_l_theta_raster <- raster(matrix(k_snow_l_theta_s, ncol = 237, nrow = 237))
  
  extent(k_snow_l_theta_raster) <- modis_ref
  
  
  R_l_theta_s = f_iso + f_vol*k_vol_l_theta_s + f_geo*k_geo_l_theta_s + f_snow*k_snow_l_theta_s
  
  R_l_theta_s[R_l_theta_s < 0] <- NA
  
  plot(R_l_theta_s)
  
  
  #Calculation of BRDF reflectance (R_l) integrated for both the illumination and observation hemispheres.

  rm(theta_s)


  k_vol <- function(theta_s, theta_v, phi) {

    cos_E = ((sin(theta_s))*(sin(theta_v))*(cos(phi))) + ((cos(theta_s))*(cos(theta_v)))

    E = (acos(cos_E))*(180/pi)


    k_vol_omega = (((((((pi/2) - (E*(pi/180)))*cos_E) + sin(E*(pi/180)))/(cos(theta_s) + cos(theta_v))) - (pi/4))*(sin(theta_v))*(cos(theta_v)))*(sin(theta_s))*(cos(theta_s))

  }



  k_geo <- function(theta_s, theta_v, phi) {

    theta_st = atan(tan(theta_s))

    theta_vt = atan(tan(theta_v))

    efe = sqrt(((tan(theta_st))^2) + ((tan(theta_vt))^2) - (2*(tan(theta_st)*tan(theta_vt)*cos(phi))))

    cos_k = 2*((sqrt((efe^2) + ((tan(theta_st)*tan(theta_vt)*sin(phi))^2)))/((1/cos(theta_st))+(1/cos(theta_vt))))

    cos_k[cos_k < -1] <- -1
    cos_k[cos_k > 1] <- 1

    k = (acos(cos_k))

    over = (1/pi)*(k - (sin(k)*cos_k))*((1/cos(theta_st)) + (1/cos(theta_vt)))

    cos_Et = (sin(theta_st)*sin(theta_vt)*cos(phi)) + (cos(theta_st)*cos(theta_vt))


    k_geo_omega = ((over - ((1/cos(theta_st)) + (1/cos(theta_vt)) - ((1/2)*(1 + cos_Et)*((1/cos(theta_vt))*(1/cos(theta_st))))))*(sin(theta_v))*(cos(theta_v)))*(sin(theta_s))*(cos(theta_s))

  }
  
  
  k_snow <- function(theta_s, theta_v, phi) {
    
    cos_E = ((sin(theta_s))*(sin(theta_v))*(cos(phi))) + ((cos(theta_s))*(cos(theta_v)))
    
    E = (acos(cos_E))*(180/pi)
    
    P_E = (11.1*(exp(-0.087*(180 - E)))) + (1.1*(exp(-0.014*(180 - E))))
    
    R_0 = (k1 + (k2*(cos(theta_s) + cos(theta_v))) + (k3*(cos(theta_s)*cos(theta_v))) + P_E)/(4*(cos(theta_s) + cos(theta_v)))
    
    
    k_snow_omega = (((R_0*(1-(alpha*(cos(E*(pi/180)))*(exp(-cos(E*(pi/180))))))) + ((0.4076*alpha) - 1.1081))*(sin(theta_v))*(cos(theta_v)))*(sin(theta_s))*(cos(theta_s))
    
  }
  
  
  k_vol_l <- (integral3(k_vol, 0, pi/2, 0, pi/2, 0, 2*pi, reltol = 1e-2))*(2/pi)
  
  k_geo_l <- (integral3(k_geo, 0, pi/2, 0, pi/2, 0, 2*pi, reltol = 1e-2))*(2/pi)
  
  k_snow_l <- (integral3(k_snow, 0, pi/2, 0, pi/2, 0, 2*pi, reltol = 1e-2))*(2/pi)
  
  
  R_l = f_iso + f_vol*k_vol_l + f_geo*k_geo_l + f_snow*k_snow_l
  
  R_l[R_l < 0] <- NA
  
  plot(R_l)
  
  
  
  #Calculation of BSA AN ratio (a_bsa)
  
  a_bsa = R_l_theta_s/R_omega
  
  plot(a_bsa)
  
  writeRaster(a_bsa, filename = paste0("a_bsa_", img_date, "_B", band), format = "GTiff", overwrite = T)
  
  
  #Calculation of WSA AN ratio (a_wsa)
  
  a_wsa = R_l/R_omega
  
  plot(a_wsa)
  
  writeRaster(a_wsa, filename = paste0("a_wsa_", img_date, "_B", band), format = "GTiff", overwrite = T)
  
  Sys.sleep(0.0001)
  setTxtProgressBar(pb, w)
  
}

