#This code is used to retrieve MODIS BRDF using the snow kernel model, 
#which is based on Jiao et al. (2019) "Development of a snow kernel to 
#better model the anisotropic reflectance of pure snow in a kernel-driven 
#BRDF model framework" - Remote Sensing of Environment.

#Author: Andre Bertoncini

#Institution: Centre for Hydrology - University of Saskatchewan


library(raster)
library(rgdal)
library(reshape2)
library(ExtDist)
library(pracma)
library(MASS)
#library(MODIS)


#Set target date (same as Sentinel-2 image)

Dm <- "08"
M <- "08"
Y <- "2018"

img_date <- paste0(Y, M, Dm)
modis_date <- paste0(Y, ".", M, ".", Dm)
Dm_j <- (as.numeric(Dm) - 32) + (floor(275*(as.numeric(M)/9))) + (2*floor(3/(as.numeric(M)+1)))
Dm_j <- ifelse(Y == "2016" | Y == "2020", Dm_j + 1, Dm_j) #for leap years: 2016 and 2020
date_j <- paste0(Y, Dm_j)


#Set Sentinel-2 sensor (can be S2A or S2B)

sensor <- "S2A" or "S2B"


#Uncomment lines below to download MODIS data through the MODIS package.

#begin_m <- "2019.08.28"
#end_m <- "2019.09.12"


#Retrieve MODIS data

AOI_sent <- extent(399960, 509760, 5690220, 5800020)

AOI <- extent(421290, 496851.5, 5759734, 5800008.8)

modis_ext <- extent(-120, -113, 50, 53)


#modis_products <- getProduct()

#modis_tile <- getTile(modis_ext)

#modis_raw <- getHdf("MOD09GA", begin = begin_m, end = end_m, tileH = 10, tileV = 3, 
#                    quiet = TRUE)


#modis_geo <- runMrt(product = "MOD09GA", collection = "006", begin = begin_m, 
#                    end = end_m, tileH = 10, tileV = 3, 
#                    SDSstring = "0 1 1 1 0 1 1 0 0 0 0 1 1 1 1 1 1 1",
#                    outProj = "UTM", datum = "WGS84", zone = 11, job = "MOD", mosaic = FALSE)



#Generate 16-day MODIS surface reflectance timeseries

#If using MODIS data on your computer, set the working directory to your MOD09GA images.
#Remember that you need 16 images to generate BRDF information (the 9th day is the target date).

setwd("/path to MODIS images")

modis_list <- list()

for (i in 0:15) {
  
  date_list <- as.character((as.numeric(date_j) - 8) + i)  
  
  modis_list <- append(modis_list, list.files(getwd(), pattern = paste0("MOD09GA.*.", date_list, ".*.tif$"), recursive = T))
  
}

modis_list <- unlist(modis_list)



for (w in c(1,2,3,4,6,7)) {
  
  band <- w
  
  
  setwd("/path to MODIS images")  
  
  
  #Retrieve surface reflectance
  
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
  
  
  #Generate cloud and smoke masks for MODIS
  
  modis_mask_i <- stack()
  
  for (i in 0:15) {
    
    modis_QA_500 <- raster(modis_list[i*12 + 5])
    
    modis_QA_500 <- as.vector(crop(modis_QA_500, AOI_sent))
    
    bit_info <- vector()
    
    for (i in 1:length(modis_QA_500)) {
      
      bit_info[i] <- as.integer(intToBits(modis_QA_500[i])[1])
      
    }
    
    mask_clouds <- t(matrix(bit_info, nrow = 118, ncol = 118))
    
    modis_bits <- raster(mask_clouds)
    
    extent(modis_bits) <- modis_ref
    
    modis_bits[modis_bits == 1 | modis_bits == 10] <- NA
    modis_bits[modis_bits == 0 | modis_bits == 11] <- 1
    
    modis_bits <- resample(modis_bits, modis_ref, method = "ngb")
    
    plot(modis_bits)
    
    modis_mask_i <- stack(modis_mask_i, modis_bits)
    
  }
  
  modis_clear <- ref_i*modis_mask_i
  
  plot(modis_clear)
  
  
  #Retrieve observation-illumination geometry from MODIS
  
  modis_theta_v_i <- stack()
  
  for (i in 0:15) {
    
    modis_theta_v <- raster(modis_list[i*12 + 2])
    
    modis_theta_v <- crop(modis_theta_v, AOI_sent)
    
    modis_theta_v[modis_theta_v == -32767] <- NA
    
    modis_theta_v <- resample(modis_theta_v, modis_ref, method = "ngb")*(pi/180)*0.01
    
    plot(modis_theta_v)
    
    modis_theta_v_i <- stack(modis_theta_v_i, modis_theta_v)
    
  }
  
  modis_theta_v_i <- modis_theta_v_i*modis_mask_i
  
  plot(modis_theta_v_i)
  
  
  modis_theta_s_i <- stack()
  
  for (i in 0:15) {
    
    modis_theta_s <- raster(modis_list[i*12 + 4])
    
    modis_theta_s <- crop(modis_theta_s, AOI_sent)
    
    modis_theta_s[modis_theta_s == -32767] <- NA
    
    modis_theta_s <- resample(modis_theta_s, modis_ref, method = "ngb")*(pi/180)*0.01
    
    plot(modis_theta_s)
    
    modis_theta_s_i <- stack(modis_theta_s_i, modis_theta_s)
    
  }
  
  modis_theta_s_i <- modis_theta_s_i*modis_mask_i
  
  plot(modis_theta_s_i)
  
  
  modis_phi_i <- stack()
  
  for (i in 0:15) {
    
    phi_v <- raster(modis_list[i*12 + 1])
    
    phi_v <- crop(phi_v, AOI_sent)*0.01
    
    phi_v[phi_v == -327.67] <- NA
    phi_v[phi_v > 180] <- NA
    phi_v[phi_v < -180] <- NA
    
    phi_v_p <- phi_v
    
    phi_v_n <- phi_v
    
    phi_v_p[phi_v <= 0] <- NA
    
    phi_v_n[phi_v > 0] <- NA
    
    phi_v_n <- phi_v_n + 360
    
    phi_v_t <- cover(phi_v_n, phi_v_p)
    
    
    phi_s <- raster(modis_list[i*12 + 3])
    
    phi_s <- crop(phi_s, AOI_sent)*0.01
    
    phi_s[phi_s == -327.67] <- NA
    phi_s[phi_s > 180] <- NA
    phi_s[phi_s < -180] <- NA
    
    phi_s_p <- phi_s
    
    phi_s_n <- phi_s
    
    phi_s_p[phi_s <= 0] <- NA
    
    phi_s_n[phi_s > 0] <- NA
    
    phi_s_n <- phi_s_n + 360
    
    phi_s_t <- cover(phi_s_n, phi_s_p)
    
    
    modis_phi <- abs(phi_v_t - phi_s_t)
    
    modis_phi_f <- modis_phi
    
    modis_phi_f[modis_phi < 180] <- NA
    
    modis_phi[modis_phi >= 180] <- NA
    
    modis_phi_f <- abs(modis_phi_f - 360)
    
    modis_phi_final <- cover(modis_phi, modis_phi_f) 
    
    
    modis_phi <- resample(modis_phi_final, modis_ref, method = "ngb")*(pi/180)
    
    plot(modis_phi)
    
    modis_phi_i <- stack(modis_phi_i, modis_phi)
    
  }
  
  modis_phi_i <- modis_phi_i*modis_mask_i
  
  plot(modis_phi_i)
  
  
  #Correct to Sentinel-2 Spectral Response

  modis_landcover <- raster("MCD12Q1 MODIS Landcover image from the same year.tif")

  modis_landcover_res <- resample(modis_landcover, modis_clear[[1]], method = "ngb")

  plot(modis_landcover_res)

  sbaf <- read.csv("Table with SBAF values.csv", sep = ";") #Remember that the SBAF changes depending on the Sentinel-2 sensor (S2A or 2SB) and that the table should follow the format used in lines 294-305

  sbaf_bands <- c("Red", "NIR", "Blue", "Green", "none", "SWIR1", "SWIR2")

  
  #You might need to add or remove lines below (between lines 294 and 305) depending on the landcover classes available in your study area
  
  modis_clear[modis_landcover_res == 1] <- modis_clear[modis_landcover_res == 1]*(sbaf$sbaf[which(sbaf$class_n == 1 & sbaf$band == sbaf_bands[w])])
  modis_clear[modis_landcover_res == 4] <- modis_clear[modis_landcover_res == 4]*(sbaf$sbaf[which(sbaf$class_n == 4 & sbaf$band == sbaf_bands[w])])
  modis_clear[modis_landcover_res == 5] <- modis_clear[modis_landcover_res == 5]*(sbaf$sbaf[which(sbaf$class_n == 5 & sbaf$band == sbaf_bands[w])])
  modis_clear[modis_landcover_res == 6] <- modis_clear[modis_landcover_res == 6]*(sbaf$sbaf[which(sbaf$class_n == 6 & sbaf$band == sbaf_bands[w])])
  modis_clear[modis_landcover_res == 7] <- modis_clear[modis_landcover_res == 7]*(sbaf$sbaf[which(sbaf$class_n == 7 & sbaf$band == sbaf_bands[w])])
  modis_clear[modis_landcover_res == 8] <- modis_clear[modis_landcover_res == 8]*(sbaf$sbaf[which(sbaf$class_n == 8 & sbaf$band == sbaf_bands[w])])
  modis_clear[modis_landcover_res == 9] <- modis_clear[modis_landcover_res == 9]*(sbaf$sbaf[which(sbaf$class_n == 9 & sbaf$band == sbaf_bands[w])])
  modis_clear[modis_landcover_res == 10] <- modis_clear[modis_landcover_res == 10]*(sbaf$sbaf[which(sbaf$class_n == 10 & sbaf$band == sbaf_bands[w])])
  modis_clear[modis_landcover_res == 11] <- modis_clear[modis_landcover_res == 11]*(sbaf$sbaf[which(sbaf$class_n == 11 & sbaf$band == sbaf_bands[w])])
  modis_clear[modis_landcover_res == 15] <- modis_clear[modis_landcover_res == 15]*(sbaf$sbaf[which(sbaf$class_n == 15 & sbaf$band == sbaf_bands[w])])
  modis_clear[modis_landcover_res == 16] <- modis_clear[modis_landcover_res == 16]*(sbaf$sbaf[which(sbaf$class_n == 16 & sbaf$band == sbaf_bands[w])])
  modis_clear[modis_landcover_res == 17] <- modis_clear[modis_landcover_res == 17]*(sbaf$sbaf[which(sbaf$class_n == 17 & sbaf$band == sbaf_bands[w])])

  modis_clear[modis_clear < 0] <- NA
  modis_clear[modis_clear > 1] <- NA

  plot(modis_clear)
  
  
  #Input percentage of pixel coverage
  
  setwd("/path to files of percentage of MODIS pixel coverage")
  
  list_perc <- list.files(getwd())
  
  target_file <- which(list.files(getwd()) == paste0("MOD09GA_PercCov_", img_date, ".tif"))
  
  list_perc_sixteen <- list_perc[(target_file-8):(target_file+7)]
  
  
  perc_stack <- stack()
  
  for (p in 1:16) {
    
    perc_cov <- raster(list_perc_sixteen[p])  
    
    perc_cov_proj <- projectRaster(perc_cov, modis_phi_i)
    
    perc_cov_crop <- resample(perc_cov_proj, modis_phi_i, method = "ngb")
    
    perc_stack <- stack(perc_stack, perc_cov_crop)
    
  }
  
  
  #Compute isotropic, volumetric, geometric, and snow kernels
  
  theta_s_vec <- as.vector(modis_theta_s_i)
  theta_v_vec <- as.vector(modis_theta_v_i)
  phi_vec <- as.vector(modis_phi_i)
  ref_vec <- as.vector(modis_clear)
  perc_vec <- as.vector(perc_stack)
  
  
  #This part of the script solves a weighted constrained linear system to retrieve the BRDF parameters
  #Note that 56169 is the number of MODIS pixels, so you will need to change this for your domain
  
  setwd("/path to ouput MODIS BRDF kernel parameters and evaluation")
  
  
  k1 <- 1.247
  k2 <- 1.186
  k3 <- 5.157
  alpha <- 0.3
  
  probs <- c(seq(0.5, 0.95, 0.05625), seq(0.90, 0.55, -0.05625))
  temp_weights <- round(qLaplace(probs, 0.5, 0.219), 2)
  
  
  k_vol <- function(theta_s, theta_v, phi) {
    
    cos_E = ((sin(theta_s))*(sin(theta_v))*(cos(phi))) + ((cos(theta_s))*(cos(theta_v)))
    
    E = (acos(cos_E))*(180/pi)
    
    
    k_vol_omega = (((((((pi/2) - (E*(pi/180)))*cos_E) + sin(E))/(cos(theta_s) + cos(theta_v))) - (pi/4))*(sin(theta_v))*(cos(theta_v)))*(sin(theta_s))*(cos(theta_s))
    
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
    
    
    k_snow_omega = (((R_0*(1-(alpha*(cos(E*(pi/180)))*(exp(-cos(E*(pi/180))))))) + (0.4076*alpha) - 1.1081)*(sin(theta_v))*(cos(theta_v)))*(sin(theta_s))*(cos(theta_s))
    
  }
  
  
  k_vol_l <- (integral3(k_vol, 0, pi/2, 0, pi/2, 0, 2*pi, reltol = 1e-2))*(2/pi)
  
  k_geo_l <- (integral3(k_geo, 0, pi/2, 0, pi/2, 0, 2*pi, reltol = 1e-2))*(2/pi)
  
  k_snow_l <- (integral3(k_snow, 0, pi/2, 0, pi/2, 0, 2*pi, reltol = 1e-2))*(2/pi)
  
  
  f_iso_final <- vector()
  f_vol_final <- vector()
  f_geo_final <- vector()
  f_snow_final <- vector()
  rmse_final <- vector()
  wod_wdr_final <- vector()
  wod_wsa_final <- vector()
  
  
  pb <- txtProgressBar(min = 1, max = length(theta_s_vec), style = 3)
  
  Sys.time()
  
  for (j in 1:length(theta_s_vec)) {
    
    
    theta_s <- theta_s_vec[c(j, j+56169, j+(2*56169), j+(3*56169), j+(4*56169), j+(5*56169), j+(6*56169),
                             j+(7*56169), j+(8*56169), j+(9*56169), j+(10*56169), j+(11*56169), j+(12*56169), 
                             j+(13*56169), j+(14*56169), j+(15*56169))]
    
    theta_v <- theta_v_vec[c(j, j+56169, j+(2*56169), j+(3*56169), j+(4*56169), j+(5*56169), j+(6*56169),
                             j+(7*56169), j+(8*56169), j+(9*56169), j+(10*56169), j+(11*56169), j+(12*56169), 
                             j+(13*56169), j+(14*56169), j+(15*56169))]
    
    phi <- phi_vec[c(j, j+56169, j+(2*56169), j+(3*56169), j+(4*56169), j+(5*56169), j+(6*56169),
                     j+(7*56169), j+(8*56169), j+(9*56169), j+(10*56169), j+(11*56169), j+(12*56169), 
                     j+(13*56169), j+(14*56169), j+(15*56169))]
    
    ref <- ref_vec[c(j, j+56169, j+(2*56169), j+(3*56169), j+(4*56169), j+(5*56169), j+(6*56169),
                     j+(7*56169), j+(8*56169), j+(9*56169), j+(10*56169), j+(11*56169), j+(12*56169), 
                     j+(13*56169), j+(14*56169), j+(15*56169))]
    
    perc <- perc_vec[c(j, j+56169, j+(2*56169), j+(3*56169), j+(4*56169), j+(5*56169), j+(6*56169),
                       j+(7*56169), j+(8*56169), j+(9*56169), j+(10*56169), j+(11*56169), j+(12*56169), 
                       j+(13*56169), j+(14*56169), j+(15*56169))]
    
    perc <- ifelse(is.na(perc), 0, perc)
    
    
    if(sum(is.na(theta_s)) >= 12 | sum(is.na(ifelse(ref < 0, NA, ref))) >= 12) {
      
      f_iso_final[j] <- NA
      f_vol_final[j] <- NA
      f_geo_final[j] <- NA
      f_snow_final[j] <- NA
      rmse_final[j] <- NA
      wod_wdr_final[j] <- NA
      wod_wsa_final[j] <- NA
      
    } else {
      
      #Calculation of BRDF reflectance
      
      k1 <- 1.247
      k2 <- 1.186
      k3 <- 5.157
      alpha <- 0.3
      
      
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
      
      
      k_snow_omega = (R_0*(1-(alpha*(cos(E*(pi/180)))*(exp(-cos(E*(pi/180))))))) + (0.4076*alpha) - 1.1081
      
      k_vol_omega = (((((pi/2) - (E*(pi/180)))*cos_E) + sin(E))/(cos(theta_s) + cos(theta_v))) - (pi/4)
      
      k_geo_omega = over - ((1/cos(theta_st)) + (1/cos(theta_vt)) - ((1/2)*(1 + cos_Et)*((1/cos(theta_vt))*(1/cos(theta_st)))))
      
      k_snow_omega <- matrix(ifelse(ref < 0, NA, k_snow_omega))
      
      k_vol_omega <- matrix(ifelse(ref < 0, NA, k_vol_omega))
      
      k_geo_omega <- matrix(ifelse(ref < 0, NA, k_geo_omega))
      
      
      cov_weights <- perc/100
      
      mean_weights <- (temp_weights + cov_weights)/2
      
      
      #Observation matrix
      
      b <- matrix(ifelse(ref < 0, NA, ref))
      
      weights_eval <- na.omit(ifelse(is.na(b), NA, mean_weights))
      
      b <- na.omit(b)
      
      #Kernel matrix
      
      A <- matrix(c(rep(1,length(b)), na.omit(k_vol_omega), na.omit(k_geo_omega), na.omit(k_snow_omega)), ncol = 4, nrow = length(b))
      
      #Weight matrix
      
      if (length(na.omit(ifelse(is.na(as.vector(b)), NA, mean_weights))) == 0) {
        
        W <- 1
        
      } else {
        
        W <- diag(na.omit(ifelse(is.na(as.vector(b)), NA, mean_weights)), nrow = length(b), ncol = length(b))
        
      }
      
      
      #Solve constrained weighted linear system
      
      x <- lsqlincon(C = W %*% A, d = W %*% matrix(b), lb = c(0,0,0,0))
      
      solution <- A %*% matrix(x)
      
      
      #Calculate RMSE
      
      rmse <- sqrt(sum(((b - solution)^2)*(weights_eval))/(nrow(b)-length(x)))
      
      
      #WoD Calculation
      
      theta_s <- pi/4 # 45 degrees
      
      theta_v <- 0 #nadir
      
      phi <- 0 #nadir
      
      
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
      
      
      k_snow_omega = (R_0*(1-(alpha*(cos(E*(pi/180)))*(exp(-cos(E*(pi/180))))))) + (0.4076*alpha) - 1.1081
      
      k_vol_omega = (((((pi/2) - (E*(pi/180)))*cos_E) + sin(E))/(cos(theta_s) + cos(theta_v))) - (pi/4)
      
      k_geo_omega = over - ((1/cos(theta_st)) + (1/cos(theta_vt)) - ((1/2)*(1 + cos_Et)*((1/cos(theta_vt))*(1/cos(theta_st)))))
      
      
      
      #WoD-WDR
      
      U_wdr <- c(1, k_vol_omega, k_geo_omega, k_snow_omega)*mean_weights[9]
      
      WoD_wdr <- t(matrix(U_wdr)) %*% ginv(t(A) %*% W %*% A) %*% matrix(U_wdr)
      
      
      #WoD-WSA
      
      U_wsa <- c(1, k_vol_l, k_geo_l, k_snow_l)*mean_weights[9]
      
      WoD_wsa <- t(matrix(U_wsa)) %*% ginv(t(A) %*% W %*% A) %*% matrix(U_wsa)
      
      
      
      f_iso_final[j] <- x[1]
      f_vol_final[j] <- x[2]
      f_geo_final[j] <- x[3]
      f_snow_final[j] <- x[4]
      rmse_final[j] <- rmse
      wod_wdr_final[j] <- WoD_wdr
      wod_wsa_final[j] <- WoD_wsa
      
    }
    
    
    Sys.sleep(0.0001)
    setTxtProgressBar(pb, j)
    
  }
  
  
  #Reconstruct rasters
  
  f_iso_matrix <- t(matrix(f_iso_final, nrow = 237, ncol = 237))
  
  f_vol_matrix <- t(matrix(f_vol_final, nrow = 237, ncol = 237))
  
  f_geo_matrix <- t(matrix(f_geo_final, nrow = 237, ncol = 237))
  
  f_snow_matrix <- t(matrix(f_snow_final, nrow = 237, ncol = 237))
  
  rmse_matrix <- t(matrix(rmse_final, nrow = 237, ncol = 237))
  
  wod_wdr_matrix <- t(matrix(wod_wdr_final, nrow = 237, ncol = 237))
  
  wod_wsa_matrix <- t(matrix(wod_wsa_final, nrow = 237, ncol = 237)) 
  
  
  f_iso_raster <- raster(f_iso_matrix)
  
  extent(f_iso_raster) <- modis_ref
  
  
  f_vol_raster <- raster(f_vol_matrix)
  
  extent(f_vol_raster) <- modis_ref
  
  
  f_geo_raster <- raster(f_geo_matrix)
  
  extent(f_geo_raster) <- modis_ref
  
  
  f_snow_raster <- raster(f_snow_matrix)
  
  extent(f_snow_raster) <- modis_ref
  
  
  rmse_raster <- raster(rmse_matrix)
  
  extent(rmse_raster) <- modis_ref
  
  
  wod_wdr_raster <- raster(wod_wdr_matrix)
  
  extent(wod_wdr_raster) <- modis_ref
  
  
  wod_wsa_raster <- raster(wod_wsa_matrix)
  
  extent(wod_wsa_raster) <- modis_ref
  
  
  writeRaster(f_iso_raster, filename = paste0("/path to ouput MODIS BRDF kernel parameters and evaluation/f_iso_", img_date, "_B", band), format = "GTiff", overwrite = T)
  
  writeRaster(f_vol_raster, filename = paste0("/path to ouput MODIS BRDF kernel parameters and evaluation/f_vol_", img_date, "_B", band), format = "GTiff", overwrite = T)
  
  writeRaster(f_geo_raster, filename = paste0("/path to ouput MODIS BRDF kernel parameters and evaluation/f_geo_", img_date, "_B", band), format = "GTiff", overwrite = T)
  
  writeRaster(f_snow_raster, filename = paste0("/path to ouput MODIS BRDF kernel parameters and evaluation/f_snow_", img_date, "_B", band), format = "GTiff", overwrite = T)
  
  writeRaster(rmse_raster, filename = paste0("/path to ouput MODIS BRDF kernel parameters and evaluation/rmse_", img_date, "_B", band), format = "GTiff", overwrite = T)
  
  writeRaster(wod_wdr_raster, filename = paste0("/path to ouput MODIS BRDF kernel parameters and evaluation/wod_wdr_", img_date, "_B", band), format = "GTiff", overwrite = T)
  
  writeRaster(wod_wsa_raster, filename = paste0("/path to ouput MODIS BRDF kernel parameters and evaluation/wod_wsa_", img_date, "_B", band), format = "GTiff", overwrite = T)
  
}

Sys.time()

