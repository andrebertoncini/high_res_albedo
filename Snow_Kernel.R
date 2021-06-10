#This code is used to retrieve MODIS BRDF using the snow kernel model, 
#which is based on Jiao et al. (2019) "Development of a snow kernel to 
#better model the anisotropic reflectance of pure snow in a kernel-driven 
#BRDF model framework" - Remote Sensing of Environment.

#Author: Andre Bertoncini

#Institution: Centre for Hydrology - University of Saskatchewan


library(raster)
library(rgdal)
library(reshape2)
#library(MODIS)


#Set target date (same as the Sentinel-2 image) 

Dm <- "22"
M <- "07"
Y <- "2019"

img_date <- paste0(Y, M, Dm)
modis_date <- paste0(Y, ".", M, ".", Dm)
Dm_j <- (as.numeric(Dm) - 32) + (floor(275*(as.numeric(M)/9))) + (2*floor(3/(as.numeric(M)+1)))
Dm_j <- ifelse(Y == "2016" | Y == "2020", Dm_j + 1, Dm_j) #for leap years: 2016 and 2020
date_j <- paste0(Y, Dm_j)

#Uncomment lines below to download MODIS data through the MODIS package.

#begin_m <- "2019.08.28"
#end_m <- "2019.09.12"


#Retrieve MODIS data

AOI_sent <- extent(399960, 509760, 5690220, 5800020) # Sentinel-2 extent in UTM

modis_ext <- extent(-120, -113, 50, 53) # MODIS extent in Geographic coordinates


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
  
  
  #Compute snow kernels
  
  
  theta_s_vec <- as.vector(modis_theta_s_i)
  theta_v_vec <- as.vector(modis_theta_v_i)
  phi_vec <- as.vector(modis_phi_i)
  ref_vec <- as.vector(modis_clear)
  
  input_data <- data.frame(theta_s_vec, theta_v_vec, phi_vec, ref_vec)
  
  write.table(input_data, file = paste0("/path to output data/SnowKernel_Data_B", w, ".csv"), sep = ";", dec = ".")
  
}


#This part of the script perform the optimization to retrieve BRDF information
#Note that 56169 is the number of MODIS pixels, so you will need to change this for your domain.

setwd("/path of output data from previous step")

for (i in c(1,2,3,4,6,7)) {

band <- i

input_data <- read.table(paste0("SnowKernel_Data_B", band, ".csv"), header = T, sep = ";", dec = ".")

theta_s_vec <- input_data$theta_s_vec
theta_v_vec <- input_data$theta_v_vec
phi_vec <- input_data$phi_vec
ref_vec <- input_data$ref_vec

k1 <- 1.247
k2 <- 1.186
k3 <- 5.157
alpha <- 0.3

f_iso_final <- vector()
f_snow_final <- vector()

pb <- txtProgressBar(min = 1, max = length(theta_s_vec), style = 3)


for (j in 1:length(theta_s_vec)) {
  
  tryCatch({
    
    
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
    
    
    cos_E = ((sin(theta_s))*(sin(theta_v))*(cos(phi))) + ((cos(theta_s))*(cos(theta_v)))
    
    E = (acos(cos_E))*(180/pi)
    
    P_E = (11.1*(exp(-0.087*(180 - E)))) + (1.1*(exp(-0.014*(180 - E))))
    
    R_0 = (k1 + (k2*(cos(theta_s) + cos(theta_v))) + (k3*(cos(theta_s)*cos(theta_v))) + P_E)/(4*(cos(theta_s) + cos(theta_v)))
    
    
    k_snow_omega = (R_0*(1-(alpha*(cos(E*(pi/180)))*(exp(-cos(E*(pi/180))))))) + (0.4076*alpha) - 1.1081
    
    
    if (sum(is.na(k_snow_omega)) == 16) {
      f_iso_final[j] <- NA
      f_snow_final[j] <- NA
    }
    
    else {
      
      
      nls_result <- nls(ref ~ y + z*k_snow_omega, start=list(y = 0.884,z = 1.176), control = list(maxiter = 1000))
      
      
      f_iso_final[j] <- unname(coef(nls_result))[1]
      f_snow_final[j] <- unname(coef(nls_result))[2]
      
    }
    
  }, error = function(e) {
    f_iso_final[j] <- NA
    f_snow_final[j] <- NA})
  
  
  Sys.sleep(0.0001)
  setTxtProgressBar(pb, j)
  
}

output_data <- data.frame(f_iso_final, f_snow_final)

write.table(output_data, file = paste0("SnowKernel_Output_B", band, ".csv"), sep = ";", dec = ".")

}


#Reconstruct raster

setwd("/path of output data from previous step")


for (i in c(1,2,3,4,6,7)) {
  
  band <- i  
  
  output_data <- read.table(paste0("SnowKernel_Output_B", band, ".csv"), header = T, sep = ";", dec = ".")
  
  f_iso_final <- output_data$f_iso_final
  
  f_snow_final <- output_data$f_snow_final
  
  
  f_iso_matrix <- t(matrix(f_iso_final, nrow = 237, ncol = 237))
  
  f_snow_matrix <- t(matrix(f_snow_final, nrow = 237, ncol = 237))
  
  
  f_iso_raster <- raster(f_iso_matrix)
  
  extent(f_iso_raster) <- modis_ref
  
  
  f_snow_raster <- raster(f_snow_matrix)
  
  extent(f_snow_raster) <- modis_ref
  
  
  plot(f_iso_raster)
  
  plot(f_snow_raster)
  
  
  writeRaster(f_iso_raster, filename = paste0("/path to output images/f_iso_", img_date, "_B", band), format = "GTiff", overwrite = T)
  
  writeRaster(f_snow_raster, filename = paste0("/path to output images/f_snow_", img_date, "_B", band), format = "GTiff", overwrite = T)
  
  
}

