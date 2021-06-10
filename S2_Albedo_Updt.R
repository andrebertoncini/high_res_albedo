#This code is used to retrieve Sentinel-2 albedo based on the MODIS BRDF information
#calculated previously in the codes Snow_Kernel and AN_Ratio. This code is based on
#Li et al. (2018) "Preliminary assessment of 20-m surface albedo retrievals from sentinel-2A 
#surface reflectance and MODIS/VIIRS surface anisotropy measures" - Remote Sensing of Environment.

#Author: Andre Bertoncini

#Institution: Centre for Hydrology - University of Saskatchewan


library(raster)
library(rgdal)
library(doParallel)
library(foreach)
library(sp)


#Set target date

Dm <- "28"
M <- "08"
Y <- "2019"

img_date <- paste0(Y, M, Dm)
modis_date <- paste0(Y, ".", M, ".", Dm)
Dm_j <- (as.numeric(Dm) - 32) + (floor(275*(as.numeric(M)/9))) + (2*floor(3/(as.numeric(M)+1)))
Dm_j <- ifelse(Y == "2016" | Y == "2020", Dm_j + 1, Dm_j) #for leap years: 2016 and 2020
date_j <- paste0(Y, Dm_j)


#Load Sentinel-2 Bands
#Red dots need to be changed

list_folders <- list.dirs("/path to Sentinel-2 folders")

folder <- grep(pattern = paste0(img_date), list_folders)

setwd(list_folders[[folder[13]]])

S2_list <- list.files(path = getwd(), pattern = paste0("T11UMT_", img_date, ".*.jp2$"), recursive = T)

AOI <- extent(421290, 496851.5, 5759734, 5800008.8)

AOI_sent <- extent(399960, 509760, 5690220, 5800020)

#Stack bands in the right order

for (i in 1:11){
  
  S2_bands <- readGDAL(S2_list[[i]])
  
  S2_bands <- raster(S2_bands)
  
  extent(S2_bands) <- AOI_sent
  projection(S2_bands) <- CRS("+proj=utm +zone=11 +datum=WGS84 +units=m
                              +no_defs +ellps=WGS84 +towgs84=0,0,0")
  
  S2_bands <- crop(S2_bands, AOI)
  
  assign(paste("S2_surfref_B", i, sep = ""), S2_bands)
  
}


S2_surfref <- stack(S2_surfref_B2, S2_surfref_B3, S2_surfref_B4, 
                    S2_surfref_B10, S2_surfref_B8, S2_surfref_B9)*0.0001

plot(S2_surfref)

rm(list=ls()[! ls() %in% c("S2_surfref_B1", "S2_surfref_B11", "S2_surfref", "img_date", "Dm_j", "Dm", "M", "Y", "j", "date_j")])


#Mask surfaces that should not be processed

mask <- S2_surfref_B11
cloudless_mask <- raster(paste0("/path to Sentinel-2 cloudless masks/Sentinel_2_CloudMask_", img_date, "_AOI.tif"))

cloudless_mask_res <- resample(cloudless_mask, S2_surfref, method = "ngb")

S2_surfref[[1]][mask == 0] <- NA
S2_surfref[[2]][mask == 0] <- NA
S2_surfref[[3]][mask == 0] <- NA
S2_surfref[[4]][mask == 0] <- NA
S2_surfref[[5]][mask == 0] <- NA
S2_surfref[[6]][mask == 0] <- NA

S2_surfref[[1]][mask == 1] <- NA
S2_surfref[[2]][mask == 1] <- NA
S2_surfref[[3]][mask == 1] <- NA
S2_surfref[[4]][mask == 1] <- NA
S2_surfref[[5]][mask == 1] <- NA
S2_surfref[[6]][mask == 1] <- NA

S2_surfref[[1]][mask == 2] <- NA
S2_surfref[[2]][mask == 2] <- NA
S2_surfref[[3]][mask == 2] <- NA
S2_surfref[[4]][mask == 2] <- NA
S2_surfref[[5]][mask == 2] <- NA
S2_surfref[[6]][mask == 2] <- NA

S2_surfref[[1]][cloudless_mask_res == 1] <- NA
S2_surfref[[2]][cloudless_mask_res == 1] <- NA
S2_surfref[[3]][cloudless_mask_res == 1] <- NA
S2_surfref[[4]][cloudless_mask_res == 1] <- NA
S2_surfref[[5]][cloudless_mask_res == 1] <- NA
S2_surfref[[6]][cloudless_mask_res == 1] <- NA

plot(S2_surfref)


for (j in 6:20) {
  
  #j = 16
  
  #Perform kmeans classification based on Sentinel-2 surface reflectance bands
  
  setwd("/path to albedo outputs")
  
  
  nr <- getValues(S2_surfref)
  
  set.seed(99)
  
  kmncluster <- kmeans(na.omit(nr), centers = j, iter.max = 500, nstart = 5, algorithm="Lloyd")
  
  
  npixel <- length(nr[,1])
  
  cluster_num <- as.vector(nr)
  
  cluster_num[which(!is.na(as.vector(nr)))] <- kmncluster$cluster
  
  cluster_b1 <- cluster_num[1:npixel]
  cluster_b2 <- cluster_num[(npixel + 1):(npixel*2)]
  cluster_b3 <- cluster_num[((npixel*2) + 1):(npixel*3)]
  cluster_b4 <- cluster_num[((npixel*3) + 1):(npixel*4)]
  cluster_b5 <- cluster_num[((npixel*4) + 1):(npixel*5)]
  cluster_b6 <- cluster_num[((npixel*5) + 1):(npixel*6)]
  
  
  cluster_band <- round(((cluster_b1 + cluster_b2 + cluster_b3
                          + cluster_b4 + cluster_b5 + cluster_b6)/6), digits = 0)
  
  S2_cluster <- setValues(S2_surfref[[1]], cluster_band)
  
  plot(S2_cluster)
  
  
  #Select MODIS pixels that are covered more than 60% by the same class
  
  modis_grid <- shapefile("shapefile of MODIS grid.shp")
  
{
    
    cluster_1 <- S2_cluster
    cluster_1[cluster_1 != 1] <- NA
    
    cluster_2 <- S2_cluster
    cluster_2[cluster_2 != 2] <- NA
    
    cluster_3 <- S2_cluster
    cluster_3[cluster_3 != 3] <- NA
    
    cluster_4 <- S2_cluster
    cluster_4[cluster_4 != 4] <- NA
    
    cluster_5 <- S2_cluster
    cluster_5[cluster_5 != 5] <- NA
    
    cluster_6 <- S2_cluster
    cluster_6[cluster_6 != 6] <- NA
    
    cluster_7 <- S2_cluster
    cluster_7[cluster_7 != 7] <- NA
    
    cluster_8 <- S2_cluster
    cluster_8[cluster_8 != 8] <- NA
    
    cluster_9 <- S2_cluster
    cluster_9[cluster_9 != 9] <- NA
    
    cluster_10 <- S2_cluster
    cluster_10[cluster_10 != 10] <- NA
    
    cluster_11 <- S2_cluster
    cluster_11[cluster_11 != 11] <- NA
    
    cluster_12 <- S2_cluster
    cluster_12[cluster_12 != 12] <- NA
    
    cluster_13 <- S2_cluster
    cluster_13[cluster_13 != 13] <- NA
    
    cluster_14 <- S2_cluster
    cluster_14[cluster_14 != 14] <- NA
    
    cluster_15 <- S2_cluster
    cluster_15[cluster_15 != 15] <- NA
    
    cluster_16 <- S2_cluster
    cluster_16[cluster_16 != 16] <- NA
    
    cluster_17 <- S2_cluster
    cluster_17[cluster_17 != 17] <- NA
    
    cluster_18 <- S2_cluster
    cluster_18[cluster_18 != 18] <- NA
    
    cluster_19 <- S2_cluster
    cluster_19[cluster_19 != 19] <- NA
    
    cluster_20 <- S2_cluster
    cluster_20[cluster_20 != 20] <- NA
    
    cluster_list <- list(cluster_1, cluster_2, cluster_3, cluster_4, cluster_5, cluster_6, cluster_7, cluster_8,
                         cluster_9, cluster_10, cluster_11, cluster_12, cluster_13, cluster_14, cluster_15, cluster_16,
                         cluster_17, cluster_18, cluster_19, cluster_20)
    
  }


UseCores <- detectCores() -1

cl <- makeCluster(UseCores)
registerDoParallel(cl)

modis_stats <- foreach(k = 1:j, .combine = "cbind") %dopar% {
  
  library(raster)
  
  extract(cluster_list[[k]], modis_grid, fun=function(x, ...) length(na.omit(x)))
  
}

stopCluster(cl)


modis_stats_df <- data.frame(modis_stats)

modis_class_60 <- matrix(nr = nrow(modis_stats_df), nc = j)

for (v in 1:j) {
  
  modis_class_60[,v] <- ifelse(modis_stats_df[,v] > 375, v, NA)
  
}

modis_grid@data$Id <- 1:12000

modis_one_col <- as.data.frame(rowMeans(modis_class_60, na.rm = T, dims = 1))

colnames(modis_one_col) <- "maj_cluster"

modis_one_col$maj_cluster <- ifelse(is.nan(modis_one_col$maj_cluster), NA, modis_one_col$maj_cluster)  

modis_spatial <- addAttrToGeom(x = modis_grid, y = modis_one_col, match.ID = F)

plot(modis_spatial, col=modis_spatial@data$maj_cluster, border="light gray", lwd=0.5)

modis_aggr <- aggregate(modis_spatial, by = "maj_cluster")

modis_aggr <- modis_aggr[-c(j+1), 1]

plot(modis_aggr, col=modis_aggr@data$maj_cluster, border="light gray", lwd=0.5)


#Calculate Albedo-to-Nadir ratios


setwd("/path to ouputs of the AN_ratio script")

bsa_list <- list.files(getwd(), pattern = paste0("a_bsa_", img_date, ".*.tif$"), recursive = T)

bsa_stack <- stack(bsa_list)

bsa_stack_1 <- stack(bsa_stack[[3]], bsa_stack[[4]], bsa_stack[[1]], bsa_stack[[2]], bsa_stack[[5]], bsa_stack[[6]])


wsa_list <- list.files(getwd(), pattern = paste0("a_wsa_", img_date, ".*.tif$"), recursive = T)

wsa_stack <- stack(wsa_list)

wsa_stack_1 <- stack(wsa_stack[[3]], wsa_stack[[4]], wsa_stack[[1]], wsa_stack[[2]], wsa_stack[[5]], wsa_stack[[6]])



for (l in 1:6) {
  
  band <- l
  
  modis_an_bsa <- bsa_stack_1[[l]]
  
  modis_extract <- extract(modis_an_bsa, modis_aggr, fun = mean, na.rm = T)
  
  ratio_table <- data.frame(modis_aggr$maj_cluster, modis_extract)
  
  names(ratio_table) <- c("Cluster", "Ratio")
  
  
  write.table(ratio_table, file = paste0("/path to Albedo outputs/Ratios_", img_date, "_bsa_B", band, ".csv"), sep = ";", dec = ".")
  
  
  
  modis_an_wsa <- wsa_stack_1[[l]]
  
  modis_extract <- extract(modis_an_wsa, modis_aggr, fun = mean, na.rm = T)
  
  ratio_table <- data.frame(modis_aggr$maj_cluster, modis_extract)
  
  names(ratio_table) <- c("Cluster", "Ratio")
  
  
  write.table(ratio_table, file = paste0("/path to Albedo outputs/Ratios_", img_date, "_wsa_B", band, ".csv"), sep = ";", dec = ".")
  
}


#Apply ratios to Sentinel-2 bands from 6 to 20 cluster classification classes. This is done for
#sensitivity analysis purposes.

setwd("/path to Albedo outputs")

rm(list=ls()[! ls() %in% c("S2_cluster", "S2_surfref", "img_date", "Dm_j", "Dm", "M", "Y", "j", "date_j")])

AOI2 <- extent(421290, 496851.5, 5759734, 5800008.8)

S2_surfref <- crop(S2_surfref, AOI2)
S2_cluster <- crop(S2_cluster, AOI2)


albedo_function_c6 <- function(S2_cluster, S2_surfref) {       
  
  for (i in 1:6){
    
    band = i
    
    ratio_table <- read.csv(paste0("Ratios_", img_date, "_bsa_B", band, ".csv"), sep = ";", header = T)
    
    cluster_final <- as.matrix(S2_cluster)
    albedo_dir <- as.matrix(S2_surfref[[i]])
    
    Albedo_Hemispherical <- ifelse(cluster_final == 1, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 1]),
                                   ifelse(cluster_final == 2, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 2]),
                                          ifelse(cluster_final == 3, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 3]),
                                                 ifelse(cluster_final == 4, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 4]),
                                                        ifelse(cluster_final == 5, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 5]),
                                                               ifelse(cluster_final == 6, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 6]),
                                                                      NA))))))
    
    
    albedo_bsa <- S2_surfref[[i]]
    
    albedo_bsa[] <- Albedo_Hemispherical
    
    writeRaster(albedo_bsa, filename = paste0("albedo_bsa", img_date, "_B", band, ".tif"), format = "GTiff", overwrite = T)
    
    
  }
  
  
  plot(albedo_bsa)
  
  
  for (i in 1:6){
    
    band = i
    
    ratio_table <- read.csv(paste0("Ratios_", img_date, "_wsa_B", band, ".csv"), sep = ";", header = T)
    
    cluster_final <- as.matrix(S2_cluster)
    albedo_dir <- as.matrix(S2_surfref[[i]])
    
    Albedo_Hemispherical <- ifelse(cluster_final == 1, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 1]),
                                   ifelse(cluster_final == 2, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 2]),
                                          ifelse(cluster_final == 3, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 3]),
                                                 ifelse(cluster_final == 4, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 4]),
                                                        ifelse(cluster_final == 5, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 5]),
                                                               ifelse(cluster_final == 6, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 6]),
                                                                      NA))))))
    
    
    albedo_bsa <- S2_surfref[[i]]
    
    albedo_bsa[] <- Albedo_Hemispherical
    
    writeRaster(albedo_bsa, filename = paste0("albedo_wsa", img_date, "_B", band, ".tif"), format = "GTiff", overwrite = T)
    
    
  }
  
  
  plot(albedo_bsa)
  
}
albedo_function_c7 <- function(S2_cluster, S2_surfref) {
  
  for (i in 1:6){
    
    band = i
    
    ratio_table <- read.csv(paste0("Ratios_", img_date, "_bsa_B", band, ".csv"), sep = ";", header = T)
    
    cluster_final <- as.matrix(S2_cluster)
    albedo_dir <- as.matrix(S2_surfref[[i]])
    
    Albedo_Hemispherical <- ifelse(cluster_final == 1, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 1]),
                                   ifelse(cluster_final == 2, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 2]),
                                          ifelse(cluster_final == 3, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 3]),
                                                 ifelse(cluster_final == 4, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 4]),
                                                        ifelse(cluster_final == 5, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 5]),
                                                               ifelse(cluster_final == 6, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 6]),
                                                                      ifelse(cluster_final == 7, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 7]),
                                                                             NA)))))))
    
    
    albedo_bsa <- S2_surfref[[i]]
    
    albedo_bsa[] <- Albedo_Hemispherical
    
    writeRaster(albedo_bsa, filename = paste0("albedo_bsa", img_date, "_B", band, ".tif"), format = "GTiff", overwrite = T)
    
    
  }
  
  
  plot(albedo_bsa)
  
  
  
  for (i in 1:6){
    
    band = i
    
    ratio_table <- read.csv(paste0("Ratios_", img_date, "_wsa_B", band, ".csv"), sep = ";", header = T)
    
    cluster_final <- as.matrix(S2_cluster)
    albedo_dir <- as.matrix(S2_surfref[[i]])
    
    Albedo_Hemispherical <- ifelse(cluster_final == 1, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 1]),
                                   ifelse(cluster_final == 2, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 2]),
                                          ifelse(cluster_final == 3, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 3]),
                                                 ifelse(cluster_final == 4, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 4]),
                                                        ifelse(cluster_final == 5, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 5]),
                                                               ifelse(cluster_final == 6, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 6]),
                                                                      ifelse(cluster_final == 7, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 7]),
                                                                             NA)))))))
    
    
    albedo_bsa <- S2_surfref[[i]]
    
    albedo_bsa[] <- Albedo_Hemispherical
    
    writeRaster(albedo_bsa, filename = paste0("albedo_wsa", img_date, "_B", band, ".tif"), format = "GTiff", overwrite = T)
    
    
  }
  
  
  plot(albedo_bsa)
  
} 
albedo_function_c8 <- function(S2_cluster, S2_surfref) {
  
  for (i in 1:6){
    
    band = i
    
    ratio_table <- read.csv(paste0("Ratios_", img_date, "_bsa_B", band, ".csv"), sep = ";", header = T)
    
    cluster_final <- as.matrix(S2_cluster)
    albedo_dir <- as.matrix(S2_surfref[[i]])
    
    Albedo_Hemispherical <- ifelse(cluster_final == 1, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 1]),
                                   ifelse(cluster_final == 2, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 2]),
                                          ifelse(cluster_final == 3, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 3]),
                                                 ifelse(cluster_final == 4, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 4]),
                                                        ifelse(cluster_final == 5, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 5]),
                                                               ifelse(cluster_final == 6, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 6]),
                                                                      ifelse(cluster_final == 7, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 7]),
                                                                             ifelse(cluster_final == 8, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 8]),
                                                                                    NA))))))))
    
    
    albedo_bsa <- S2_surfref[[i]]
    
    albedo_bsa[] <- Albedo_Hemispherical
    
    writeRaster(albedo_bsa, filename = paste0("albedo_bsa", img_date, "_B", band, ".tif"), format = "GTiff", overwrite = T)
    
    
  }
  
  
  plot(albedo_bsa)
  
  
  
  for (i in 1:6){
    
    band = i
    
    ratio_table <- read.csv(paste0("Ratios_", img_date, "_wsa_B", band, ".csv"), sep = ";", header = T)
    
    cluster_final <- as.matrix(S2_cluster)
    albedo_dir <- as.matrix(S2_surfref[[i]])
    
    Albedo_Hemispherical <- ifelse(cluster_final == 1, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 1]),
                                   ifelse(cluster_final == 2, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 2]),
                                          ifelse(cluster_final == 3, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 3]),
                                                 ifelse(cluster_final == 4, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 4]),
                                                        ifelse(cluster_final == 5, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 5]),
                                                               ifelse(cluster_final == 6, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 6]),
                                                                      ifelse(cluster_final == 7, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 7]),
                                                                             ifelse(cluster_final == 8, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 8]),
                                                                                    NA))))))))
    
    
    albedo_bsa <- S2_surfref[[i]]
    
    albedo_bsa[] <- Albedo_Hemispherical
    
    writeRaster(albedo_bsa, filename = paste0("albedo_wsa", img_date, "_B", band, ".tif"), format = "GTiff", overwrite = T)
    
    
  }
  
  
  plot(albedo_bsa)  
  
} 
albedo_function_c9 <- function(S2_cluster, S2_surfref) {
  
  for (i in 1:6){
    
    band = i
    
    ratio_table <- read.csv(paste0("Ratios_", img_date, "_bsa_B", band, ".csv"), sep = ";", header = T)
    
    cluster_final <- as.matrix(S2_cluster)
    albedo_dir <- as.matrix(S2_surfref[[i]])
    
    Albedo_Hemispherical <- ifelse(cluster_final == 1, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 1]),
                                   ifelse(cluster_final == 2, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 2]),
                                          ifelse(cluster_final == 3, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 3]),
                                                 ifelse(cluster_final == 4, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 4]),
                                                        ifelse(cluster_final == 5, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 5]),
                                                               ifelse(cluster_final == 6, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 6]),
                                                                      ifelse(cluster_final == 7, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 7]),
                                                                             ifelse(cluster_final == 8, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 8]),
                                                                                    ifelse(cluster_final == 9, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 9]),
                                                                                           NA)))))))))
    
    
    albedo_bsa <- S2_surfref[[i]]
    
    albedo_bsa[] <- Albedo_Hemispherical
    
    writeRaster(albedo_bsa, filename = paste0("albedo_bsa", img_date, "_B", band, ".tif"), format = "GTiff", overwrite = T)
    
    
  }
  
  
  plot(albedo_bsa)
  
  
  
  for (i in 1:6){
    
    band = i
    
    ratio_table <- read.csv(paste0("Ratios_", img_date, "_wsa_B", band, ".csv"), sep = ";", header = T)
    
    cluster_final <- as.matrix(S2_cluster)
    albedo_dir <- as.matrix(S2_surfref[[i]])
    
    Albedo_Hemispherical <- ifelse(cluster_final == 1, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 1]),
                                   ifelse(cluster_final == 2, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 2]),
                                          ifelse(cluster_final == 3, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 3]),
                                                 ifelse(cluster_final == 4, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 4]),
                                                        ifelse(cluster_final == 5, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 5]),
                                                               ifelse(cluster_final == 6, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 6]),
                                                                      ifelse(cluster_final == 7, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 7]),
                                                                             ifelse(cluster_final == 8, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 8]),
                                                                                    ifelse(cluster_final == 9, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 9]),
                                                                                           NA)))))))))
    
    
    albedo_bsa <- S2_surfref[[i]]
    
    albedo_bsa[] <- Albedo_Hemispherical
    
    writeRaster(albedo_bsa, filename = paste0("albedo_wsa", img_date, "_B", band, ".tif"), format = "GTiff", overwrite = T)
    
    
  }
  
  
  plot(albedo_bsa)  
  
}
albedo_function_c10 <- function(S2_cluster, S2_surfref) {
  
  for (i in 1:6){
    
    band = i
    
    ratio_table <- read.csv(paste0("Ratios_", img_date, "_bsa_B", band, ".csv"), sep = ";", header = T)
    
    cluster_final <- as.matrix(S2_cluster)
    albedo_dir <- as.matrix(S2_surfref[[i]])
    
    Albedo_Hemispherical <- ifelse(cluster_final == 1, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 1]),
                                   ifelse(cluster_final == 2, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 2]),
                                          ifelse(cluster_final == 3, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 3]),
                                                 ifelse(cluster_final == 4, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 4]),
                                                        ifelse(cluster_final == 5, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 5]),
                                                               ifelse(cluster_final == 6, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 6]),
                                                                      ifelse(cluster_final == 7, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 7]),
                                                                             ifelse(cluster_final == 8, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 8]),
                                                                                    ifelse(cluster_final == 9, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 9]),
                                                                                           ifelse(cluster_final == 10, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 10]),
                                                                                                  NA))))))))))
    
    
    albedo_bsa <- S2_surfref[[i]]
    
    albedo_bsa[] <- Albedo_Hemispherical
    
    writeRaster(albedo_bsa, filename = paste0("albedo_bsa", img_date, "_B", band, ".tif"), format = "GTiff", overwrite = T)
    
    
  }
  
  
  plot(albedo_bsa)
  
  
  
  for (i in 1:6){
    
    band = i
    
    ratio_table <- read.csv(paste0("Ratios_", img_date, "_wsa_B", band, ".csv"), sep = ";", header = T)
    
    cluster_final <- as.matrix(S2_cluster)
    albedo_dir <- as.matrix(S2_surfref[[i]])
    
    Albedo_Hemispherical <- ifelse(cluster_final == 1, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 1]),
                                   ifelse(cluster_final == 2, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 2]),
                                          ifelse(cluster_final == 3, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 3]),
                                                 ifelse(cluster_final == 4, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 4]),
                                                        ifelse(cluster_final == 5, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 5]),
                                                               ifelse(cluster_final == 6, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 6]),
                                                                      ifelse(cluster_final == 7, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 7]),
                                                                             ifelse(cluster_final == 8, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 8]),
                                                                                    ifelse(cluster_final == 9, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 9]),
                                                                                           ifelse(cluster_final == 10, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 10]),
                                                                                                  NA))))))))))
    
    
    albedo_bsa <- S2_surfref[[i]]
    
    albedo_bsa[] <- Albedo_Hemispherical
    
    writeRaster(albedo_bsa, filename = paste0("albedo_wsa", img_date, "_B", band, ".tif"), format = "GTiff", overwrite = T)
    
    
  }
  
  
  plot(albedo_bsa)  
  
  
}
albedo_function_c11 <- function(S2_cluster, S2_surfref) {
  
  for (i in 1:6){
    
    band = i
    
    ratio_table <- read.csv(paste0("Ratios_", img_date, "_bsa_B", band, ".csv"), sep = ";", header = T)
    
    cluster_final <- as.matrix(S2_cluster)
    albedo_dir <- as.matrix(S2_surfref[[i]])
    
    Albedo_Hemispherical <- ifelse(cluster_final == 1, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 1]),
                                   ifelse(cluster_final == 2, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 2]),
                                          ifelse(cluster_final == 3, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 3]),
                                                 ifelse(cluster_final == 4, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 4]),
                                                        ifelse(cluster_final == 5, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 5]),
                                                               ifelse(cluster_final == 6, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 6]),
                                                                      ifelse(cluster_final == 7, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 7]),
                                                                             ifelse(cluster_final == 8, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 8]),
                                                                                    ifelse(cluster_final == 9, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 9]),
                                                                                           ifelse(cluster_final == 10, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 10]),
                                                                                                  ifelse(cluster_final == 11, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 11]),
                                                                                                         NA)))))))))))
    
    
    albedo_bsa <- S2_surfref[[i]]
    
    albedo_bsa[] <- Albedo_Hemispherical
    
    writeRaster(albedo_bsa, filename = paste0("albedo_bsa", img_date, "_B", band, ".tif"), format = "GTiff", overwrite = T)
    
    
  }
  
  
  plot(albedo_bsa)
  
  
  
  for (i in 1:6){
    
    band = i
    
    ratio_table <- read.csv(paste0("Ratios_", img_date, "_wsa_B", band, ".csv"), sep = ";", header = T)
    
    cluster_final <- as.matrix(S2_cluster)
    albedo_dir <- as.matrix(S2_surfref[[i]])
    
    Albedo_Hemispherical <- ifelse(cluster_final == 1, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 1]),
                                   ifelse(cluster_final == 2, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 2]),
                                          ifelse(cluster_final == 3, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 3]),
                                                 ifelse(cluster_final == 4, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 4]),
                                                        ifelse(cluster_final == 5, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 5]),
                                                               ifelse(cluster_final == 6, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 6]),
                                                                      ifelse(cluster_final == 7, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 7]),
                                                                             ifelse(cluster_final == 8, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 8]),
                                                                                    ifelse(cluster_final == 9, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 9]),
                                                                                           ifelse(cluster_final == 10, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 10]),
                                                                                                  ifelse(cluster_final == 11, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 11]),
                                                                                                         NA)))))))))))
    
    
    albedo_bsa <- S2_surfref[[i]]
    
    albedo_bsa[] <- Albedo_Hemispherical
    
    writeRaster(albedo_bsa, filename = paste0("albedo_wsa", img_date, "_B", band, ".tif"), format = "GTiff", overwrite = T)
    
    
  }
  
  
  plot(albedo_bsa)  
  
  
}
albedo_function_c12 <- function(S2_cluster, S2_surfref) {
  
  for (i in 1:6){
    
    band = i
    
    ratio_table <- read.csv(paste0("Ratios_", img_date, "_bsa_B", band, ".csv"), sep = ";", header = T)
    
    cluster_final <- as.matrix(S2_cluster)
    albedo_dir <- as.matrix(S2_surfref[[i]])
    
    Albedo_Hemispherical <- ifelse(cluster_final == 1, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 1]),
                                   ifelse(cluster_final == 2, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 2]),
                                          ifelse(cluster_final == 3, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 3]),
                                                 ifelse(cluster_final == 4, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 4]),
                                                        ifelse(cluster_final == 5, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 5]),
                                                               ifelse(cluster_final == 6, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 6]),
                                                                      ifelse(cluster_final == 7, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 7]),
                                                                             ifelse(cluster_final == 8, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 8]),
                                                                                    ifelse(cluster_final == 9, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 9]),
                                                                                           ifelse(cluster_final == 10, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 10]),
                                                                                                  ifelse(cluster_final == 11, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 11]),
                                                                                                         ifelse(cluster_final == 12, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 12]),
                                                                                                                NA))))))))))))
    
    
    albedo_bsa <- S2_surfref[[i]]
    
    albedo_bsa[] <- Albedo_Hemispherical
    
    writeRaster(albedo_bsa, filename = paste0("albedo_bsa", img_date, "_B", band, ".tif"), format = "GTiff", overwrite = T)
    
    
  }
  
  
  plot(albedo_bsa)
  
  
  
  for (i in 1:6){
    
    band = i
    
    ratio_table <- read.csv(paste0("Ratios_", img_date, "_wsa_B", band, ".csv"), sep = ";", header = T)
    
    cluster_final <- as.matrix(S2_cluster)
    albedo_dir <- as.matrix(S2_surfref[[i]])
    
    Albedo_Hemispherical <- ifelse(cluster_final == 1, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 1]),
                                   ifelse(cluster_final == 2, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 2]),
                                          ifelse(cluster_final == 3, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 3]),
                                                 ifelse(cluster_final == 4, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 4]),
                                                        ifelse(cluster_final == 5, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 5]),
                                                               ifelse(cluster_final == 6, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 6]),
                                                                      ifelse(cluster_final == 7, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 7]),
                                                                             ifelse(cluster_final == 8, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 8]),
                                                                                    ifelse(cluster_final == 9, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 9]),
                                                                                           ifelse(cluster_final == 10, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 10]),
                                                                                                  ifelse(cluster_final == 11, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 11]),
                                                                                                         ifelse(cluster_final == 12, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 12]),
                                                                                                                NA))))))))))))
    
    
    albedo_bsa <- S2_surfref[[i]]
    
    albedo_bsa[] <- Albedo_Hemispherical
    
    writeRaster(albedo_bsa, filename = paste0("albedo_wsa", img_date, "_B", band, ".tif"), format = "GTiff", overwrite = T)
    
    
  }
  
  
  plot(albedo_bsa)  
  
  
}
albedo_function_c13 <- function(S2_cluster, S2_surfref) {
  
  for (i in 1:6){
    
    band = i
    
    ratio_table <- read.csv(paste0("Ratios_", img_date, "_bsa_B", band, ".csv"), sep = ";", header = T)
    
    cluster_final <- as.matrix(S2_cluster)
    albedo_dir <- as.matrix(S2_surfref[[i]])
    
    Albedo_Hemispherical <- ifelse(cluster_final == 1, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 1]),
                                   ifelse(cluster_final == 2, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 2]),
                                          ifelse(cluster_final == 3, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 3]),
                                                 ifelse(cluster_final == 4, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 4]),
                                                        ifelse(cluster_final == 5, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 5]),
                                                               ifelse(cluster_final == 6, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 6]),
                                                                      ifelse(cluster_final == 7, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 7]),
                                                                             ifelse(cluster_final == 8, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 8]),
                                                                                    ifelse(cluster_final == 9, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 9]),
                                                                                           ifelse(cluster_final == 10, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 10]),
                                                                                                  ifelse(cluster_final == 11, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 11]),
                                                                                                         ifelse(cluster_final == 12, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 12]),
                                                                                                                ifelse(cluster_final == 13, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 13]),
                                                                                                                       NA)))))))))))))
    
    
    albedo_bsa <- S2_surfref[[i]]
    
    albedo_bsa[] <- Albedo_Hemispherical
    
    writeRaster(albedo_bsa, filename = paste0("albedo_bsa", img_date, "_B", band, ".tif"), format = "GTiff", overwrite = T)
    
    
  }
  
  
  plot(albedo_bsa)
  
  
  
  for (i in 1:6){
    
    band = i
    
    ratio_table <- read.csv(paste0("Ratios_", img_date, "_wsa_B", band, ".csv"), sep = ";", header = T)
    
    cluster_final <- as.matrix(S2_cluster)
    albedo_dir <- as.matrix(S2_surfref[[i]])
    
    Albedo_Hemispherical <- ifelse(cluster_final == 1, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 1]),
                                   ifelse(cluster_final == 2, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 2]),
                                          ifelse(cluster_final == 3, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 3]),
                                                 ifelse(cluster_final == 4, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 4]),
                                                        ifelse(cluster_final == 5, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 5]),
                                                               ifelse(cluster_final == 6, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 6]),
                                                                      ifelse(cluster_final == 7, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 7]),
                                                                             ifelse(cluster_final == 8, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 8]),
                                                                                    ifelse(cluster_final == 9, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 9]),
                                                                                           ifelse(cluster_final == 10, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 10]),
                                                                                                  ifelse(cluster_final == 11, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 11]),
                                                                                                         ifelse(cluster_final == 12, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 12]),
                                                                                                                ifelse(cluster_final == 13, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 13]),
                                                                                                                       NA)))))))))))))
    
    
    albedo_bsa <- S2_surfref[[i]]
    
    albedo_bsa[] <- Albedo_Hemispherical
    
    writeRaster(albedo_bsa, filename = paste0("albedo_wsa", img_date, "_B", band, ".tif"), format = "GTiff", overwrite = T)
    
    
  }
  
  
  plot(albedo_bsa)
  
  
}
albedo_function_c14 <- function(S2_cluster, S2_surfref) {
  
  for (i in 1:6){
    
    band = i
    
    ratio_table <- read.csv(paste0("Ratios_", img_date, "_bsa_B", band, ".csv"), sep = ";", header = T)
    
    ratio_table <- ratio_table[!is.na(ratio_table$Cluster), ]
    
    cluster_final <- as.matrix(S2_cluster)
    albedo_dir <- as.matrix(S2_surfref[[i]])
    
    Albedo_Hemispherical <- ifelse(cluster_final == 1, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 1]),
                                   ifelse(cluster_final == 2, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 2]),
                                          ifelse(cluster_final == 3, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 3]),
                                                 ifelse(cluster_final == 4, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 4]),
                                                        ifelse(cluster_final == 5, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 5]),
                                                               ifelse(cluster_final == 6, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 6]),
                                                                      ifelse(cluster_final == 7, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 7]),
                                                                             ifelse(cluster_final == 8, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 8]),
                                                                                    ifelse(cluster_final == 9, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 9]),
                                                                                           ifelse(cluster_final == 10, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 10]),
                                                                                                  ifelse(cluster_final == 11, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 11]),
                                                                                                         ifelse(cluster_final == 12, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 12]),
                                                                                                                ifelse(cluster_final == 13, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 13]),
                                                                                                                       ifelse(cluster_final == 14, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 14]),
                                                                                                                              NA))))))))))))))
    
    
    albedo_bsa <- S2_surfref[[i]]
    
    albedo_bsa[] <- Albedo_Hemispherical
    
    writeRaster(albedo_bsa, filename = paste0("albedo_bsa", img_date, "_B", band, ".tif"), format = "GTiff", overwrite = T)
    
    
  }
  
  
  plot(albedo_bsa)
  
  
  
  for (i in 1:6){
    
    band = i
    
    ratio_table <- read.csv(paste0("Ratios_", img_date, "_wsa_B", band, ".csv"), sep = ";", header = T)
    
    ratio_table <- ratio_table[!is.na(ratio_table$Cluster), ]
    
    cluster_final <- as.matrix(S2_cluster)
    albedo_dir <- as.matrix(S2_surfref[[i]])
    
    Albedo_Hemispherical <- ifelse(cluster_final == 1, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 1]),
                                   ifelse(cluster_final == 2, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 2]),
                                          ifelse(cluster_final == 3, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 3]),
                                                 ifelse(cluster_final == 4, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 4]),
                                                        ifelse(cluster_final == 5, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 5]),
                                                               ifelse(cluster_final == 6, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 6]),
                                                                      ifelse(cluster_final == 7, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 7]),
                                                                             ifelse(cluster_final == 8, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 8]),
                                                                                    ifelse(cluster_final == 9, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 9]),
                                                                                           ifelse(cluster_final == 10, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 10]),
                                                                                                  ifelse(cluster_final == 11, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 11]),
                                                                                                         ifelse(cluster_final == 12, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 12]),
                                                                                                                ifelse(cluster_final == 13, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 13]),
                                                                                                                       ifelse(cluster_final == 14, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 14]),
                                                                                                                              NA))))))))))))))
    
    
    albedo_bsa <- S2_surfref[[i]]
    
    albedo_bsa[] <- Albedo_Hemispherical
    
    writeRaster(albedo_bsa, filename = paste0("albedo_wsa", img_date, "_B", band, ".tif"), format = "GTiff", overwrite = T)
    
    
  }
  
  
  plot(albedo_bsa)
  
  
}
albedo_function_c15 <- function(S2_cluster, S2_surfref) {
  
  for (i in 1:6){
    
    band = i
    
    ratio_table <- read.csv(paste0("Ratios_", img_date, "_bsa_B", band, ".csv"), sep = ";", header = T)
    
    ratio_table <- ratio_table[!is.na(ratio_table$Cluster), ]
    
    cluster_final <- as.matrix(S2_cluster)
    albedo_dir <- as.matrix(S2_surfref[[i]])
    
    Albedo_Hemispherical <- ifelse(cluster_final == 1, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 1]),
                                   ifelse(cluster_final == 2, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 2]),
                                          ifelse(cluster_final == 3, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 3]),
                                                 ifelse(cluster_final == 4, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 4]),
                                                        ifelse(cluster_final == 5, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 5]),
                                                               ifelse(cluster_final == 6, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 6]),
                                                                      ifelse(cluster_final == 7, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 7]),
                                                                             ifelse(cluster_final == 8, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 8]),
                                                                                    ifelse(cluster_final == 9, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 9]),
                                                                                           ifelse(cluster_final == 10, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 10]),
                                                                                                  ifelse(cluster_final == 11, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 11]),
                                                                                                         ifelse(cluster_final == 12, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 12]),
                                                                                                                ifelse(cluster_final == 13, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 13]),
                                                                                                                       ifelse(cluster_final == 14, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 14]),
                                                                                                                              ifelse(cluster_final == 15, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 15]),
                                                                                                                                     NA)))))))))))))))
    
    
    albedo_bsa <- S2_surfref[[i]]
    
    albedo_bsa[] <- Albedo_Hemispherical
    
    writeRaster(albedo_bsa, filename = paste0("albedo_bsa", img_date, "_B", band, ".tif"), format = "GTiff", overwrite = T)
    
    
  }
  
  
  plot(albedo_bsa)
  
  
  
  for (i in 1:6){
    
    band = i
    
    ratio_table <- read.csv(paste0("Ratios_", img_date, "_wsa_B", band, ".csv"), sep = ";", header = T)
    
    ratio_table <- ratio_table[!is.na(ratio_table$Cluster), ]
    
    cluster_final <- as.matrix(S2_cluster)
    albedo_dir <- as.matrix(S2_surfref[[i]])
    
    Albedo_Hemispherical <- ifelse(cluster_final == 1, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 1]),
                                   ifelse(cluster_final == 2, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 2]),
                                          ifelse(cluster_final == 3, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 3]),
                                                 ifelse(cluster_final == 4, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 4]),
                                                        ifelse(cluster_final == 5, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 5]),
                                                               ifelse(cluster_final == 6, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 6]),
                                                                      ifelse(cluster_final == 7, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 7]),
                                                                             ifelse(cluster_final == 8, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 8]),
                                                                                    ifelse(cluster_final == 9, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 9]),
                                                                                           ifelse(cluster_final == 10, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 10]),
                                                                                                  ifelse(cluster_final == 11, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 11]),
                                                                                                         ifelse(cluster_final == 12, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 12]),
                                                                                                                ifelse(cluster_final == 13, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 13]),
                                                                                                                       ifelse(cluster_final == 14, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 14]),
                                                                                                                              ifelse(cluster_final == 15, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 15]),
                                                                                                                                     NA)))))))))))))))
    
    
    albedo_bsa <- S2_surfref[[i]]
    
    albedo_bsa[] <- Albedo_Hemispherical
    
    writeRaster(albedo_bsa, filename = paste0("albedo_wsa", img_date, "_B", band, ".tif"), format = "GTiff", overwrite = T)
    
    
  }
  
}  
albedo_function_c16 <- function(S2_cluster, S2_surfref) {
  
  for (i in 1:6){
    
    band = i
    
    ratio_table <- read.csv(paste0("Ratios_", img_date, "_bsa_B", band, ".csv"), sep = ";", header = T)
    
    ratio_table <- ratio_table[!is.na(ratio_table$Cluster), ]
    
    cluster_final <- as.matrix(S2_cluster)
    albedo_dir <- as.matrix(S2_surfref[[i]])
    
    Albedo_Hemispherical <- ifelse(cluster_final == 1, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 1]),
                                   ifelse(cluster_final == 2, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 2]),
                                          ifelse(cluster_final == 3, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 3]),
                                                 ifelse(cluster_final == 4, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 4]),
                                                        ifelse(cluster_final == 5, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 5]),
                                                               ifelse(cluster_final == 6, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 6]),
                                                                      ifelse(cluster_final == 7, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 7]),
                                                                             ifelse(cluster_final == 8, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 8]),
                                                                                    ifelse(cluster_final == 9, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 9]),
                                                                                           ifelse(cluster_final == 10, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 10]),
                                                                                                  ifelse(cluster_final == 11, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 11]),
                                                                                                         ifelse(cluster_final == 12, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 12]),
                                                                                                                ifelse(cluster_final == 13, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 13]),
                                                                                                                       ifelse(cluster_final == 14, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 14]),
                                                                                                                              ifelse(cluster_final == 15, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 15]),
                                                                                                                                     ifelse(cluster_final == 16, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 16]),
                                                                                                                                            NA))))))))))))))))
    
    
    albedo_bsa <- S2_surfref[[i]]
    
    albedo_bsa[] <- Albedo_Hemispherical
    
    writeRaster(albedo_bsa, filename = paste0("albedo_bsa", img_date, "_B", band, ".tif"), format = "GTiff", overwrite = T)
    
    
  }
  
  
  plot(albedo_bsa)
  
  
  
  for (i in 1:6){
    
    band = i
    
    ratio_table <- read.csv(paste0("Ratios_", img_date, "_wsa_B", band, ".csv"), sep = ";", header = T)
    
    ratio_table <- ratio_table[!is.na(ratio_table$Cluster), ]
    
    cluster_final <- as.matrix(S2_cluster)
    albedo_dir <- as.matrix(S2_surfref[[i]])
    
    Albedo_Hemispherical <- ifelse(cluster_final == 1, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 1]),
                                   ifelse(cluster_final == 2, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 2]),
                                          ifelse(cluster_final == 3, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 3]),
                                                 ifelse(cluster_final == 4, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 4]),
                                                        ifelse(cluster_final == 5, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 5]),
                                                               ifelse(cluster_final == 6, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 6]),
                                                                      ifelse(cluster_final == 7, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 7]),
                                                                             ifelse(cluster_final == 8, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 8]),
                                                                                    ifelse(cluster_final == 9, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 9]),
                                                                                           ifelse(cluster_final == 10, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 10]),
                                                                                                  ifelse(cluster_final == 11, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 11]),
                                                                                                         ifelse(cluster_final == 12, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 12]),
                                                                                                                ifelse(cluster_final == 13, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 13]),
                                                                                                                       ifelse(cluster_final == 14, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 14]),
                                                                                                                              ifelse(cluster_final == 15, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 15]),
                                                                                                                                     ifelse(cluster_final == 16, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 16]),
                                                                                                                                            NA))))))))))))))))
    
    
    albedo_bsa <- S2_surfref[[i]]
    
    albedo_bsa[] <- Albedo_Hemispherical
    
    writeRaster(albedo_bsa, filename = paste0("albedo_wsa", img_date, "_B", band, ".tif"), format = "GTiff", overwrite = T)
    
    
  }
  
  
  plot(albedo_bsa)  
  
  
}
albedo_function_c17 <- function(S2_cluster, S2_surfref) {
  
  for (i in 1:6){
    
    band = i
    
    ratio_table <- read.csv(paste0("Ratios_", img_date, "_bsa_B", band, ".csv"), sep = ";", header = T)
    
    ratio_table  <- ratio_table[!is.na(ratio_table$Cluster), ]
    
    cluster_final <- as.matrix(S2_cluster)
    albedo_dir <- as.matrix(S2_surfref[[i]])
    
    Albedo_Hemispherical <- ifelse(cluster_final == 1, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 1]),
                                   ifelse(cluster_final == 2, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 2]),
                                          ifelse(cluster_final == 3, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 3]),
                                                 ifelse(cluster_final == 4, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 4]),
                                                        ifelse(cluster_final == 5, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 5]),
                                                               ifelse(cluster_final == 6, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 6]),
                                                                      ifelse(cluster_final == 7, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 7]),
                                                                             ifelse(cluster_final == 8, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 8]),
                                                                                    ifelse(cluster_final == 9, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 9]),
                                                                                           ifelse(cluster_final == 10, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 10]),
                                                                                                  ifelse(cluster_final == 11, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 11]),
                                                                                                         ifelse(cluster_final == 12, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 12]),
                                                                                                                ifelse(cluster_final == 13, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 13]),
                                                                                                                       ifelse(cluster_final == 14, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 14]),
                                                                                                                              ifelse(cluster_final == 15, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 15]),
                                                                                                                                     ifelse(cluster_final == 16, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 16]),
                                                                                                                                            ifelse(cluster_final == 17, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 17]),
                                                                                                                                                   NA)))))))))))))))))
    
    
    albedo_bsa <- S2_surfref[[i]]
    
    albedo_bsa[] <- Albedo_Hemispherical
    
    writeRaster(albedo_bsa, filename = paste0("albedo_bsa", img_date, "_B", band, ".tif"), format = "GTiff", overwrite = T)
    
    
  }
  
  
  plot(albedo_bsa)
  
  
  
  for (i in 1:6){
    
    band = i
    
    ratio_table <- read.csv(paste0("Ratios_", img_date, "_wsa_B", band, ".csv"), sep = ";", header = T)
    
    ratio_table  <- ratio_table[!is.na(ratio_table$Cluster), ]
    
    cluster_final <- as.matrix(S2_cluster)
    albedo_dir <- as.matrix(S2_surfref[[i]])
    
    Albedo_Hemispherical <- ifelse(cluster_final == 1, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 1]),
                                   ifelse(cluster_final == 2, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 2]),
                                          ifelse(cluster_final == 3, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 3]),
                                                 ifelse(cluster_final == 4, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 4]),
                                                        ifelse(cluster_final == 5, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 5]),
                                                               ifelse(cluster_final == 6, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 6]),
                                                                      ifelse(cluster_final == 7, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 7]),
                                                                             ifelse(cluster_final == 8, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 8]),
                                                                                    ifelse(cluster_final == 9, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 9]),
                                                                                           ifelse(cluster_final == 10, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 10]),
                                                                                                  ifelse(cluster_final == 11, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 11]),
                                                                                                         ifelse(cluster_final == 12, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 12]),
                                                                                                                ifelse(cluster_final == 13, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 13]),
                                                                                                                       ifelse(cluster_final == 14, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 14]),
                                                                                                                              ifelse(cluster_final == 15, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 15]),
                                                                                                                                     ifelse(cluster_final == 16, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 16]),
                                                                                                                                            ifelse(cluster_final == 17, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 17]),
                                                                                                                                                   NA)))))))))))))))))
    
    
    albedo_bsa <- S2_surfref[[i]]
    
    albedo_bsa[] <- Albedo_Hemispherical
    
    writeRaster(albedo_bsa, filename = paste0("albedo_wsa", img_date, "_B", band, ".tif"), format = "GTiff", overwrite = T)
    
    
  }
  
  
  plot(albedo_bsa)  
  
  
}
albedo_function_c18 <- function(S2_cluster, S2_surfref) {
  
  for (i in 1:6){
    
    band = i
    
    ratio_table <- read.csv(paste0("Ratios_", img_date, "_bsa_B", band, ".csv"), sep = ";", header = T)
    
    ratio_table  <- ratio_table[!is.na(ratio_table$Cluster), ]
    
    cluster_final <- as.matrix(S2_cluster)
    albedo_dir <- as.matrix(S2_surfref[[i]])
    
    Albedo_Hemispherical <- ifelse(cluster_final == 1, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 1]),
                                   ifelse(cluster_final == 2, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 2]),
                                          ifelse(cluster_final == 3, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 3]),
                                                 ifelse(cluster_final == 4, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 4]),
                                                        ifelse(cluster_final == 5, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 5]),
                                                               ifelse(cluster_final == 6, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 6]),
                                                                      ifelse(cluster_final == 7, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 7]),
                                                                             ifelse(cluster_final == 8, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 8]),
                                                                                    ifelse(cluster_final == 9, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 9]),
                                                                                           ifelse(cluster_final == 10, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 10]),
                                                                                                  ifelse(cluster_final == 11, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 11]),
                                                                                                         ifelse(cluster_final == 12, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 12]),
                                                                                                                ifelse(cluster_final == 13, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 13]),
                                                                                                                       ifelse(cluster_final == 14, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 14]),
                                                                                                                              ifelse(cluster_final == 15, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 15]),
                                                                                                                                     ifelse(cluster_final == 16, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 16]),
                                                                                                                                            ifelse(cluster_final == 17, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 17]),
                                                                                                                                                   ifelse(cluster_final == 18, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 18]),
                                                                                                                                                          NA))))))))))))))))))
    
    
    albedo_bsa <- S2_surfref[[i]]
    
    albedo_bsa[] <- Albedo_Hemispherical
    
    writeRaster(albedo_bsa, filename = paste0("albedo_bsa", img_date, "_B", band, ".tif"), format = "GTiff", overwrite = T)
    
    
  }
  
  
  plot(albedo_bsa)
  
  
  
  for (i in 1:6){
    
    band = i
    
    ratio_table <- read.csv(paste0("Ratios_", img_date, "_wsa_B", band, ".csv"), sep = ";", header = T)
    
    ratio_table  <- ratio_table[!is.na(ratio_table$Cluster), ]
    
    cluster_final <- as.matrix(S2_cluster)
    albedo_dir <- as.matrix(S2_surfref[[i]])
    
    Albedo_Hemispherical <- ifelse(cluster_final == 1, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 1]),
                                   ifelse(cluster_final == 2, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 2]),
                                          ifelse(cluster_final == 3, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 3]),
                                                 ifelse(cluster_final == 4, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 4]),
                                                        ifelse(cluster_final == 5, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 5]),
                                                               ifelse(cluster_final == 6, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 6]),
                                                                      ifelse(cluster_final == 7, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 7]),
                                                                             ifelse(cluster_final == 8, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 8]),
                                                                                    ifelse(cluster_final == 9, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 9]),
                                                                                           ifelse(cluster_final == 10, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 10]),
                                                                                                  ifelse(cluster_final == 11, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 11]),
                                                                                                         ifelse(cluster_final == 12, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 12]),
                                                                                                                ifelse(cluster_final == 13, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 13]),
                                                                                                                       ifelse(cluster_final == 14, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 14]),
                                                                                                                              ifelse(cluster_final == 15, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 15]),
                                                                                                                                     ifelse(cluster_final == 16, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 16]),
                                                                                                                                            ifelse(cluster_final == 17, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 17]),
                                                                                                                                                   ifelse(cluster_final == 18, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 18]),
                                                                                                                                                          NA))))))))))))))))))
    
    
    albedo_bsa <- S2_surfref[[i]]
    
    albedo_bsa[] <- Albedo_Hemispherical
    
    writeRaster(albedo_bsa, filename = paste0("albedo_wsa", img_date, "_B", band, ".tif"), format = "GTiff", overwrite = T)
    
    
  }
  
  
  plot(albedo_bsa)  
  
  
}
albedo_function_c19 <- function(S2_cluster, S2_surfref) {
  
  for (i in 1:6){
    
    band = i
    
    ratio_table <- read.csv(paste0("Ratios_", img_date, "_bsa_B", band, ".csv"), sep = ";", header = T)
    
    ratio_table  <- ratio_table[!is.na(ratio_table$Cluster), ]
    
    cluster_final <- as.matrix(S2_cluster)
    albedo_dir <- as.matrix(S2_surfref[[i]])
    
    Albedo_Hemispherical <- ifelse(cluster_final == 1, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 1]),
                                   ifelse(cluster_final == 2, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 2]),
                                          ifelse(cluster_final == 3, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 3]),
                                                 ifelse(cluster_final == 4, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 4]),
                                                        ifelse(cluster_final == 5, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 5]),
                                                               ifelse(cluster_final == 6, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 6]),
                                                                      ifelse(cluster_final == 7, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 7]),
                                                                             ifelse(cluster_final == 8, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 8]),
                                                                                    ifelse(cluster_final == 9, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 9]),
                                                                                           ifelse(cluster_final == 10, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 10]),
                                                                                                  ifelse(cluster_final == 11, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 11]),
                                                                                                         ifelse(cluster_final == 12, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 12]),
                                                                                                                ifelse(cluster_final == 13, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 13]),
                                                                                                                       ifelse(cluster_final == 14, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 14]),
                                                                                                                              ifelse(cluster_final == 15, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 15]),
                                                                                                                                     ifelse(cluster_final == 16, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 16]),
                                                                                                                                            ifelse(cluster_final == 17, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 17]),
                                                                                                                                                   ifelse(cluster_final == 18, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 18]),
                                                                                                                                                          ifelse(cluster_final == 19, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 19]),
                                                                                                                                                                 NA)))))))))))))))))))
    
    
    albedo_bsa <- S2_surfref[[i]]
    
    albedo_bsa[] <- Albedo_Hemispherical
    
    writeRaster(albedo_bsa, filename = paste0("albedo_bsa", img_date, "_B", band, ".tif"), format = "GTiff", overwrite = T)
    
    
  }
  
  
  plot(albedo_bsa)
  
  
  
  for (i in 1:6){
    
    band = i
    
    ratio_table <- read.csv(paste0("Ratios_", img_date, "_wsa_B", band, ".csv"), sep = ";", header = T)
    
    ratio_table  <- ratio_table[!is.na(ratio_table$Cluster), ]
    
    cluster_final <- as.matrix(S2_cluster)
    albedo_dir <- as.matrix(S2_surfref[[i]])
    
    Albedo_Hemispherical <- ifelse(cluster_final == 1, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 1]),
                                   ifelse(cluster_final == 2, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 2]),
                                          ifelse(cluster_final == 3, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 3]),
                                                 ifelse(cluster_final == 4, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 4]),
                                                        ifelse(cluster_final == 5, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 5]),
                                                               ifelse(cluster_final == 6, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 6]),
                                                                      ifelse(cluster_final == 7, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 7]),
                                                                             ifelse(cluster_final == 8, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 8]),
                                                                                    ifelse(cluster_final == 9, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 9]),
                                                                                           ifelse(cluster_final == 10, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 10]),
                                                                                                  ifelse(cluster_final == 11, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 11]),
                                                                                                         ifelse(cluster_final == 12, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 12]),
                                                                                                                ifelse(cluster_final == 13, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 13]),
                                                                                                                       ifelse(cluster_final == 14, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 14]),
                                                                                                                              ifelse(cluster_final == 15, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 15]),
                                                                                                                                     ifelse(cluster_final == 16, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 16]),
                                                                                                                                            ifelse(cluster_final == 17, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 17]),
                                                                                                                                                   ifelse(cluster_final == 18, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 18]),
                                                                                                                                                          ifelse(cluster_final == 19, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 19]),
                                                                                                                                                                 NA)))))))))))))))))))
    
    
    albedo_bsa <- S2_surfref[[i]]
    
    albedo_bsa[] <- Albedo_Hemispherical
    
    writeRaster(albedo_bsa, filename = paste0("albedo_wsa", img_date, "_B", band, ".tif"), format = "GTiff", overwrite = T)
    
    
  }
  
  
  plot(albedo_bsa)  
  
  
}
albedo_function_c20 <- function(S2_cluster, S2_surfref) {
  
  for (i in 1:6){
    
    band = i
    
    ratio_table <- read.csv(paste0("Ratios_", img_date, "_bsa_B", band, ".csv"), sep = ";", header = T)
    
    ratio_table  <- ratio_table[!is.na(ratio_table$Cluster), ]
    
    cluster_final <- as.matrix(S2_cluster)
    albedo_dir <- as.matrix(S2_surfref[[i]])
    
    Albedo_Hemispherical <- ifelse(cluster_final == 1, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 1]),
                                   ifelse(cluster_final == 2, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 2]),
                                          ifelse(cluster_final == 3, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 3]),
                                                 ifelse(cluster_final == 4, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 4]),
                                                        ifelse(cluster_final == 5, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 5]),
                                                               ifelse(cluster_final == 6, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 6]),
                                                                      ifelse(cluster_final == 7, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 7]),
                                                                             ifelse(cluster_final == 8, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 8]),
                                                                                    ifelse(cluster_final == 9, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 9]),
                                                                                           ifelse(cluster_final == 10, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 10]),
                                                                                                  ifelse(cluster_final == 11, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 11]),
                                                                                                         ifelse(cluster_final == 12, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 12]),
                                                                                                                ifelse(cluster_final == 13, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 13]),
                                                                                                                       ifelse(cluster_final == 14, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 14]),
                                                                                                                              ifelse(cluster_final == 15, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 15]),
                                                                                                                                     ifelse(cluster_final == 16, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 16]),
                                                                                                                                            ifelse(cluster_final == 17, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 17]),
                                                                                                                                                   ifelse(cluster_final == 18, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 18]),
                                                                                                                                                          ifelse(cluster_final == 19, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 19]),
                                                                                                                                                                 ifelse(cluster_final == 20, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 20]),
                                                                                                                                                                        NA))))))))))))))))))))
    
    
    albedo_bsa <- S2_surfref[[i]]
    
    albedo_bsa[] <- Albedo_Hemispherical
    
    writeRaster(albedo_bsa, filename = paste0("albedo_bsa", img_date, "_B", band, ".tif"), format = "GTiff", overwrite = T)
    
    
  }
  
  
  plot(albedo_bsa)
  
  
  
  for (i in 1:6){
    
    band = i
    
    ratio_table <- read.csv(paste0("Ratios_", img_date, "_wsa_B", band, ".csv"), sep = ";", header = T)
    
    ratio_table  <- ratio_table[!is.na(ratio_table$Cluster), ]
    
    cluster_final <- as.matrix(S2_cluster)
    albedo_dir <- as.matrix(S2_surfref[[i]])
    
    Albedo_Hemispherical <- ifelse(cluster_final == 1, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 1]),
                                   ifelse(cluster_final == 2, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 2]),
                                          ifelse(cluster_final == 3, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 3]),
                                                 ifelse(cluster_final == 4, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 4]),
                                                        ifelse(cluster_final == 5, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 5]),
                                                               ifelse(cluster_final == 6, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 6]),
                                                                      ifelse(cluster_final == 7, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 7]),
                                                                             ifelse(cluster_final == 8, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 8]),
                                                                                    ifelse(cluster_final == 9, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 9]),
                                                                                           ifelse(cluster_final == 10, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 10]),
                                                                                                  ifelse(cluster_final == 11, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 11]),
                                                                                                         ifelse(cluster_final == 12, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 12]),
                                                                                                                ifelse(cluster_final == 13, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 13]),
                                                                                                                       ifelse(cluster_final == 14, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 14]),
                                                                                                                              ifelse(cluster_final == 15, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 15]),
                                                                                                                                     ifelse(cluster_final == 16, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 16]),
                                                                                                                                            ifelse(cluster_final == 17, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 17]),
                                                                                                                                                   ifelse(cluster_final == 18, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 18]),
                                                                                                                                                          ifelse(cluster_final == 19, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 19]),
                                                                                                                                                                 ifelse(cluster_final == 20, albedo_dir*(ratio_table$Ratio[ratio_table$Cluster == 20]),
                                                                                                                                                                        NA))))))))))))))))))))
    
    
    albedo_bsa <- S2_surfref[[i]]
    
    albedo_bsa[] <- Albedo_Hemispherical
    
    writeRaster(albedo_bsa, filename = paste0("albedo_wsa", img_date, "_B", band, ".tif"), format = "GTiff", overwrite = T)
    
    
  }
  
  
  plot(albedo_bsa)  
  
  
}



eval(parse(text=paste("albedo_function_c", j, "(S2_cluster, S2_surfref)", sep = "")))


#Create masks based on snow and ice presence

rm(list=ls()[! ls() %in% c("img_date", "S2_surfref", "AOI2", "Dm_j", "Dm", "M", "Y", "j", "date_j")])

files1 <- list.files(getwd(), pattern = paste0("albedo_bsa", img_date, ".*.tif"))

albedo_bsa <- stack(files1)

files2 <- list.files(getwd(), pattern = paste0("albedo_wsa", img_date, ".*.tif"))

albedo_wsa <- stack(files2)


#Input solar irradiance data

E_in_data <- read.csv("Measured_Solar_Irradiance.csv", sep = ",", header = T)

E_in_meas <- E_in_data$E_in_meas[E_in_data$julian_day == as.numeric(date_j)]


SolarIrradiance <- function(latitude, J) {
  
  #Input variables and parameters 
  lat <- latitude*0.0174533 #rad 
  long <- -115.1389 #deg
  SM <- -105 #deg longitude for local time zone
  Gsc <- 1367 #W.m^-2 Solar Constant
  
  
  #Declination angle
  
  delta = (0.409*sin(((2*pi/365)*J) - 1.39))
  
  delta_deg = delta*57.2958
  
  #Equation of time
  
  t = ((2*pi*J)/366) + 4.8718
  
  E =  (5.0323 -430.847*cos(t) +12.5024*cos(2*t) +18.25*cos(3*t) -100.976*sin(t)
        +595.275*sin(2*t) +3.6858*sin(3*t) -12.47*sin(4*t))/60
  
  #Solar noon and solar time
  
  SN = 12 - (E/60) - ((SM - long)/15)
  
  LST = 0:23 #local standard time
  
  #Hour angle
  
  h = (LST-SN)*(pi/12)
  
  #Zenith angle
  
  cos_Z = (sin(lat)*sin(delta)) + (cos(lat)*cos(delta)*cos(h))
  
  #Hourly extraterrestrial radiation
  
  dr = 1 + 0.033*cos(((2*pi)/365)*J)
  
  Ra_hor = Gsc*dr*cos_Z
  
  Ra_hor = ifelse(Ra_hor < 0, 0, Ra_hor)
  
  return(Ra_hor[13])
  
}

#Calculate directional blue-sky albedo conditionally based on the presence of snow (Li et al., 2018)

E_in <- SolarIrradiance(52.156971, Dm_j)


kt <- E_in_meas/E_in


albedo_bsa_int <- (-0.0001 + (albedo_bsa[[1]]*-0.1992) + (albedo_bsa[[2]]*2.3002)
                   + (albedo_bsa[[3]]*-1.9121) + (albedo_bsa[[4]]*0.6715)
                   + (albedo_bsa[[5]]*-2.2728) + (albedo_bsa[[6]]*1.9341))

albedo_wsa_int <- (-0.0001 + (albedo_wsa[[1]]*-0.1992) + (albedo_wsa[[2]]*2.3002)
                   + (albedo_wsa[[3]]*-1.9121) + (albedo_wsa[[4]]*0.6715)
                   + (albedo_wsa[[5]]*-2.2728) + (albedo_wsa[[6]]*1.9341))


f_dif <- 1.1 - (1.09*kt) #from Ellis and Pomeroy (2007)

blue_sky_albedo = ((1 - f_dif)*albedo_bsa_int) + (f_dif*albedo_wsa_int)

plot(blue_sky_albedo)


#Create masks based on snow and ice presence

S2_surfref <- crop(S2_surfref, AOI2)

snow_mask <- S2_surfref[[1]]

snow_mask[] <- NA
snow_mask[(S2_surfref[[2]] - S2_surfref[[5]])/(S2_surfref[[2]] + S2_surfref[[5]]) > 0.4] <- 1
snow_mask[S2_surfref[[4]] < 0.1] <- NA


S2_albedo_snow <- blue_sky_albedo*snow_mask


plot(S2_albedo_snow, main = paste0("Sentinel-2 Blue-sky Albedo - ", img_date), col = rev(topo.colors(12)), zlim = c(0,1.2))


writeRaster(S2_albedo_snow, filename = paste0("blue_sky_albedo_snow_", img_date, "_c", j), format = "GTiff", overwrite = T)

rm(list=ls()[! ls() %in% c("img_date", "S2_surfref", "AOI2", "Dm_j", "Dm", "M", "Y", "j", "date_j")])

}

Sys.time()
