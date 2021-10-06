#Calculates net shortwave radiation corrected for slope and aspect using the remote sensing
#albedo retrievals from the S2_Albedo_Updt script. This code is based on Allen et al. (2006) 
#"Analytical integrated functions for daily solar radiation on slopes" - Agricultural and Forest Meteorology

#Author: Andre Bertoncini

#Institution: Centre for Hydrology - University of Saskatchewan


library(raster)

setwd("/path to S2_Albedo_Updt outputs")

observation_vars <- read.table("/radiation observations.csv", header = T, sep = ",", dec = ".")


files <- list.files(getwd(), pattern = "blue_sky_albedo_snow_.*.tif$")

s2_albedo_stack <- stack(files)


#The code below will loop through all your albedo images

for (i in 1:19) {
  
  s2_albedo_stack[[i]][s2_albedo_stack[[i]] > 1] <- NA
  s2_albedo_stack[[i]][s2_albedo_stack[[i]] < 0] <- NA
  
  
  #Parameters
  
  J <- as.numeric(substr(observation_vars$julian_day[i],5,7))
  Ta <- observation_vars$Ta_meas[i] #observed air temperature in deg. C
  RH <- observation_vars$RH_meas[i] #observed relative humidity in %
  T0 <- 273.15 #K
  Rv <- 461 #J.K^-1.kg^-1
  e0 <- 0.6113 #kPa
  Kt <- c(1,1,1,1,1,1,1,1,0.75,0.75,1,0.75,1,1,1,1,1,1,1) #empirical turbidity coefficient from Allen et al. (2006)
  SM <- -105 #deg longitude for local time zone
  
  
  srtm <- raster("/path to SRTM elevation in UTM.tif")
  srtm_geog <- raster("/path to SRTM in geographic coordinates.tif")
  
  plot(srtm)
  
  lat <- srtm_geog
  long <- srtm_geog
  xy <- coordinates(srtm_geog)
  lat[] <- xy[,2]
  long[] <- xy[,1]
  
  lat <- mean(as.vector(lat))
  long <- mean(as.vector(long))
  alt <- mean(as.vector(srtm), na.rm = T)
  
  
  #Calculates spatially-distributed Cos Z based on SRTM elevation
  
  #Hour angle
  
  t = ((2*pi*J)/366) + 4.8718
  
  E =  (5.0323 -430.847*cos(t) +12.5024*cos(2*t) +18.25*cos(3*t) -100.976*sin(t)
        +595.275*sin(2*t) +3.6858*sin(3*t) -12.47*sin(4*t))/60
  
  #Solar noon and solar time
  
  SN = 12 - (E/60) - ((SM - long)/15)
  
  LST = 13 #local standard time
  
  #Hour angle
  
  h = (LST-SN)*(pi/12)
  
  #Calculates slope
  
  s <- terrain(srtm, opt = "slope", unit = "radians")
  
  plot(s)
  
  #Calculates aspect
  
  gamma <- terrain(srtm, opt = "aspect", unit = "radians")
  
  gamma <- gamma - pi
  
  plot(gamma)
  
  
  #Calculates Cos Z corrected for slope and aspect
  
  delta = (0.409*sin(((2*pi/365)*J) - 1.39))
  
  cos_Z_slope = ((sin(delta)*sin(lat*(pi/180))*cos(s))
                 -(sin(delta)*cos(lat*(pi/180))*sin(s)*cos(gamma))
                 +(cos(delta)*cos(lat*(pi/180))*cos(s)*cos(h))
                 +(cos(delta)*sin(lat*(pi/180))*sin(s)*cos(gamma)*cos(h))
                 +(cos(delta)*sin(gamma)*sin(s)*sin(h)))
  
  
  plot(cos_Z_slope)
  
  
  #Calculates net shortwave radiation
  
  E_in_meas <- observation_vars$E_in_meas[i] #incoming shortwave radiation from station (image, reference image)
  
  alb <- s2_albedo_stack[[i]] #albedo images
  
  
  #Function to calculate point-based theoretical solar irradiance with the purpose of calculating station transmissivity
  
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
    
    LST = 13 #local standard time
    
    #Hour angle
    
    h = (LST-SN)*(pi/12)
    
    #Zenith angle
    
    cos_Z = (sin(lat)*sin(delta)) + (cos(lat)*cos(delta)*cos(h))
    
    #Hourly extraterrestrial radiation
    
    dr = 1 + 0.033*cos(((2*pi)/365)*J)
    
    Ra_hor = Gsc*dr*cos_Z
    
    Ra_hor = ifelse(Ra_hor < 0, 0, Ra_hor)
    
    return(Ra_hor)
    
  }
  
  E_in <- SolarIrradiance(lat, J)
  
  
  transm <- E_in_meas/E_in
  
  transm
  
  kb = 1.56*transm - 0.55
  
  kd = transm - kb
  
  
  #Function to calculate spatially-distributed net shortwave radiation corrected for slope and aspect
  
  SolarIrrSlope <- function(latitude, J, cos_Z_slope) {
    
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
    
    LST = 13 #local standard time
    
    #Hour angle
    
    h = (LST-SN)*(pi/12)
    
    #Hourly extraterrestrial radiation
    
    dr = 1 + 0.033*cos(((2*pi)/365)*J)
    
    Ra_hor = Gsc*dr*cos_Z_slope
    
    return(Ra_hor)
    
  }
  
  E_in_slope <- SolarIrrSlope(lat, J, cos_Z_slope)
  
  E_in_slope[E_in_slope < 0] <- 0
  
  L <- ifelse(Ta < 0, 2.83e+6, 2.5e+6)
  
  partvap <- function(e0, L, Rv, T0, Ta){
    es <- e0*exp((L/Rv)*((1/T0)-(1/(Ta))))
    es
  }
  
  es <- partvap(e0, L, Rv, T0, Ta+273.15)
  
  ea = (RH/100)*es
  
  P_alt = 101.3*exp(-alt/8200)
  
  W = 0.14*ea*P_alt + 2.1
  
  
  #Calculates clearness index
  
  elev_ang =  ((sin(lat*(pi/180)))*(sin(delta))) + ((cos(lat*(pi/180)))*(cos(delta))*(cos(h)))
  
  kbo = 0.98*exp(((-0.00146*P_alt)/(Kt[i]*sin(elev_ang)))-(0.075*((W/sin(elev_ang))^0.4)))
  
  fb = (kbo/kb)*(E_in_slope/E_in_meas)
  
  fi = 0.75 + (0.25*(cos(s))) - ((0.5*s)/pi)
  
  fia = (1 - kb)*((1+((kb/(kb+kd))^0.5)*(sin(s/2)^3))*fi) + fb*kb
  
  albedo_slope <- resample(alb, cos_Z_slope)
  
  Rs = E_in_meas*((fb*(kb/transm)) + (fia*(kd/transm)) + (albedo_slope*(1 - fi)))
  
  plot(Rs)
  
  
  Rs_hor = Rs/(cos(s))
  
  plot(Rs_hor)
  
  
  NSW_radiation = Rs_hor*(1 - albedo_slope)
  
  plot(NSW_radiation)
  
  writeRaster(NSW_radiation, filename = paste0("/path to SW_Forcing outputs/SW_Forcing_", i), 
              format = "GTiff", overwrite = T)
  
}

