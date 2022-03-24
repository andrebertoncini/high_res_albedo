# High Spatial Resolution Albedo
Scripts to generate high spatial resolution albedo from Sentinel-2 surface reflectance and MODIS BRDF information. The inputs are Sentinel-2 surface reflectance images; MODIS MOD09GA surface
reflectance product; and observed incoming/outgoing shortwave radiation, air temperature, and relative humidity. The meteorological information can also be estimated by the user
from another sources if a station is not available. The scripts need to be run in this order: RTLSR_Kernel or Snow_Kernel, AN_Ratio, and S2_Albedo_Updt. You can choose to use the RTLSR_Kernel or the Snow_Kernel to model the BRDF, but remember to check performance for your area to choose which model to use. You can use the Albedo_NSW_Radiation
script if you wish to calculate the net shortwave radiation of snow and ice surfaces.
