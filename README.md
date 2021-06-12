# snow_ice_rs_albedo
Scripts to generate snow and ice albedo from Sentinel-2 surface reflectance and MODIS BRDF information. The inputs are Sentinel-2 surface reflectance images; MODIS MOD09GA surface
reflectance product; and observed incoming/outgoing shortwave radiation, air temperature, and relative humidity. The meteorological information can also be estimated by the user
from another sources if a station is not available. The scripts need to be run in this order: Snow_Kernel, AN_Ratio, and S2_Albedo_Updt. You can use the Albedo_NSW_Radiation
script if you wish to calculate the net shortwave radiation of snow and ice surfaces.
