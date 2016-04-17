
#### Read MODIS data for Aerosol Optical Depth ##############
#### MODIS Collection 6 (3Km resolution AOD) ################

## download data from here:
# aftp://karafede:Password07@nrt3.modaps.eosdis.nasa.gov/allData/6/MOD04_L2/ftp://ladsftp.nascom.nasa.gov/allData/6/MOD04_L2/

setwd("C:/RICARDO-AEA/SATELLITE_STUFF/AOD_MOD04_L2_3K/MOD04_AOD_2016_103/")
# setwd("C:/SATELLITE_STUFF/AOD_MOD04_L2_3K/MOD04_AOD_2016_103/")

# source("http://bioconductor.org/biocLite.R")
# biocLite("rhdf5")

# library(rhdf5)
# library(maps)
# library(RNetCDF)
library(RCurl)
library(rgdal)
library(gdalUtils)
library(raster)
library(threadr)
library(readr)
library(MODISTools)
library(mapview)
library(ncdf4)


filenames <- list.files(pattern = "\\.hdf$")
# dir_WW <- "C:/RICARDO-AEA/SATELLITE_STUFF/AOD_MOD04_L2_3K/Admin_shp"
dir_WW <- "C:/SATELLITE_STUFF/AOD_MOD04_L2_3K/Admin_shp"
shp <- readOGR(dsn = dir_WW, layer = "admin01")

# ----- Transform to EPSG 4326 - WGS84 (required)
shp <- spTransform(shp, CRS("+init=epsg:4326"))
crs <- projection(shp) ### get projections from shp file
plot(shp)

#######################################################################

# Get a list of sds names

sds <- get_subdatasets(filenames[1])  ## it is going to open QGIS internally and using GDAL libraries
# Isolate the name of the first sds
name <- sds[64]  # AOD_550_Dark_Target_Deep_Blue_Combined
filename <- rasterTmpFile()
extension(filename) <- 'tif'
gdal_translate(sds[64], dst_dataset = filename)
# Load the Geotiff created into R
r <- raster(filename)
# crs <- projection(r) ### get projections from r file (which projections are these?....no projection!!!!)
# projection(r) <- CRS("+proj=longlat +datum=WGS84")
plot(r)

# plot(shp, add=TRUE, lwd=2)

# prova_raster <- writeRaster(r,
#                               filename="prova_raster.nc",
#                               format="CDF", overwrite=TRUE)

lon <- sds[71]
filename <- rasterTmpFile()
extension(filename) <- 'tif'
gdal_translate(sds[71], dst_dataset = filename)
# Load the Geotiff created into R
lon <- raster(filename)
plot(lon)


lat <- sds[72]
filename <- rasterTmpFile()
extension(filename) <- 'tif'
gdal_translate(sds[72], dst_dataset = filename)
# Load the Geotiff created into R
lat <- raster(filename)
plot(lat)


### Exctract poitns from raster ######################################

r <- rasterToPoints(r)
head(r)
colnames(r) <- c("Lon", "Lat", "values")
r <- as.data.frame (r)
r <- subset(r, !is.na(values) & values>0)
write.csv(r, file = "r.csv", row.names=FALSE)


# make a function to generate all .csv files with Lon, lat and value

extract_HDF <- function (file) {      ## this is the filenames 
  
  # get list of field names
  nome <- str_sub(file, start = 1, end = -9)
  sds <- get_subdatasets(file)
  name <- sds[64]
  filename <- rasterTmpFile()
  extension(filename) <- 'tif'
  # get raster data
  gdal_translate(sds[64], dst_dataset = filename)
  r <- raster(filename)
  # extract values from raster
  r <- rasterToPoints(r)
  colnames(r) <- c("Lon", "Lat", "values")
  r <- as.data.frame (r)
  r <- subset(r, !is.na(values) & values>0)
  write.csv(r, file = paste(nome,".csv", sep = ""), row.names=FALSE)
  
}


# AAA <- extract_HDF(filenames[7])
# Apply function to all files
BBB <- lapply(filenames[1:5], extract_HDF)
BBB <- lapply(filenames, extract_HDF)


######################################################################################
######################################################################################

# change working directory where .csv file of world wide tiles are stored (10km resolution)
setwd("C:/RICARDO-AEA/SATELLITE_STUFF/AOD_MOD04_L2_3K/MOD04_AOD_2016_103")

filenames_tiles <- list.files(pattern = "\\.csv$")

LAT = NULL
LON = NULL
aod = NULL

# filenames_tiles <- filenames_tiles[1:5]

## Bind all data together 
for (i in 1:length(filenames_tiles)) {
  lon <- read_csv(filenames_tiles[i])[,1]
  lat <- read_csv(filenames_tiles[i])[,2]
  AOD <- read_csv(filenames_tiles[i])[,3]
  LON = rbind(LON, data.frame(lon))
  LAT = rbind(LAT, data.frame(lat))
  aod = rbind(aod, data.frame(AOD))
}

DATA <- cbind(LON, LAT, aod)
write.csv(DATA, "AOD_MOD04_10km_WW.csv")
head(DATA)

#crop to desired area
# Saudi Arabia

# DATA_SAUDI <- subset(DATA, lon <= 180 & lon >= -180 &
#                 lat >= 20 & lat <= 50)
# head(DATA_SAUDI)


#### create a raster for AOD at World Wide Level #######

coordinates(DATA) <- ~ lon + lat
# coerce to SpatialPixelsDataFrame
gridded(DATA) <- TRUE
raster_AOD <- raster(DATA)
projection(raster_AOD) <- CRS("+proj=longlat +datum=WGS84")
plot(raster_AOD)
plot(shp)

#### crop the raster over the shp file ###########################
raster_PM25_UK_AIR_cropped <- crop(raster_PM25_UK_AIR, extent(shp))
raster_PM25_UK_AIR_cropped <- mask(raster_PM25_UK_AIR_cropped, shp)

plot(raster_PM25_UK_AIR_cropped)
# plot(shp, add=TRUE, lwd=2)

PM25_UK_AIR_nc <- writeRaster(raster_PM25_UK_AIR_cropped,
                              filename="PM25_UK_AIR.nc",
                              format="CDF", overwrite=TRUE) 
PM25_UK_AIR_nc <- raster("PM25_UK_AIR.nc")





