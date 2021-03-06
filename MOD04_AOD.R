
#### Read MODIS data for Aerosol Optical Depth ##############
#### MODIS Collection 6 (3Km resolution AOD) ################

# setwd("C:/RICARDO-AEA/SATELLITE_STUFF/AOD_MOD04_L2_3K/MOD04_AOD_2016_103/")
# setwd("C:/SATELLITE_STUFF/AOD_MOD04_L2_3K/MOD04_AOD_2016_103/")
setwd("E:/Ricardo-AEA/SATELLITE_STUFF/AOD_MOD04_L2_3K/MOD04_AOD_3K_2016_124/")
# new data from v� october 2016
# setwd("E:/Ricardo-AEA/SATELLITE_STUFF/AOD_MOD04_L2_3K/MOD04_AOD_3K_2016_277_3Oct2016")
setwd("E:/Ricardo-AEA/SATELLITE_STUFF/AOD_MOD04_L2_3K/MOD04_AOD_10K_2016_277_3Oct2016")

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
library(dplyr)
library(gissr)
library(gstat)
library(htmlwidgets)


filenames <- list.files(pattern = "\\.hdf$")
# dir_WW <- "C:/RICARDO-AEA/SATELLITE_STUFF/AOD_MOD04_L2_3K/Admin_shp"
dir_WW <- "C:/SATELLITE_STUFF/AOD_MOD04_L2_3K/Admin_shp"
shp <- readOGR(dsn = dir_WW, layer = "admin01")

# ----- Transform to EPSG 4326 - WGS84 (required)
shp <- spTransform(shp, CRS("+init=epsg:4326"))
crs <- projection(shp) ### get projections from shp file
# plot(shp)


#######################################################################
# to CHECK!!!!! jusr run once to check!
 # Get a list of sds names

sds <- get_subdatasets("MOD04_3K_prova.hdf") # chose a trial file to see the names
sds <- get_subdatasets("MODIS_prova.hdf") # chose a trial file to see the names
  
# Isolate the name of the first sds
# name <- sds[64]  # AOD_550_Dark_Target_Deep_Blue_Combined 10km
name <- sds[12]  # AOT Image_Optical_Depth_Land_And_Ocean 3km
filename <- rasterTmpFile()
extension(filename) <- 'tif'
# gdal_translate(sds[12], dst_dataset = filename)
gdal_translate(sds[64], dst_dataset = filename)
# Load the Geotiff created into R
r <- raster(filename)
plot(r)


# lon <- sds[71]  #longitude  #10km
lon <- sds[51]  #longitude  #3km
filename <- rasterTmpFile()
extension(filename) <- 'tif'
# gdal_translate(sds[51], dst_dataset = filename) # 3km
gdal_translate(sds[71], dst_dataset = filename) # 10km
# Load the Geotiff created into R
lon <- raster(filename)
plot(lon)


# lat <- sds[72]  #longitude 10km
lat <- sds[52] # latitude 3km
filename <- rasterTmpFile()
extension(filename) <- 'tif'
gdal_translate(sds[72], dst_dataset = filename) #10km
# gdal_translate(sds[52], dst_dataset = filename) #3km
# Load the Geotiff created into R
lat <- raster(filename)
plot(lat)


# ### Exctract poitns from raster ######################################
# 
# # data values for AOD_550_Dark_Target_Deep_Blue_Combined
# values <- rasterToPoints(r)  # NA values are not converted
# head(values)
# colnames(values) <- c("x", "y", "values")
# values <- as.data.frame (values)
# # r <- subset(r, !is.na(values) & values>0)
# 
# # data values for longitude
# longitude <- rasterToPoints(lon)  
# head(longitude)
# colnames(longitude) <- c("x", "y", "lon")
# longitude <- as.data.frame (longitude)
# 
# # data values for longitude
# latitude <- rasterToPoints(lat)  
# head(latitude)
# colnames(latitude) <- c("x", "y", "lat")
# latitude <- as.data.frame (latitude)
# 
# 
# # Join  lat, lon 
# Lat_Lon <- latitude %>% 
#   inner_join(longitude, c("x", "y")) 
# 
# Lat_Lon_Values <- Lat_Lon %>% 
#   inner_join(values, c("x", "y")) 
# 
# Lat_Lon_Values <- merge(Lat_Lon, values, all=TRUE)
# 
# Tile <- Lat_Lon_Values %>%
#   select(lon, lat, values)
# # Tile[is.na(Tile)] <- 0
# Tile <- na.omit(Tile)



################################################################################
################################################################################

# make a function to generate all .csv files with Lon, lat and value

extract_HDF <- function (file) {      ## this is the filenames 

  # get list of field names
  nome <- str_sub(file, start = 1, end = -9)
  sds <- get_subdatasets(file)
  
  AOD <- sds[64]  ## field value AOD AOD_550_Dark_Target_Deep_Blue_Combined # 10km
#  AOD <- sds[12]  ## # AOT Image_Optical_Depth_Land_And_Ocean #3km
  filename <- rasterTmpFile()
  extension(filename) <- 'tif'
  # get raster data
#  gdal_translate(sds[12], dst_dataset = filename) #3km
  gdal_translate(sds[64], dst_dataset = filename) #10km
  r <- raster(filename)
  values <- rasterToPoints(r)  
  colnames(values) <- c("x", "y", "values")
  values <- as.data.frame (values)
  
  lon <- sds[71]  #longitude # 10km
#  lon <- sds[51]  #longitude #3km
  filename <- rasterTmpFile()
  extension(filename) <- 'tif'
#  gdal_translate(sds[51], dst_dataset = filename) #3km
  gdal_translate(sds[71], dst_dataset = filename) #10km
  lon <- raster(filename)
  # data values for longitude
  longitude <- rasterToPoints(lon)  
  colnames(longitude) <- c("x", "y", "lon")
  longitude <- as.data.frame (longitude)
  
  
  lat <- sds[72] # latitude  #10km
#  lat <- sds[52] # latitude  #3km
  filename <- rasterTmpFile()
  extension(filename) <- 'tif'
 # gdal_translate(sds[52], dst_dataset = filename) #3km
  gdal_translate(sds[72], dst_dataset = filename) #10km
  # Load the Geotiff created into R
  lat <- raster(filename)
  # data values for longitude
  latitude <- rasterToPoints(lat)  
  colnames(latitude) <- c("x", "y", "lat")
  latitude <- as.data.frame (latitude)
  
  # Join  lat, lon 
  Lat_Lon <- latitude %>% 
    inner_join(longitude, c("x", "y"))
  
   Lat_Lon_Values <- Lat_Lon %>% 
        inner_join(values, c("x", "y"))  
#  Lat_Lon_Values <- merge(Lat_Lon, values, all=TRUE)
    
  MODIS_data <- Lat_Lon_Values %>%
    select(lon, lat, values)
  MODIS_data <- na.omit(MODIS_data)

  write.csv(MODIS_data, file = paste(nome,".csv", sep = ""), row.names=FALSE)
  
}  


# total <- length(filenames)
# pb <- txtProgressBar(min = 0, max = total, style = 3)
# for(i in 1:total){
  
# BBB <- lapply(filenames[1:5], extract_HDF)
BBB <- lapply(filenames, extract_HDF)

# Sys.sleep(0.1)
# # update progress bar
# setTxtProgressBar(pb, i)
# }
# close(pb)

######################################################################################
######################################################################################

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

MODIS04_data <- cbind(LON, LAT, aod)
MODIS04_data <- subset(MODIS04_data, !is.na(values) & !lat == -999 & !lon == -999)

write.csv(MODIS04_data, "AOD_MOD04_10km_WW_03Oct2016.csv")
# write.csv(MODIS04_data, "AOD_MOD04_3K_3km_WW.csv")

head(MODIS04_data)


################################################################################
# subset data for a selected region
# Saudi Arabia & EU

# MODIS04_data <- read.csv("C:/RICARDO-AEA/SATELLITE_STUFF/AOD_MOD04_L2_3K/MOD04_AOD_2016_103/AOD_MOD04_10km_WW.csv")
# MODIS04_data <- read.csv("C:/SATELLITE_STUFF/AOD_MOD04_L2_3K/MOD04_AOD_2016_103/AOD_MOD04_10km_WW.csv")

# MODIS04_data <- read.csv("E:/Ricardo-AEA/SATELLITE_STUFF/AOD_MOD04_L2_3K/MOD04_AOD_3K_2016_124/AOD_MOD04_3K_3km_WW.csv")[2:4]
MODIS04_data <- read.csv("E:/Ricardo-AEA/SATELLITE_STUFF/AOD_MOD04_L2_3K/MOD04_AOD_10K_2016_277_3Oct2016/AOD_MOD04_10km_WW_03Oct2016.csv")[2:4]
DATA_EU_UAE <- subset(MODIS04_data, lon <= 75 & lon >= -22 & lat >= -1.7 & lat <= 66.5)
# DATA_EMIRATES <- subset(MODIS04_data, lon <= 80 & lon >= 10 & lat >= 10 & lat <= 55)

head(DATA_EU_UAE)

# write.csv(DATA_EMIRATES, "AOD_MOD04_10km_EMIRATES.csv")
# write.csv(DATA_EMIRATES, "AOD_MOD04_3K_3km_EMIRATES.csv")
# write.csv(DATA_EU_UAE, "AOD_MOD04_3K_3km_EMIRATES.csv")
write.csv(DATA_EU_UAE, "AOD_MOD04_10km_03Oct2016_EMIRATES.csv")
# DATA_EMIRATES <- read.csv("E:/Ricardo-AEA/SATELLITE_STUFF/AOD_MOD04_L2_3K/MOD04_AOD_3K_2016_124/AOD_MOD04_3K_3km_EMIRATES.csv")
DATA_EMIRATES <- read.csv("E:/Ricardo-AEA/SATELLITE_STUFF/AOD_MOD04_L2_3K/MOD04_AOD_10K_2016_277_3Oct2016/AOD_MOD04_10km_03Oct2016_EMIRATES.csv")[-1]
head(DATA_EMIRATES)


#############################################

# crs <- projection(shp) ### get projections from shp file
# 
# DATA_EMIRATES <- SpatialPointsDataFrame(DATA_EMIRATES[,1:2], DATA_EMIRATES, proj4string=CRS(crs)) 
# pixels <- SpatialPixelsDataFrame(DATA_EMIRATES, tolerance = 0.916421, DATA_EMIRATES@data)
# raster_DATA_EMIRATES <- raster(pixels[,'values'])
# plot(raster_DATA_EMIRATES)
# plot(shp, add=TRUE, lwd=1)


#########################################################
##################################################################################
##################################################################################

# Function for converting a data frame to a SpatialPointsDataFrame
# from the GISSR package from Stuart

data_frame_to_points <- function (df, latitude = "lat", 
                                  longitude = "lon", 
                                  projection = "+proj=longlat +datum=WGS84 +no_defs") {
  
  # Catch for dplyr's data frame class
  df <- threadr::base_df(df)
  
  # Make sp points object
  sp::coordinates(df) <- c(longitude, latitude)
  
  # Reassign
  sp <- df
  
  # Give the object a projection
  sp <- sp_transform(sp, projection, warn = FALSE)
  
  # Return
  sp
  
}


##################################################################################
##################################################################################
# use GISSR package from Stuart
# DATA_EMIRATES <- data_frame_to_points(DATA_EMIRATES, "lat", "lon")

############  Make a data interpolation on a regular grid ########################

min(DATA_EMIRATES$lat)
max(DATA_EMIRATES$lat)
min(DATA_EMIRATES$lon)
max(DATA_EMIRATES$lon)

DATA_EMIRATES$x <- DATA_EMIRATES$lon
DATA_EMIRATES$y <- DATA_EMIRATES$lat

coordinates(DATA_EMIRATES) = ~x + y  ## Set spatial coordinates to create a Spatial object:
plot(DATA_EMIRATES)

x.range <- as.numeric(c(-23, 75))  # min/max longitude of the interpolation area
y.range <- as.numeric(c(-1.7, 67))  # min/max latitude of the interpolation area

grd <- expand.grid(x = seq(from = x.range[1], to = x.range[2], by = 0.1),
                   y = seq(from = y.range[1], to = y.range[2], by = 0.1))  # expand points to grid
coordinates(grd) <- ~x + y
gridded(grd) <- TRUE

plot(grd, cex = 1.5, col = "grey")
points(DATA_EMIRATES, pch = 1, col = "red", cex = 1)

idw <- idw(formula = values ~ 1, locations = DATA_EMIRATES, 
           newdata = grd)  # apply idw model for the data (interpolation)

idw.output = as.data.frame(idw)  # output is defined as a data table
names(idw.output)[1:3] <- c("Lon", "Lat", "values")  # give names to the modelled variables

write.csv(idw.output, file = "C:/RICARDO-AEA/SATELLITE_STUFF/AOD_MOD04_L2_3K/AOD_MOD04_10km_03Oct2016_EMIRATES_Interp.csv", row.names=FALSE)


############################################################################################
############################################################################################
############################################################################################
############################################################################################

# Load interpolated data at resolution of 10km

# DATA_EMIRATES_interp <- read.csv("C:/RICARDO-AEA/SATELLITE_STUFF/AOD_MOD04_L2_3K/MOD04_AOD_2016_103/DATA_EMIRATES_Interp.csv")
# DATA_EMIRATES_interp <- read.csv("C:/SATELLITE_STUFF/AOD_MOD04_L2_3K/MOD04_AOD_2016_103/DATA_EMIRATES_Interp.csv")

# Load interpolated data at resolution of 3km 
# DATA_EMIRATES_interp <- read.csv("E:/Ricardo-AEA/SATELLITE_STUFF/AOD_MOD04_L2_3K/MOD04_AOD_3K_2016_124/DATA_EMIRATES_Interp_3K.csv")
DATA_EMIRATES_interp <- read.csv("E:/Ricardo-AEA/SATELLITE_STUFF/AOD_MOD04_L2_3K/MOD04_AOD_10K_2016_277_3Oct2016/AOD_MOD04_10km_03Oct2016_EMIRATES.csv_Interp.csv")


# make raster with interpolated data

coordinates(DATA_EMIRATES_interp) <- ~ Lon + Lat
# coerce to SpatialPixelsDataFrame
gridded(DATA_EMIRATES_interp) <- TRUE
raster_DATA_EMIRATES_interp <- raster(DATA_EMIRATES_interp)
projection(raster_DATA_EMIRATES_interp) <- CRS("+proj=longlat +datum=WGS84")
plot(raster_DATA_EMIRATES_interp)
# mapview(raster_DATA_EMIRATES_interp)

writeRaster(raster_DATA_EMIRATES_interp, "Data_Emirates.tif", overwrite = TRUE)

Data_Emirates_tif <- raster("Data_Emirates.tif")
plot(Data_Emirates_tif)

MIN_data <- min(minValue(Data_Emirates_tif))
MAX_data <- max(maxValue(Data_Emirates_tif))    


## create an unique legend and colorbar for Data Emirates

rast_pal_EMIRATES <- colorNumeric(c("#ffffff", "#b7b700", "#e50000"), 
                              c(MIN_data,MAX_data),
                              na.color = "transparent")


  map <- leaflet() %>% 
    setView(45, 25, 5) %>%
    addTiles(group = "OSM (default)") %>%
    addProviderTiles("OpenStreetMap.Mapnik", group = "Road map") %>%
    addProviderTiles("Thunderforest.Landscape", group = "Topographical") %>%
    addProviderTiles("Esri.WorldImagery", group = "Satellite") %>%
    addProviderTiles("Stamen.TonerLite", group = "Toner Lite") %>%
    addRasterImage(Data_Emirates_tif, 
                   colors = rast_pal_EMIRATES, 
                   opacity = 0.6,
                   group = "AOD_EMIRATES") %>%
  addLegend("bottomright", pal = rast_pal_EMIRATES, values = c(MIN_data, MAX_data),
        #    title = "<br><strong>AOD from MODIS (10km) 12 April 2016 </strong>",
            title = "<br><strong>AOD from MODIS (3km) 03 may 2016 </strong>",
            labFormat = labelFormat(prefix = ""),
            opacity = 0.6) %>%
    addLayersControl(
      baseGroups = c("Road map", "Topographical", "Satellite", "Toner Lite"),
      overlayGroups = c("AOD_EMIRATES"),
      options = layersControlOptions(collapsed = TRUE)) 
  map
  
  
# save map in html format
#   saveWidget(map,
#              file="AOD_EMIRATES_2016_04_12_3km.html",
#              selfcontained = FALSE)


  # save map in html format
  saveWidget(map,
             file="AOD_EMIRATES_2016_05_03_3km.html",
             selfcontained = FALSE)
  
  
  