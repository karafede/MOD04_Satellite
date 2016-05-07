
# install.packages("devtools")
 #library(devtools)
# devtools::install_github("hadley/httr")

library("httr")
library(RCurl)
library(stringr)
library(plyr)
library(curl)

# setInternet2(TRUE)
# options('download.file.method'='curl')
# options(RCurlOptions = list(cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl")))   

# setwd("C:/RICARDO-AEA/AOD_MOD04_L2_3K")
setwd("C:/SATELLITE_STUFF/AOD_MOD04_L2_3K")
setwd("E:/Ricardo-AEA/SATELLITE_STUFF/AOD_MOD04_L2_3K")

# url = "ftp://disc2.nascom.nasa.gov/data/TRMM/Gridded/Derived_Products/3B42_V6/Daily/2009/" # this works
# url = "ftp://ladsftp.nascom.nasa.gov/allData/6/MOD04_L2/2016/090/"  # this works
url = 'ftp://nrt3.modaps.eosdis.nasa.gov/allData/6/MOD04_3K/2016/105/' # 3km resolution
# download data for April 12 2016
url = 'ftp://karafede:Password07@nrt3.modaps.eosdis.nasa.gov/allData/6/MOD04_L2/2016/103/' #10km resolution, 103 rd day of the year
url = 'ftp://karafede:Password07@nrt3.modaps.eosdis.nasa.gov/allData/6/MOD04_3K/2016/124/' #3km resolution, 124 th day of the year

filenames = getURL(url, ftp.use.epsv = FALSE, ftplistonly = TRUE, crlf = TRUE) 
filenames = paste(url, strsplit(filenames, "\r*\n")[[1]], sep = "") 
filenames_hdf <- unlist(str_extract_all(filenames, ".+(.hdf$)"))
# mapply(download.file, filenames, , method="curl") 
# start downloading data in the main directory
mapply(download.file, filenames_hdf,basename(filenames_hdf), method = 'curl') 

# move file into the directory "MOD04_AOD_2016_103"
origindir <- c("C:/SATELLITE_STUFF/AOD_MOD04_L2_3K")
targetdir <- c("C:/SATELLITE_STUFF/AOD_MOD04_L2_3K/MOD04_AOD_2016_103")

origindir <- c("E:/Ricardo-AEA/SATELLITE_STUFF/AOD_MOD04_L2_3K")
targetdir <- c("E:/Ricardo-AEA/SATELLITE_STUFF/AOD_MOD04_L2_3K/MOD04_AOD_3K_2016_124")


filenames = getURL(url, ftp.use.epsv = FALSE, ftplistonly = TRUE, crlf = TRUE) 
filenames_trunc <- strsplit(filenames, "\r*\n")[[1]]
filenames_trunc_hdf <- unlist(str_extract_all(filenames_trunc, ".+(.hdf$)"))
filestocopy <- filenames_trunc_hdf

file.copy(from = filestocopy, to = targetdir, 
          overwrite = TRUE, recursive = FALSE, 
          copy.mode = TRUE)
file.remove(filestocopy)

######################################################################################
######################################################################################
######################################################################################
########### OLD STUFF ###########################################
###### do not run stuff below ###################################

# FTP 
# Download the files within a directory. 
# ftp = 'ftp://ladsftp.nascom.nasa.gov/allData/6/MOD04_L2/2016/090/'
ftp = 'ftp://karafede:Password07@nrt3.modaps.eosdis.nasa.gov/allData/6/MOD04_L2/2016/103/' 
filenames = getURL(ftp, ftp.use.epsv = FALSE, dirlistonly = TRUE) 


# filenames = paste(ftp, strsplit(filenames, "\r*\n")[[1]], sep = "") 
# or
filenames <- strsplit(filenames, "\r*\n")[[1]]
filenames_hdf <- unlist(str_extract_all(filenames, ".+(.hdf)"))
filenames_hdf <- as.data.frame(filenames_hdf)
str(filenames_hdf)


dl <- lapply(filenames_hdf, function(x) {
  path <- ftp
  dest <- 'prova'
  try(download.file(path, dest))
})


try(download.file(path, dest))

setInternet2(use = TRUE)
download.file(url='ftp://karafede:Password07@nrt3.modaps.eosdis.nasa.gov/allData/6/MOD04_L2/2016/103/',
              destfile='localfile.zip', method='auto')

###  This works ! #################################################################################
# url = "ftp://disc2.nascom.nasa.gov/data/TRMM/Gridded/Derived_Products/3B42_V6/Daily/2009/" # this works
# url = "ftp://ladsftp.nascom.nasa.gov/allData/6/MOD04_L2/2016/090/"  # this works
# url = 'ftp://karafede:Password07@nrt3.modaps.eosdis.nasa.gov/allData/6/MOD04_L2/2016/104/'
filenames = getURL(url, ftp.use.epsv = FALSE, ftplistonly = TRUE, crlf = TRUE) 
filenames = paste(url, strsplit(filenames, "\r*\n")[[1]], sep = "") 
mapply(download.file, filenames, basename(filenames), method="curl") 
####################################################################################################










con = getCurlHandle( ftp.use.epsv = FALSE) 

ftp_files = getURL(ftp, dirlistonly = TRUE) 
filenames = paste(ftp, strsplit(filenames, "\r*\n")[[1]])
# filenames <- strsplit(ftp_files, "\r*\n")[[1]]
filenames_hdf <- unlist(str_extract_all(filenames, ".+(.hdf)"))







download.file(url='https://s3.amazonaws.com/tripdata/201307-citibike-tripdata.zip',
              destfile='localfile.zip', method='auto')



# url = 'ftp://ladsftp.nascom.nasa.gov/allData/6/MOD04_L2/2016/090/' 
# url = 'ftp://nrt3.modaps.eosdis.nasa.gov/allData/6/MOD04_L2/2016/103/'
ftp = 'ftp://karafede:Password07@nrt3.modaps.eosdis.nasa.gov/allData/6/MOD04_L2/2016/103/' 
# filenames = getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE) 
ftp_files = getURL(ftp, dirlistonly = TRUE) 


# Deal with newlines as \n or \r\n. (BDR) 
# Or alternatively, instruct libcurl to change \n's to \r\n's for us 
# with crlf = TRUE 
#  filenames = getURL(url, ftp.use.epsv = FALSE, ftplistonly = TRUE, 
#  crlf = TRUE) 

# filenames = paste(url, strsplit(filenames, "\r*\n")[[1]], sep = "") 
filenames <- strsplit(ftp_files, "\r*\n")[[1]]
filenames_hdf <- unlist(str_extract_all(filenames, ".+(.hdf)"))
filenames_hdf[1:3]
filenames_list <- as.list(filenames_hdf)

downloadFTP <- function(filename, folder, handle) {
  dir.create(folder, showWarnings = FALSE)
  fileurl <- str_c (ftp, filename)
  if (!file.exists(str_c (folder, "/", filename))) {
    content <- try(getURLContent(fileurl, curl = handle))
    write(content, str_c (folder, "/", filename))
    Sys.sleep(1)
  } 
} 

HEAD("http://google.com/")
handle_find("http://google.com/")
handle_reset("http://google.com/")
HEAD("http://google.com/")

HEAD("ftp://karafede:Password07@nrt3.modaps.eosdis.nasa.gov/allData/6/MOD04_L2/2016/103/")
AAA <- handle_find("ftp://karafede:Password07@nrt3.modaps.eosdis.nasa.gov/allData/6/MOD04_L2/2016/103/")
AAA <- handle_reset("ftp://karafede:Password07@nrt3.modaps.eosdis.nasa.gov/allData/6/MOD04_L2/2016/103/")
HEAD("ftp://karafede:Password07@nrt3.modaps.eosdis.nasa.gov/allData/6/MOD04_L2/2016/103/")

handle_a = getCurlHandle( ftp.use.epsv = TRUE)

## download files
l_ply(filenames_hdf, downloadFTP, 
       folder = "cran_task",
       handle = handle_a)



con = getCurlHandle( ftp.use.epsv = FALSE)

filenames[1]


contents = sapply(filenames[1:length(filenames)], function(x) 
  try(RCurl::getURL(x, curl = con)))   
names(contents) = filenames[1:length(contents)] 
# }

#extracting the filenames from the names(contents) 
#Used to extract the exact file names from the URL 

file_name <- unlist(lapply(names(contents),function(i){substring(i,59)})) 
file_name

#creating the file names 

physicalFileName <- 
  unlist(lapply(seq_along(1:length(file_name)),function(i){paste("./MOD04_AOD_2016_90/",file_name[i],sep="",collapse="")})) 
physicalFileName[1] 

for(i in 1:length(contents)) 
{ 
  file.create(physicalFileName[i]) 
  write(contents[i],physicalFileName[i]) 
} 
