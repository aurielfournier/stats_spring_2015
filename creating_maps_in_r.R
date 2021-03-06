# maps in R from Quan Tech Class
# 2 19 2015
# https://github.com/Robinlovelace/Creating-maps-in-R
# 


x <- c("ggplot2","spdep","ggmap", "rgdal", "rgeos", "maptools", "dplyr", "tidyr") # the packages
#install.packages(x) # warning: this may take a number of minutes
lapply(x, library, character.only = TRUE) # load the required packages


library(rgdal)
lnd_sport <- readOGR(dsn = "data", layer = "london_sport")
#dsn = folder of where the data is (can also be connection to a database)
#layer in the geospatial database or the shapefile (ours is a shapefile)

lnd_sport@polygons
  
  # the @ symbol gives you access to the slots
  # there is a data, polygons, coordinates and other slots
  # hit tab after the @ symbol to see what is avaiable

lnd_sport@data[lnd_sport$Partic_Per<15,]

sel <- lnd_sport$Partic_Per>25

plot(lnd_sport[sel,],col="lightblue", add=T)

EPSG <- make_EPSG()
head(EPSG)


lnd_sport@proj4string

# to convert we would need the EPSG code for the new projection

EPSG[grep("WGS 84$", EPSG$note),]

# so 4326 is the code for WGS 84

#now to reproject it

lnd_sport_wgs84 <- spTransform(lnd_sport, CRS("+init=epsg:4326"))
lnd_sport_wgs84@bbox
# see now it is in lat long

crime_data <- read.csv("data/mps-recordedcrime-borough.csv")
summary(crime_data$CrimeType)

crime_theft <- crime_data[crime_data$CrimeType=="Theft & Handling",]
crime_ag <- aggregate(CrimeCount ~ Borough, FUN=sum, data=crime_theft)

  # do these names match?
lnd_sport$name %in% crime_ag$Borough

crime_ag$Borough
# there is a null value!! not good (City of London is missing)

#page 11 of https://github.com/Robinlovelace/Creating-maps-in-R/blob/master/intro-spatial-rl.pdf
