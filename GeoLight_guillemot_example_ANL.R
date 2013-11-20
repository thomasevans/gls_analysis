# Written by Tom Evans <thomas.jude.evans@gmail.com>, but based on code from
# tamara.emmenegger@vogelwarte.ch from a workshop at Lund University in October 2013.

# Here I use the GeoLight package to analyse geolocator data from a common guillemot using the GeoLight R package.

# 1. Packages --------
install.packages(c("GeoLight", "maps", "adehabitatHR"))
library(GeoLight)

# 2. Set working-directory --------
setwd("D:/Dropbox/R_projects/gls_analysis/examples/aak970_2009_2011")
warning("Change to your R working directory - e.g. where you have stored the geolocator data.")

# 3. Read in data  --------
gls_data <- read.delim("20110615_aak970_gls_download_000_thresh10.trn", sep = ",")
names(gls_data)   <-  c("datetime", "type", "certainty")

# Convert to datetime object
gls_data$datetime     <-    as.POSIXct(as.character(gls_data$datetime), format = "%d/%m/%y %H:%M", "UTC")

# 4. Get data into GeoLight format   -----

# Set up empty vectors for transition pairs
tFirst <- NULL
tSecond <- NULL
ttype <- NULL


n <- 1   # Keep count of number of transition pairs
# i <- 2    # For testing the for loop below

for(i in 1:(length(gls_data$type)-1)){
  
  # Calculate time difference between consecutive transitions
  t.dif <- difftime(gls_data$datetime[i+1], gls_data$datetime[i], 
                    units = "mins")
  t.dif <- as.numeric(t.dif)
  
  # If time difference is < 24h and consecutive transitions are of
  # different type (one of sunset + one of sunrise) add to dataframe
  if(t.dif < 1440  &  gls_data$type[i] != gls_data$type[i+1] ){
    
    # First transtion time
    tFirst[n]     <-  as.character(gls_data$datetime[i])
    # Second transition time
    tSecond[n]    <-  as.character(gls_data$datetime[i+1])
    
    # Type of transition pair
    if(gls_data$type[i] == "Sunrise") {ttype[n] <- 1} else{
      ttype[n] <- 2
    }

    n <- n + 1
    
  }
}

# Make a dataframe of transitions
tab <- cbind(tFirst,tSecond,ttype)
tab <- as.data.frame(tab)
names(tab) <- c("tFirst", "tSecond", "type")

# Correct format (e.g. get datetimes as datetime, not character)
tab$tFirst     <-    as.POSIXct(tab$tFirst, "UTC")
tab$tSecond    <-    as.POSIXct(tab$tSecond, "UTC")
tab$type       <-    as.numeric(tab$type)
# Remove row numbers
row.names(tab) <-    NULL

# Check that this looks ok, view first few rows
head(tab)


# 5. Get subset of data (year 2010-11) -----
tab.original <- tab

tab  <- tab.original[700:1200,]      # Year 2, 2010-2011 - excludes most of breeding period
# tab  <- tab.original[50:530,]        # Year 1, 2009-2010 - excludes most of breeding period

# # view dates for first and final days in dataframe 'tab'
# tab.original[700,]
# tab.original[1200,]
# 
# tab.original[50,]
# tab.original[530,]


# 6. Loess filter -------------------------

# Filters

# Interquartile ranges to accept
# Here we set quite a low k value - i.e. we have a low tolerance for sunrises or sunset that fall far the loess line
k = 0.5  

# Run loessFilter function
toFilter <- loessFilter(tab$tFirst, tab$tSecond, tab$type, k = k, plot = TRUE)

# Remove outliers
tab <- tab[toFilter,]

# ?loessFilter     # Lookup help for function


# 7. Identify Stationary periods -----

# Paramaters 'quantile' and 'days' can be varied and will result in differing numbers of 'stationary periods' being identified.
# Also look at help, there are some futher paramaters.

#?changeLight

periods <- changeLight(tab$tFirst, tab$tSecond, tab$type,
                       , plot = TRUE, quantile = 0.97,
                       days = 15, rise.prob = NA) 


# 8. Calibration/ sun elevation calculation ------
# There are two methods provided, 'Hill-Ekstrom' and 'In habitat', we try both here.

# In-habitat calibration
# Run for a stationary period where your organism is at a known location (e.g. breeding site)

karlso = c(17.958419, 57.287751)     # Stora Karlsö

inHabitatCalib <- getElevation(tab$tFirst[periods$site == 1],
                               tab$tSecond[periods$site == 1],
                               tab$type[periods$site == 1],
                               known.coord = karlso, plot = TRUE) 

# Hill-Ekstrom calibration
# Important here to choose a starting angle which is close to expected value. The function runs as a itterative process, finding a local 'minima' in latitudinal variance. If you start too far from the 'true' value then you may get a wrong elevation value. I would suggest trying with a few different starting values in range -2:-7, where -2 represents high levels of shading, and -6 near 'perfect' conditions with no shading.

# Tamarra suggested that this (HillEkstromCalib2) beta version Hill-Ekstrom calibration works better than that currently in the package (HillEkstromCalib)
source("hillekstrom2.R")

HECalib <- HillEkstromCalib2(tab$tFirst, tab$tSecond, tab$type,
                            periods$site,  plot = FALSE,
                            distanceFilter = TRUE, distance = 5,
                            start.angle = -4)

# View calculated sun elevations
HECalib

# If you want to change the value for one period, e.g. if you use in habitat calibration rather than Hill-Ekstrom for breeding area
# HECalib[1] <- -4

# Sometimes you get NA (missing values) values for some areas, I think when there are too few points to reliably estimate elevation.
# Here we replace missing values with some value (-4 here)
#  HECalib[is.na(HECalib)] <- -4




# 9. Calculate coordinates -------

# Positioning - comparing different sun-elevation values

# Different elevation value for each 'stationary' period
coords <- coord(tab$tFirst,tab$tSecond,tab$type, HECalib)

# Fixed elevation value, the mean across that calculated for all sites
coords.fixed <- coord(tab$tFirst,tab$tSecond,tab$type, mean(HECalib))

# Some 'wrong' value for elevation value
coords.fixedwrong <- coord(tab$tFirst,tab$tSecond,tab$type, mean(HECalib) - 2)            



# 10. Map data - points ----------

# Make sure that 'sitemap2.R' is in your working directory
source("sitemap2.R")

# Vairable sun elevation angle
par(mfrow=c(1,1),mar=c(0.1,0.1,0.1,0.1))
siteMap2(coords, periods$site, points = TRUE,
         xlim = c(14,23), ylim = c(53,61), xlab = "Longitude",
         ylab = "Latitude", lwd = 1, lty = 0, pch = 20, cex = 1, add = FALSE)

# A fixed sun elevation angle
par(mfrow=c(1,1),mar=c(0.1,0.1,0.1,0.1))
siteMap2(coords.fixed, periods$site, points = TRUE,
        xlim = c(14,23), ylim = c(53,61), xlab = "Longitude",
        ylab = "Latitude", lwd = 1, lty = 0, pch = 20, cex = 1, add = FALSE)

# A wrong elevation value
par(mfrow=c(1,1),mar=c(0.1,0.1,0.1,0.1))
siteMap2(coords.fixedwrong, periods$site, points = TRUE,
         xlim = c(14,23), ylim = c(53,61), xlab = "Longitude",
         ylab = "Latitude", lwd = 1, lty = 0, pch = 20, cex = 1, add = FALSE)




# 11a. Map data - Kernels ------

library("adehabitatHR")

par(mfrow=c(1,1),mar=c(0.1,0.1,0.1,0.1))

# 'empty' map
siteMap2(coords, periods$site, points = FALSE,
         xlim = c(12,23), ylim = c(52,62), xlab = "Longitude",
         ylab = "Latitude", lwd = 1, lty = 0, pch = 20, cex = 1, add = FALSE)

# Get a vector of some colours
cols <- rainbow(length(levels(as.factor(periods$site)))-1)

# For each period plot a kernels
for(i in 1:length(unique(periods$site))){
  # A filter - only data from this period
  f <- !is.na(coords[,1]) & !is.na(coords[,2])  & periods$site == i #time period....
  
  # Get in format needed for adehabitatHR
  xy <- t(rbind(coords[f,1],coords[f,2]))
  xy <- as.data.frame(xy)
  names(xy) <- c("x", "y")
  coordinates(xy) <- c("x","y")
  
  # Make kernel object
  kud <- kernelUD(xy)
  
  # Plot different kerel density bands (5, 25, 50 %)
  plot(getverticeshr(kud, 5), col = NA, add = TRUE, lwd = 2,  border = cols[i])
  plot(getverticeshr(kud, 25), col = NA, add = TRUE, lty = 4, lwd = 2,  border = cols[i])
  plot(getverticeshr(kud, 50), col = NA, add = TRUE, lty = 3, lwd = 1,  border = cols[i])

}




# 11b. Map data - Kernels - no equinox period ------

library("adehabitatHR")

par(mfrow=c(1,1),mar=c(0.1,0.1,0.1,0.1))

# 'empty' map
siteMap2(coords, periods$site, points = FALSE,
         xlim = c(12,23), ylim = c(52,62), xlab = "Longitude",
         ylab = "Latitude", lwd = 1, lty = 0, pch = 20, cex = 1, add = FALSE)

# Get a vector of some colours
cols <- rainbow(length(levels(as.factor(periods$site)))-1)

# For each period plot a kernels
for(i in c(1,2,4)){
  # A filter - only data from this period
  f <- !is.na(coords[,1]) & !is.na(coords[,2])  & periods$site == i #time period....
  
  # Get in format needed for adehabitatHR
  xy <- t(rbind(coords[f,1],coords[f,2]))
  xy <- as.data.frame(xy)
  names(xy) <- c("x", "y")
  coordinates(xy) <- c("x","y")
  
  # Make kernel object
  kud <- kernelUD(xy)
  
  # Plot different kerel density bands (5, 25, 50 %)
  plot(getverticeshr(kud, 5), col = NA, add = TRUE, lwd = 2,  border = cols[i])
  plot(getverticeshr(kud, 25), col = NA, add = TRUE, lty = 4, lwd = 2,  border = cols[i])
  plot(getverticeshr(kud, 50), col = NA, add = TRUE, lty = 3, lwd = 1,  border = cols[i])
  
}