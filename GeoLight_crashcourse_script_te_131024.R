###########################################################
###  Geolocator data analysis script (workshop 2013)    ###
###########################################################
# 
# Author(s)          : tamara.emmenegger@vogelwarte.ch
# Data               : GDL data (.gle files)
#
# start -------------------------------------------------------------------

# 1) Prep -------------------------

# Packages
install.packages("GeoLight")
library(GeoLight)

# debugged functions (preliminary beta version!)
# HillEkstromCalib function
HillEkstromCalib <- function (tFirst, tSecond, type, site, start.angle = -6, distanceFilter = FALSE, 
                              distance, plot = TRUE) {
  tFirst <- as.POSIXct(as.character(tFirst), "UTC")
  tSecond <- as.POSIXct(as.character(tSecond), "UTC")
  sites <- as.numeric(length(levels(as.factor(site[as.numeric(site) != 0]))))
  HECalib <- rep(NA, sites)
  for (j in 1:sites) {
    b <- 0
    start <- start.angle
    repeat {
      if (start - ((b * 0.1) - 0.1) < -9) {
        HECalib[j] <- NA
        break
      }
      t0 <- var(na.omit(coord(tFirst[site == j], tSecond[site == 
                                                           j], type[site == j], start - ((b * 0.1) - 0.1), 
                              note = F)[, 2]))
      t1 <- var(na.omit(coord(tFirst[site == j], tSecond[site == 
                                                           j], type[site == j], start - (b * 0.1), note = F)[, 
                                                                                                             2]))
      t2 <- var(na.omit(coord(tFirst[site == j], tSecond[site == 
                                                           j], type[site == j], start - ((b * 0.1) + 0.1), 
                              note = F)[, 2]))
      if (sum(is.na(c(t0, t1, t2))) > 0) {
        HECalib[j] <- NA
        break
      }
      if (t0 > t1 & t1 < t2) {
        HECalib[j] <- start - (b * 0.1)
        break
      }
      if (start - ((b * 0.1) - 0.1) > 9) {
        HECalib[j] <- NA
        break
      }
      f0 <- var(na.omit(coord(tFirst[site == j], tSecond[site == 
                                                           j], type[site == j], start + ((b * 0.1) - 0.1), 
                              note = F)[, 2]))
      f1 <- var(na.omit(coord(tFirst[site == j], tSecond[site == 
                                                           j], type[site == j], start + (b * 0.1), note = F)[, 
                                                                                                             2]))
      f2 <- var(na.omit(coord(tFirst[site == j], tSecond[site == 
                                                           j], type[site == j], start + ((b * 0.1) + 0.1), 
                              note = F)[, 2]))
      if (sum(is.na(c(f0, f1, f2))) > 0) {
        HECalib[j] <- NA
        break
      }
      if (f0 > f1 & f1 < f2) {
        HECalib[j] <- start + (b * 0.1)
        break
      }
      b <- b + 1
    }
  }
  layout(matrix(seq(1, sites * 2), ncol = 2, nrow = sites, 
                byrow = T))
  par(oma = c(3, 3, 3, 2))
  for (j in 1:sites) {
    if (is.na(HECalib[j])) {
      par(mar = c(1, 1, 1, 1), bty = "n")
      plot(0, 0, cex = 0, pch = 20, col = "white", ylab = "", 
           xlab = "", xaxt = "n", yaxt = "n",cex.axis = 1)
      text(0, 0, "NA", cex = 1)
      par(mar = c(1, 1, 1, 1), bty = "n")
      plot(0, 0, cex = 0, pch = 20, col = "white", ylab = "", 
           xlab = "", xaxt = "n", yaxt = "n",cex.axis = 1)
    }
    else {
      angles <- c(seq(HECalib[j] - 2, HECalib[j] + 2, 0.2))
      latM <- matrix(ncol = length(angles), nrow = length(tFirst[site == 
                                                                   j]))
      for (i in 1:ncol(latM)) {
        latM[, i] <- coord(tFirst[site == j], tSecond[site == 
                                                        j], type[site == j], c(angles[i]), note = F)[, 
                                                                                                     2]
      }
      latT <- latM
      var1 <- rep(NA, ncol(latT))
      n1 <- rep(NA, ncol(latT))
      min <- rep(NA, ncol(latT))
      max <- rep(NA, ncol(latT))
      for (t in 1:length(var1)) {
        var1[t] <- var(na.omit(latT[, t]))
        n1[t] <- length(na.omit(latT[, t]))
        min[t] <- if (length(na.omit(latT[, t])) <= 1) 
          NA
        else min(na.omit(latT[, t]))
        max[t] <- if (length(na.omit(latT[, t])) <= 1) 
          NA
        else max(na.omit(latT[, t]))
      }
      colors <- grey.colors(length(angles))
      par(mar = c(1, 1, 1, 1))
      plot(tFirst[site == j], latT[, 1],cex.axis = 1, ylim = c(min(na.omit(min)), 
                                                  max(na.omit(max))), type = "o", cex = 0.5, pch = 20, 
           col = colors[1], ylab = "", xlab = "")
      for (p in 2:ncol(latT)) {
        lines(tFirst[site == j], latM[, p], type = "o", 
              cex = 0.5, pch = 20, col = colors[p])
      }
      lines(tFirst[site == j], coord(tFirst[site == j], 
                                     tSecond[site == j], type[site == j], HECalib[j], 
                                     note = F)[, 2], col = "tomato2", type = "o", 
            lwd = 2, cex = 1, pch = 19)
      if (j == sites) 
        mtext("Latitude", side = 2, line = 3)
      if (j == sites) 
        mtext("Date", side = 1, line = 2.8)
      par(mar = c(1, 1, 1, 1))
      plot(angles, var1, type = "o", cex = 1, pch = 20, 
           ylab = "",cex.axis = 1)
      lines(angles, var1, type = "p", cex = 0.5, pch = 7)
      abline(v = HECalib[j], lty = 2, col = "red", lwd = 1.5)
      par(new = T)
      plot(angles, n1, type = "o", xaxt = "n", yaxt = "n", 
           col = "blue", pch = 20, ylab = "",cex.axis = 1)
      points(angles, n1, type = "p", cex = 0.5, pch = 8, 
             col = "blue")
      axis(4, cex.axis = 1)
      if (j == sites) 
        mtext("Sun elevation angles", side = 1, line = 2.8)
      legend("topright", c(paste(HECalib[j], " degrees", 
                                 sep = ""), "sample size", "variance"), pch = c(-1, 
                                                                                20, 20), lty = c(2, 2, 2), lwd = c(1.5, 0, 0), 
             col = c("red", "blue", "black"), bg = "White")
    }
  }
  mtext("Hill-Ekstrom Calibration", adj=0,cex = 1.1, line = 0.8, 
        outer = T)
  return(HECalib)
}
# siteMap function (preliminary beta version!)
siteMap <- function(coord,site,points=TRUE,map.range=c("EuroAfrica","AustralAsia","America","World"),
                    xlim=NULL,ylim=NULL,xlab="Longitude",ylab="Latitude",lwd=1,lty=1,pch=1,cex=1,
                    col="black",main=NULL,add=FALSE) {
  
  nr.sites <- length(levels(as.factor(site[site!=0])))
  site <- as.factor(site)
  
  
  if(sum(map.range==c("EuroAfrica","AustralAsia","America","World"))==4) {
    range <- c(-180,180,-75,90) 
  } else
  {
    if(map.range %in% c("EuroAfrica","AustralAsia","America","World")) {
      if(sum(map.range=="EuroAfrica") ==1)    range <- c(-24,55,-55,70)
      if(sum(map.range=="AustralAsia") ==1)   range <- c(60,190,-55,78) 
      if(sum(map.range=="America") ==1)       range <- c(-170,-20,-65,78)
      if(sum(map.range=="World") ==1)        range <- c(-180,180,-75,90)
    } else {
      range <- c(-180,180,-75,90) 
    }
  }
  
  if(length(xlim)==0) range[1:2] <- range[1:2] else range[1:2] <- xlim
  if(length(ylim)==0) range[3:4] <- range[3:4] else range[3:4] <- ylim
  
  if(length(main)==0) main <- "" else main <- main
  
  
  # Colors
  colors <- rainbow(nr.sites)
  
  
  if(add==FALSE) {         
    par(oma=c(3,3,1,0.5))
    map(xlim=c(range[1],range[2]),ylim=c(range[3],range[4]),fill=T,lwd=0.01,col=c("grey90"),add=F,mar=c(rep(0.5,4)))
    map(xlim=c(range[1],range[2]),ylim=c(range[3],range[4]),interior=TRUE,col=c("darkgrey"),add=TRUE)
    mtext(xlab,side=1,line=2.2,font=3)
    mtext(ylab,side=2,line=2.5,font=3)
    map.axes()
    
    mtext(main,line=0.6,cex=1.2)
  }
  
  if(length(colors)==length(levels(site))) col <- colors else col <- colors[1:length(levels(site))]
  
  sites <- length(levels(site[site!=0]))
  site2<-unique(site)
  
  if(points==TRUE){
    for(i in 1:sites){
      points(coord[site==site2[i],],cex=cex,pch=pch,col=col[i])
    }
  }
  
  for(j in 1:sites){
    X <- na.omit(coord[site==site2[j],])
    hpts <- chull(na.omit(coord[site==site2[j],]))
    hpts <- c(hpts,hpts[1])
    lines(X[hpts,],lty=lty,lwd=lwd,col=col[j])
  }
  
  legend("bottomright",c(letters[1:nr.sites]),pch=pch, col=colors[1:nr.sites])
  
}

# Working directory
# setwd("path of filename")             ### change !!
setwd("C:\\Users\\tamarlig\\Dropbox\\homeoffice\\canmove\\04-geolight")

# 2) Read-in data -------------------------

# .gle file
filename <- "example_data.gle"     ### change !!

tab <- gleTrans(filename)
tab$tFirst <- as.POSIXct(as.character(tab$tFirst),"UTC")
tab$tSecond <- as.POSIXct(as.character(tab$tSecond),"UTC")
head(tab)

# 3) Analysis -------------------------

# Filters
k=2                                     # number of interquartile ranges before saying that a twilight event is an outlier
toFilter <- loessFilter(tab$tFirst,tab$tSecond,tab$type,k,plot=TRUE)
tab <- tab[toFilter,]

# Stationary periods
days=21                                 # lower threshold length of stationary period
quantile=0.95                           # probability threshold for stationary site selection
periods <- changeLight(tab$tFirst,tab$tSecond,tab$type,quantile,days,plot=T)  

# In-habitat calibration
known.coord=c(7.216667, 46.18333)             ### change !!
inHabitatCalib <- getElevation(tab$tFirst[periods$site==1],tab$tSecond[periods$site==1],
                              tab$type[periods$site==1],known.coord,plot=F) 

#windows()
start.angle=-5                          # initial value for sun elevation angle in the iterative process
HECalib <- HillEkstromCalib(tab$tFirst,tab$tSecond,tab$type,periods$site,start.angle,plot=T)

# Positioning
coord <- coord(tab$tFirst,tab$tSecond,tab$type, -6)             ### change !!
 

# 4) Visualization -------------------------
#windows()
par(mfrow=c(1,1),mar=c(0.1,0.1,0.1,0.1))
siteMap(coord,periods$site,points=T,map.range="EuroAfrica",
        xlim=c(-20,30),ylim=c(0,50),xlab="Longitude",
        ylab="Latitude",lwd=1,lty=0,pch=20,cex=1,add=FALSE)


# 5) Read-out -------------------------


# end ---------------------------------------------------------------------
