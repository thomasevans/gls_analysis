# Primarily written by tamara.emmenegger@vogelwarte.ch from a workshop at Lund University in October 2013.
# Adapated by Tom Evans (Lund University): thomas.jude.evans@gmail.com




# siteMap function (preliminary beta version!)
siteMap2 <- function(coord,site,points=TRUE,map.range=c("EuroAfrica","AustralAsia","America","World"),
                     xlim=NULL,ylim=NULL,xlab="Longitude",ylab="Latitude",lwd=1,lty=1,pch=1,cex=1,
                     col="black",main=NULL,add=FALSE) {
  require("maps")
  
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
  #   levels(as.factor(periods$site))
  if(length(colors)==length(levels(site))) col <- colors else col <- colors[1:length(levels(site))]
  
  sites <- length(levels(site[site != 0]))
  site2 <- unique(site[site != 0])
  #   site2 <- site2[2:length(site2)]   # Test
  col <-  rainbow(length(levels(as.factor(site)))-1)  #test
  
  
  if(points == TRUE){
    for(i in 1:sites){
      points(coord[site == site2[i],], cex = cex, pch = pch, col = col[i])
    }
  }
  
  for(j in 1:sites){
    X <- na.omit(coord[site == site2[j],])
    hpts <- chull(na.omit(coord[site == site2[j],]))
    hpts <- c(hpts,hpts[1])
    lines(X[hpts,],lty = lty, lwd = lwd, col = col[j])
  }
  
  legend("bottomleft",c(letters[1:nr.sites]), pch = pch, col = colors[1:nr.sites])
  
}