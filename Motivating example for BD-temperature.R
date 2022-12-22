##############################################################################-#
# Motivating example for studying the BD-temperature ERRs in the original      #
# scale. TS-based strategy + bs-based strategy                                 #
##############################################################################-#

library(splines)
library(dlnm)
library(mvmeta)
library(sf)
source("fun for TS.R")

################################################################################
################### load data ##################################################
################################################################################
load("data_by_province.Rdata")
data_all <- do.call(rbind,datalist)
regions <- names(datalist)
m <- length(regions)
probs <- c(10,30,60,90)/100
knots <- quantile(data_all$tm,probs = probs)
bound <- range(data_all$tm)

load("mapdata_China.Rdata")
data_province <- sapply(datalist, function(x) sum(x$count))
index <- sapply(names(data_province), function(x) which(grepl(substr(x,1,2),mapdata_province$NAME)))
data_province <- data.frame(name = mapdata_province$ID[index],counts = as.vector(data_province),
                            X = mapdata_province$X[index], Y = mapdata_province$Y[index])
data_province <- st_as_sf(data_province,coords = c("X","Y"),crs = "WGS84")
data_province <- st_transform(data_province,crs = st_crs(mapdata_iland))

names(datalist) <- data_province$name
regions <- data_province$name
ranges <- t(sapply(datalist,function(x) range(x$tm,na.rm=T)))
################################################################################
############# Parameter selection based on AIC #################################
################################################################################
## df for long-term trends-----
degree <- 2
knots <- quantile(data_all$tm,probs = probs)
df <- length(knots) + degree
ymat <- matrix(NA,m,df,dimnames=list(rownames = regions))
Slist <- vector("list",m)
names(Slist) <- regions
AICs_all <- NULL
a <- 3:20
for (j in a) {
  AICs <- 0
  for (i in seq(m)) {
    cat(regions[i]," ")
    data <- datalist[[i]]
    data$time <- 1:dim(data)[1]
    data$lograin <- log(data$rain)
    basvar <- trun_spline(x = data$tm,knots = knots,cenknot_lower = knots[2],cenknot_upper = knots[3],degree = degree)
    index <- apply(basvar, 2, function(x) any(x != 0))
    model <- glm(count ~ basvar[,index] + ns(humid,df = 3) + ns(sun,df = 3)+ns(time,df = j),family=quasipoisson(),data)
    AICs <- AICs+fqbic(model)
  }
  AICs_all <- c(AICs_all,AICs)
}
names(AICs_all) <- a
plot(a,AICs_all)

## degree and centered knots----
degree <- 2
knots <- quantile(data_all$tm,probs = probs)
df <- length(knots) + degree
ymat <- matrix(NA,m,df,dimnames=list(rownames = regions))
Slist <- vector("list",m)
names(Slist) <- regions
AICs <- 0
for (i in seq(m)) {
  cat(regions[i]," ")
  data <- datalist[[i]]
  data$time <- 1:dim(data)[1]
  data$lograin <- log(data$rain)
  basvar <- trun_spline(x = data$tm,knots = knots,cenknot_lower = knots[2],cenknot_upper = knots[3],degree = degree)
  index <- apply(basvar, 2, function(x) any(x != 0))
  model <- glm(count ~ basvar[,index] + humid + sun+ns(time,df = 6),family=quasipoisson(),data)
  AICs <- AICs+fqbic(model)
}
AICs

## knots ----
degree <- 2
probs <- c(10,30,60,90)/100
knots <- quantile(data_all$tm,probs = probs)
df <- length(knots) + degree
ymat <- matrix(NA,m,df,dimnames=list(rownames = regions))
Slist <- vector("list",m)
names(Slist) <- regions
AICs <- 0
for (i in seq(m)) {
  cat(regions[i]," ")
  data <- datalist[[i]]
  data$time <- 1:dim(data)[1]
  data$lograin <- log(data$rain)
  basvar <- trun_spline(x = data$tm,knots = knots,cenknot_lower = knots[2],cenknot_upper = knots[3],degree = degree)
  index <- apply(basvar, 2, function(x) any(x != 0))
  model <- glm(count ~ basvar[,index] + ns(humid,df = 3) + ns(sun,df = 3)+ns(time,df = 6),family=quasipoisson(),data)
  AICs <- AICs+fqbic(model)
}
AICs
## confounder------
degree <- 2
probs <- c(20,90)/100
knots <- quantile(data_all$tm,probs = probs)
df <- length(knots) + degree
ymat <- matrix(NA,m,df,dimnames=list(rownames = regions))
Slist <- vector("list",m)
names(Slist) <- regions
AICs <- 0
for (i in seq(m)) {
  cat(regions[i]," ")
  data <- datalist[[i]]
  data$time <- 1:dim(data)[1]
  data$lograin <- log(data$rain)
  basvar <- trun_spline(x = data$tm,knots = knots,cenknot_lower = knots[1],cenknot_upper = knots[2],degree = degree)
  index <- apply(basvar, 2, function(x) any(x != 0))
  model <- glm(count ~ basvar[,index] + humid + sun + ns(time,df = 6),family=quasipoisson(),data)
  AICs <- AICs+fqbic(model)
}
AICs


################################################################################
############### Description analysis ###########################################
################################################################################
## geographic distribution----------
library(tmap)
xx <- tm_shape(mapdata_province,ylim = c(500000,6500000)) + 
        tm_polygons(col = "grey100",border.col = "grey60",lwd = 0.8) +
        tm_shape(mapdata_iland)+
        tm_lines(col = "grey60",lwd = 1) + 
        tm_shape(data_province) + 
        tm_text("name",size = 0.35,col = "red")
tmap_save(xx,filename = "location.tiff",width = 10,height = 10,units = "cm",dpi = 300)

xx <- tm_shape(mapdata_province,ylim = c(500000,6500000)) + 
  tm_polygons(col = "grey100",border.col = "grey60",lwd = 1.5) +
  tm_shape(mapdata_iland)+
  tm_lines(col = "grey60",lwd = 1.5)
tmap_save(xx,filename = "location-magnifying.tiff",width = 10,height = 10,units = "cm",dpi = 300)


## temperature range -----------
tiff(filename = "temperature range.tiff",width = 12,height = 15,units = "cm",res = 300,pointsize = 8)
par(mar=c(5,6,4,1)+0.1)
layout(1)
plot(seq(-25,35,length=m),seq(m),type="n",yaxt="n",
     ylab="",xlab="Temperature",mgp=c(2.5,1,0))
axis(2,at=seq(m),labels=regions,las=1,cex.axis=0.7)
arrows(ranges[,1],seq(m),ranges[,2],seq(m),angle=90,length=0.05,code=3)
title(ylab="Province",mgp=c(4.5,1,0))
abline(v=knots,col = "red")
abline(v=bound,lty=2)
axis(3,at = c(bound,knots),labels = c("0%","100%","10%","30%","60%","90%"),tcl = -0.2,mgp = c(1,0.5,0))
dev.off()
################################################################################
############### ERRs in the original scale using TS ############################
################################################################################
#### ----
degree <- 2
knots <- quantile(data_all$tm,probs = probs)
df <- length(knots) + degree
ymat <- matrix(NA,m,df,dimnames=list(rownames = regions))
Slist <- vector("list",m)
names(Slist) <- regions
for (i in seq(m)) {
  cat(regions[i]," ")
  data <- datalist[[i]]
  data$time <- 1:dim(data)[1]
  data$lograin <- log(data$rain)
  basvar <- trun_spline(x = data$tm,knots = knots,cenknot_lower = knots[2],cenknot_upper = knots[3],degree = degree)
  model <- glm(count ~ basvar + humid + sun + ns(time,df = 6),family=quasipoisson(),data)
  a <- grepl("basvarb[0-9]",names(model$coefficients))
  
  a1 <- model$coefficients[a]
  ymat[i,] <- a1
  
  a2 <- vcov(model)[a,a]
  Slist[[i]] <- a2
}
fit2 <- mvmeta(ymat,Slist,method = "reml")
yfit2 <- blup(fit2,vcov = T)
qt <- qtest(fit2)
H2 <- max(1,qt$Q/qt$df)
I2 <- (H2-1)/H2
cen <- median(data_all$tm)
basvar_ts <- basvar
####  average ERR -----
tiff(filename = "Average ERR for BD.tiff",width = 12,height = 8,units = "cm",res = 300,pointsize = 8)
opar <- par(mfrow = c(1,1))
x <- seq(range(ranges)[1],range(ranges)[2],by = 1) %>% c(max(ranges))
a <- rr_ts(basis = basvar,coef = fit2$coefficients,vcov = fit2$vcov,at = x,cen = cen)
par(cex.axis=1,cex.lab=1,cex.main=1.2, mar=c(4,4,3.6,0.8))
plot(x,a$RR, type="l", lwd=1, col="red",ylim = c(0.2,4.5),
     xlab="Temperature",ylab="Risk Ratio (RR)", main=NULL,mgp = c(2.5,1,0))
redtrans <- rgb(255, 0, 0, 60, maxColorValue=255) 
polygon(c(rev(x),x),c(rev(a$RR_upper),a$RR_lower),col=redtrans, border = NA)
abline(h=1,lwd = 0.7)
for (i in seq(m)) {
  x <- seq(ranges[i,1],ranges[i,2],by = 1)
  a <- rr_ts(basis = basvar,coef = yfit2[[i]]$blup,vcov = yfit2[[i]]$vcov,at = x,cen = cen)
  lines(x,a$RR,lty = 2,lwd = 0.6)
}
par(opar)
dev.off()
#### Province-specific ERR -----
tiff(filename = "Province-specific ERR for BD.tiff",width = 15,height = 25,units = "cm",res = 300,pointsize = 8)
opar <- par(mfrow = c(8,4),cex.axis=1,cex.lab=1,cex.main=1.2, mar=c(3.8,4,3.5,0.8)-0.5)
x <- seq(range(ranges)[1],range(ranges)[2],by = 1) %>% c(max(ranges))
a <- rr_ts(basis = basvar,coef = fit2$coefficients,vcov = fit2$vcov,at = x,cen = cen)
plot(x,a$RR, type="l", lwd=1, col="red",ylim = c(0.2,6),
     xlab="Temperature",ylab="RR",cex.lab = 1.5,mgp = c(2.2,1,0))
title("Average",line = 0.5)
redtrans <- rgb(255, 0, 0, 60, maxColorValue=255) 
polygon(c(rev(x),x),c(rev(a$RR_upper),a$RR_lower),col=redtrans, border = NA)
abline(h=1)
abline(v = range(ranges),lty = 2,lwd = 0.8)
for (i in seq(m)) {
  a <- rr_ts(basis = basvar,coef = yfit2[[i]]$blup,vcov = yfit2[[i]]$vcov,at = x,cen = cen)
  plot(x,a$RR, type="l", lwd=1, col="red",ylim = c(0.2,6),
       xlab="Temperature",ylab="RR",cex.lab = 1.5,mgp = c(2.2,1,0))
  title(regions[i],line = 0.5)
  redtrans <- rgb(255, 0, 0, 60, maxColorValue=255) 
  polygon(c(rev(x),x),c(rev(a$RR_upper),a$RR_lower),col=redtrans, border = NA)
  abline(h=1)
  abline(v = c(ranges[i,]),lty = 2,lwd = 0.8)
}
par(opar)
dev.off()
################################################################################
############ ERRs in the relative scale using BS ###############################
################################################################################
####------
degree <- 2
knots_relative <- probs
df <- length(knots_relative) + degree
ymat <- matrix(NA,m,df,dimnames=list(rownames = regions))
Slist <- vector("list",m)
names(Slist) <- regions

for (i in seq(m)) {
  cat(regions[i]," ")
  data <- datalist[[i]]
  data$time <- 1:dim(data)[1]
  data$lograin <- log(data$rain)
  knotsi <- quantile(data$tm,knots_relative)
  basvar <- onebasis(data$tm,fun = "bs",knots = knotsi,degree = degree)
  model <- glm(count ~ basvar + humid + sun + ns(time,df = 6),family=quasipoisson(),data)
  a <- grepl("basvarb[0-9]",names(model$coefficients))
  
  a1 <- model$coefficients[a]
  ymat[i,] <- a1
  
  a2 <- vcov(model)[a,a]
  Slist[[i]] <- a2
}
fit1 <- mvmeta(ymat,Slist,method = "reml",control = list(maxiter = 1000))
yfit1 <- blup(fit1,vcov = T)
qt <- qtest(fit1)
H2 <- max(1,qt$Q/qt$df)
I2 <- (H2-1)/H2
#### plot the ERRs in the two scales -------
options(warn=-1)
tiff(filename = "ERR in two scales for BD.tiff",width = 15,height = 25,units = "cm",res = 300,pointsize = 8)
opar <- par(mfrow = c(8,4),cex.axis=1,cex.lab=1.5,cex.main=1.2, mar=c(3.8,4,3.5,0.8)-0.5)
cen <- sapply(datalist, function(x) median(x$tm)) %>% mean()
x_all <- -23:32
x <- x_all
a <- rr_ts(basis = basvar_ts,coef = fit2$coefficients,vcov = fit2$vcov,at = x,cen = cen)
plot(x,a$RR, type="l", lwd=1, col="red",ylim = c(-0.5,5),mgp = c(2.2,1,0),
     xlab="Temperature",ylab="RR")
title("Average",line = 0.5)
redtrans <- rgb(255, 0, 0, 60, maxColorValue=255) 
polygon(c(rev(x),x),c(rev(a$RR_upper),a$RR_lower),col=redtrans, border = NA)
abline(h=1)


knotsi <- sapply(datalist, function(x) quantile(x$tm,knots_relative)) %>% apply(1,mean)
bound <- colMeans(ranges)
x <- seq(bound[1],bound[2],by = 1) %>% c(max(bound))
btmean <- onebasis(x,fun="bs",degree=degree,knots=knotsi,Bound=bound)
a <- crosspred(basis = btmean,coef = yfit1[[i]]$blup,vcov = yfit1[[i]]$vcov,at = x,cen = cen,model.link="log")
lines(x,a$allRRfit, lwd=1, col="blue")
redtrans <- rgb(0, 0,255, 60, maxColorValue=255) 
polygon(c(rev(x),x),c(rev(a$allRRhigh),a$allRRlow),col=redtrans, border = NA)
reaxis <- sapply(datalist, function(x) quantile(x$tm,c(0,25,50,75,100)/100)) %>% apply(1,mean)
axis(3,at = reaxis,labels = c(0,25,50,75,100),col = "blue",tcl = 0.25,mgp = c(2,-1.5,0),lwd = 0.8)
abline(v = cen,lty = 2,lwd = 0.8)
axis(1,at = cen,labels = round(cen,1),tcl = 0.25,mgp = c(2.2,-1.5,0))
abline(v = bound,lty = 2,col = "blue",lwd = 0.8)
for (i in seq(m)) {
  x <- x_all
  cen <- median(datalist[[i]]$tm)
  a <- rr_ts(basis = basvar_ts,coef = yfit2[[i]]$blup,vcov = yfit2[[i]]$vcov,at = x,cen = cen)
  plot(x,a$RR, type="l", lwd=1, col="red",ylim = c(-0.5,5),mgp = c(2.2,1,0),
       xlab="Temperature",ylab="RR")
  title(regions[i],line = 0.5)
  redtrans <- rgb(255, 0, 0, 60, maxColorValue=255) 
  polygon(c(rev(x),x),c(rev(a$RR_upper),a$RR_lower),col=redtrans, border = NA)
  abline(h=1)
  
  knotsi <- quantile(datalist[[i]]$tm,knots_relative)
  bound <- range(datalist[[i]]$tm)
  x <- seq(bound[1],bound[2],by = 1) %>% c(max(bound))
  btmean <- onebasis(x,fun="bs",degree=degree,knots=knotsi,Bound=bound)
  a <- crosspred(basis = btmean,coef = yfit1[[i]]$blup,vcov = yfit1[[i]]$vcov,at = x,cen = cen,model.link="log")
  lines(x,a$allRRfit,lwd=1, col="blue")
  redtrans <- rgb( 0, 0,255, 60, maxColorValue=255) 
  polygon(c(rev(x),x),c(rev(a$allRRhigh),a$allRRlow),col=redtrans, border = NA)
  reaxis <- quantile(datalist[[i]]$tm,c(0,25,50,75,100)/100)
  axis(3,at = reaxis,labels = c(0,25,50,75,100),col = "blue",tcl = 0.25,mgp = c(2.2,-1.5,0),lwd = 0.8)
  abline(v = cen,lty = 2,lwd = 0.8)
  axis(1,at = cen,labels = round(cen,1),tcl = 0.25,mgp = c(2,-1.5,0))
  abline(v = bound,lty = 2,col = "blue",lwd = 0.8)
}
par(opar)
dev.off()


