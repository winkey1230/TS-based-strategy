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
degree <- 2
df <- length(knots) + degree
cen <- median(data_all$tm)


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
############ ERRs in the original scale using BS ###############################
################################################################################
ymat <- matrix(NA,m,df,dimnames=list(rownames = regions))
Slist <- vector("list",m)
names(Slist) <- regions

for (i in seq(m)) {
  cat(regions[i]," ")
  data <- datalist[[i]]
  data$time <- 1:dim(data)[1]
  data$lograin <- log(data$rain)
  basvar <- onebasis(data$tm,fun = "bs",knots = knots,degree = degree,Bound=bound)
  model <- glm(count ~ basvar + humid + sun + ns(time,df = 6),family=quasipoisson(),data)
  a <- grepl("basvarb[0-9]",names(model$coefficients))
  
  a1 <- model$coefficients[a]
  ymat[i,] <- a1
  
  a2 <- vcov(model)[a,a]
  Slist[[i]] <- a2
}
fit1 <- mvmeta(ymat,Slist,method = "reml",control = list(maxiter = 1000))
yfit1 <- blup(fit1,vcov = T)
basvar_bs <- basvar

################################################################################
############### ERRs in the original scale using TS ############################
################################################################################
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
basvar_ts <- basvar


################################################################################
############### plot ERRs in the original scale using TS and bs ################
################################################################################
options(warn=-1)
tiff(filename = "Province-specific ERRs for BD in the original scale for ts vs bs.tiff",width = 15,height = 25,units = "cm",res = 300,pointsize = 8)
opar <- par(mfrow = c(8,4),cex.axis=1,cex.lab=1,cex.main=1.2, mar=c(3.8,4,3.5,0.8)-0.5)
x <- seq(range(ranges)[1],range(ranges)[2],by = 1) %>% c(max(ranges))
a <- rr_ts(basis = basvar_ts,coef = fit2$coefficients,vcov = fit2$vcov,at = x,cen = cen)
plot(x,a$RR, type="l", lwd=1, col="red",ylim = c(0.2,6),
     xlab="Temperature",ylab="RR",cex.lab = 1.5,mgp = c(2.2,1,0))
title("Average",line = 0.5)
redtrans <- rgb(255, 0, 0, 60, maxColorValue=255) 
polygon(c(rev(x),x),c(rev(a$RR_upper),a$RR_lower),col=redtrans, border = NA)

# basvar_bs <- onebasis(x,fun="bs",degree=degree,knots=knots,Bound=bound)
a <- crosspred(basis = basvar_bs,coef = fit1$coefficients[1,],vcov = fit1$vcov,at = x,cen = cen,model.link="log")
lines(x,a$allRRfit, lwd=1, col="blue")
redtrans <- rgb(0, 0,255, 60, maxColorValue=255) 
polygon(c(rev(x),x),c(rev(a$allRRhigh),a$allRRlow),col=redtrans, border = NA)
abline(h=1)
abline(v = range(ranges),lty = 2,lwd = 0.8)
legend("top",legend = c("ts","bs"),lty = c(1,1),col = c("red","blue"),bty = "n",ncol = 1)
for (i in seq(m)) {
  a <- rr_ts(basis = basvar_ts,coef = yfit2[[i]]$blup,vcov = yfit2[[i]]$vcov,at = x,cen = cen)
  plot(x,a$RR, type="l", lwd=1, col="red",ylim = c(0.2,6),
       xlab="Temperature",ylab="RR",cex.lab = 1.5,mgp = c(2.2,1,0))
  title(regions[i],line = 0.5)
  redtrans <- rgb(255, 0, 0, 60, maxColorValue=255) 
  polygon(c(rev(x),x),c(rev(a$RR_upper),a$RR_lower),col=redtrans, border = NA)
  
  a <- crosspred(basis = basvar_bs,coef = yfit1[[i]]$blup,vcov = yfit1[[i]]$vcov,at = x,cen = cen,model.link="log")
  lines(x,a$allRRfit,lwd=1, col="blue")
  redtrans <- rgb( 0, 0,255, 60, maxColorValue=255) 
  polygon(c(rev(x),x),c(rev(a$allRRhigh),a$allRRlow),col=redtrans, border = NA)
  abline(h=1)
  abline(v = c(ranges[i,]),lty = 2,lwd = 0.8)
  legend("top",legend = c("ts","bs"),lty = c(1,1),col = c("red","blue"),bty = "n",ncol = 1)
  
}
par(opar)
dev.off()

