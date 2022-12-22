##############################################################################-#
# Using the data in Gasparrini et al.'s work to compare TS-based strategy and  #
# BS/NS-based strategy in studying ERRs in the original scale                  #
##############################################################################-#
library(dlnm) ; library(splines) ; library(xtable)
source("fun for TS.R")
################################################################################
############ load data and data process ########################################
################################################################################
data_all <- read.csv("regEngWales-from gasparrini.csv",row.names=1)
regions <- as.character(unique(data_all$regnames))
datalist <- lapply(regions, function(region) data_all[data_all$regnames==region,])
names(datalist) <- regions

m <- length(datalist)

# MOVING AVERAGE OF TMEAN OVER LAG 0-6
for(i in seq(datalist)) datalist[[i]]$tmean05 <- 
  filter(datalist[[i]]$tmean,rep(1,6)/6,side=1)

# TEMPERATURE RANGES (FOR LAG 0-5)
ranges <- t(sapply(datalist,function(x) range(x$tmean05,na.rm=T)))

# COMPUTE 25TH-75TH PERCENTILES OF META-VARIABLES

# DEFINE THE AVERAGE RANGE, CENTERING POINT, DEGREE AND TYPE OF THE SPLINE
# (THESE PARAMETERS CAN BE CHANGED BY THE USER FOR ADDITIONAL ANALYSES)
cen <- 17


# DEFINE THE KNOTS AT TEMPERATURE CORRESPONDING TO AVERAGE PERCENTILES
knotperc <- c(5,35,65,95)/100
a <- NULL
for (i in 1:m) {
  a <- c(a,as.numeric(datalist[[i]]$tmean05))
}
knots <- quantile(a,probs = knotperc,na.rm = T)
knots_ns <- quantile(a,probs = c(5,25,50,75,95)/100,na.rm = T)
bound <- range(ranges,na.rm = T)
## temperature range
tiff(filename = "temperature range for Gaspirirni's example.tiff",width = 10,height = 8,units = "cm",res = 300,pointsize = 8)
par(mar=c(5,6,4,1)+0.1)
plot(seq(-5,28,length=m),seq(m),type="n",yaxt="n",
     ylab="",xlab="Temperature")
axis(2,at=seq(m),labels=regions,las=1,cex.axis=0.7)
arrows(ranges[,1],seq(m),ranges[,2],seq(m),angle=90,length=0.05,code=3)
title(ylab="Province",mgp=c(4.5,1,0))
abline(v=knots,col = "red")
abline(v=bound,lty=2)
axis(3,at = c(bound,knots),labels = c("0%","100%","5%","35%","65%","95%"),tcl = -0.2,mgp = c(1,0.5,0))
dev.off()
################################################################################
########################### bs #################################################
################################################################################

degree <- 2
df <- length(knots) + degree
type <- "bs"

ymat <- matrix(NA,m,df,dimnames=list(regions,paste("spl",seq(df),sep="")))
Slist <- vector("list",m)
names(Slist) <- regions
for(i in seq(m)) {
  
  # PRINT ITERATION
  cat(i,"")
  
  # LOAD
  data <- datalist[[i]]
  
  # CREATE THE SPLINE
  # NB: KNOTS AND BOUNDARIES FIXED AT SAME VALUES
  btmean05 <- onebasis(data$tmean05,fun=type,degree=degree,knots=knots,
                       Bound=bound)
  
  # RUN THE MODEL
  model <- glm(death ~ btmean05 + dow + ns(time,7*14),family=quasipoisson(),data)
  
  # EXTRACT AND SAVE THE RELATED COEF AND VCOV
  predtmean05 <- crosspred(btmean05,model,cen=cen)
  ymat[i,] <- predtmean05$coef
  Slist[[i]] <- predtmean05$vcov
}
fit_bs <- mvmeta(ymat,Slist,method = "reml")
yfit_bs <- blup(fit_bs,vcov = T)
################################################################################
########################### ns #################################################
################################################################################
knots_ns <- knots
df <- length(knots_ns) + 1
type <- "ns"

ymat <- matrix(NA,m,df,dimnames=list(regions,paste("spl",seq(df),sep="")))
Slist <- vector("list",m)
names(Slist) <- regions
for(i in seq(m)) {
  
  # PRINT ITERATION
  cat(i,"")
  
  # LOAD
  data <- datalist[[i]]
  
  # CREATE THE SPLINE
  # NB: KNOTS AND BOUNDARIES FIXED AT SAME VALUES
  btmean05 <- onebasis(data$tmean05,fun=type,knots=knots_ns,
                       Bound=bound)
  
  # RUN THE MODEL
  model <- glm(death ~ btmean05 + dow + ns(time,7*14),family=quasipoisson(),data)
  
  # EXTRACT AND SAVE THE RELATED COEF AND VCOV
  predtmean05 <- crosspred(btmean05,model,cen=cen)
  ymat[i,] <- predtmean05$coef
  Slist[[i]] <- predtmean05$vcov
}
fit_ns <- mvmeta(ymat,Slist,method = "reml")
yfit_ns <- blup(fit_ns,vcov = T)


################################################################################
######################### ts ###################################################
################################################################################
degree <- 2
df <- length(knots) + degree
ymat <- matrix(NA,m,df,dimnames=list(rownames = regions))
Slist <- vector("list",m)
names(Slist) <- regions
for (i in seq(m)) {
  cat(i,"")
  data <- datalist[[i]]
  basvar <- trun_spline(x = data$tmean05,knots = knots,cenknot_lower = knots[2],cenknot_upper = knots[3],degree = degree)
  model <- glm(death ~ basvar + dow + ns(time,7*14),family=quasipoisson(),data)
  a <- grepl("basvarb[0-9]",names(model$coefficients))
  ymat[i,] <- model$coefficients[a]
  Slist[[i]] <- vcov(model)[a,a]
}
fit_ts <- mvmeta(ymat,Slist,method = "reml")
yfit_ts <- blup(fit_ts,vcov = T)
############## comparison between bs and ts -------
tiff(filename = "bs vs ts for Gaspirini's example.tiff",width = 12,height = 10,units = "cm",res = 300,pointsize = 8)
opar <- par(mfrow = c(3,4),cex.axis=1,cex.lab=1.5,cex.main=1.2, mar=c(3.8,4,3.5,0.8)-0.5)
x <- seq(range(ranges)[1],range(ranges)[2],by = 1) %>% c(max(ranges))
a <- rr_ts(basis = basvar,coef = fit_ts$coefficients,vcov = fit_ts$vcov,at = x,cen = cen)
plot(x,a$RR, type="l", lwd=1, col="red",ylim = c(0.5,2),# ylim=range(c(a$RR_upper,a$RR_lower)),
     xlab="Temperature",ylab="RR",mgp = c(2.2,1,0))
title("Average",line = 0.5)
redtrans <- rgb(255, 0, 0, 60, maxColorValue=255) 
polygon(c(rev(x),x),c(rev(a$RR_upper),a$RR_lower),col=redtrans, border = NA)
abline(h=1)
abline(v = range(ranges),lty = 2)
btmean <- onebasis(x,fun="bs",degree=degree,knots=knots,Bound=bound)
a <- crosspred(basis = btmean,coef = coef(fit_bs),vcov = vcov(fit_bs),at = x,cen = cen,model.link="log")
lines(x,a$allRRfit,col = "blue",lty = 5)
redtrans <- rgb( 0, 0,255, 60, maxColorValue=255) 
polygon(c(rev(x),x),c(rev(a$allRRhigh),a$allRRlow),col=redtrans, border = NA,density = 40)
legend("top",legend = c("ts","bs"),lty = c(1,5),col = c("red","blue"),bty = "n",ncol = 1)

for (i in seq(m)) {
  a <- rr_ts(basis = basvar,coef = yfit_ts[[i]]$blup,vcov = yfit_ts[[i]]$vcov,at = x,cen = cen)
  plot(x,a$RR, type="l", lwd=1, col="red",ylim = c(0.5,2),# ylim=range(c(a$RR_upper,a$RR_lower)),
       xlab="Temperature",ylab="RR",mgp = c(2.2,1,0))
  title(regions[i],line = 0.5)
  redtrans <- rgb(255, 0, 0, 60, maxColorValue=255) 
  polygon(c(rev(x),x),c(rev(a$RR_upper),a$RR_lower),col=redtrans, border = NA)
  abline(h=1)
  abline(v = c(ranges[i,]),lty = 2)
  a <- crosspred(basis = btmean,coef = yfit_bs[[i]]$blup,vcov = yfit_bs[[i]]$vcov,at = x,cen = cen,model.link="log")
  lines(x,a$allRRfit,col = "blue",lty = 5)
  redtrans <- rgb( 0, 0,255, 60, maxColorValue=255) 
  polygon(c(rev(x),x),c(rev(a$allRRhigh),a$allRRlow),col=redtrans, border = NA,density = 40)
  legend("top",legend = c("ts","bs"),lty = c(1,5),col = c("red","blue"),bty = "n",ncol = 1)
}
par(opar)
dev.off()

############## comparison between ns and ts-------
tiff(filename = "ns vs ts for Gaspirini's example.tiff",width = 12,height = 10,units = "cm",res = 300,pointsize = 8)
opar <- par(mfrow = c(3,4),cex.axis=1,cex.lab=1.5,cex.main=1.2, mar=c(3.8,4,3.5,0.8)-0.5)
x <- seq(range(ranges)[1],range(ranges)[2],by = 1) %>% c(max(ranges))
a <- rr_ts(basis = basvar,coef = fit_ts$coefficients,vcov = fit_ts$vcov,at = x,cen = cen)
plot(x,a$RR, type="l", lwd=1, col="red",ylim = c(0.5,2),# ylim=range(c(a$RR_upper,a$RR_lower)),
     xlab="Temperature",ylab="RR",mgp = c(2.2,1,0))
title("Average",line = 0.5)
redtrans <- rgb(255, 0, 0, 60, maxColorValue=255) 
polygon(c(rev(x),x),c(rev(a$RR_upper),a$RR_lower),col=redtrans, border = NA)
abline(h=1)
abline(v = range(ranges),lty = 2)
btmean <- onebasis(x,fun="ns",knots=knots_ns,Bound=bound)
a <- crosspred(basis = btmean,coef = coef(fit_ns),vcov = vcov(fit_ns),at = x,cen = cen,model.link="log")
lines(x,a$allRRfit,col = "blue",lty = 5)
redtrans <- rgb( 0, 0,255, 60, maxColorValue=255) 
polygon(c(rev(x),x),c(rev(a$allRRhigh),a$allRRlow),col=redtrans, border = NA)
legend("top",legend = c("ts","ns"),lty = c(1,5),col = c("red","blue"),bty = "n",ncol = 1)

for (i in seq(m)) {
  a <- rr_ts(basis = basvar,coef = yfit_ts[[i]]$blup,vcov = yfit_ts[[i]]$vcov,at = x,cen = cen)
  plot(x,a$RR, type="l", lwd=1, col="red",ylim = c(0.5,2),# ylim=range(c(a$RR_upper,a$RR_lower)),
       xlab="Temperature",ylab="RR",mgp = c(2.2,1,0))
  title(regions[i],line = 0.5)
  redtrans <- rgb(255, 0, 0, 60, maxColorValue=255) 
  polygon(c(rev(x),x),c(rev(a$RR_upper),a$RR_lower),col=redtrans, border = NA)
  abline(h=1)
  abline(v = c(ranges[i,]),lty = 2)
  a <- crosspred(basis = btmean,coef = yfit_ns[[i]]$blup,vcov = yfit_ns[[i]]$vcov,at = x,cen = cen,model.link="log")
  lines(x,a$allRRfit,col = "blue",lty = 5)
  redtrans <- rgb( 0, 0,255, 60, maxColorValue=255) 
  polygon(c(rev(x),x),c(rev(a$allRRhigh),a$allRRlow),col=redtrans, border = NA)
  legend("top",legend = c("ts","ns"),lty = c(1,5),col = c("red","blue"),bty = "n",ncol = 1)
}
par(opar)
dev.off()



