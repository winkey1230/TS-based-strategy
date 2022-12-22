##############################################################################-#
# Illustrate the basis functions for bs, ns and ts 
# Illustrate the complete multicolinearity 
##############################################################################-#
library(splines)
source("fun for TS.R")
x <- -10:110/100
knots <- c(0.2,0.4,0.6,0.8)# knots <- seq(20,80,by = 10)/100
bound <- c(0,1)

################################################################################
############## plot basis functions for ns and bs               ################
############# illustrate the problem related to singular matrix ################
################################################################################
tiff(filename = "basis functions for ns and bs.tiff",width = 15,height = 7,units = "cm",res = 300,pointsize = 6)
opar <- par(mfrow = c(1,2),mar=c(4,4,2,1))
########## ns ---------
nsvar <- ns(x,knots = knots,intercept = F,Boundary.knots = bound)
plot(range(x),range(nsvar),type = "n",main = "Basis functions for natural splines",xlab = "x",
     ylab = expression(italic(b[j](x))),font.lab = 3,mgp = c(2,1,0),cex.lab=1.5)
for (i in seq(ncol(nsvar))) {
  lines(x,nsvar[,i],col = i)
}
points(knots,rep(0,length(knots)),pch = 16)
points(bound,c(0,0),pch = 16,col = "red")
legend_lab <- vector(mode = "expression",length = ncol(nsvar))
for (i in 1:ncol(nsvar)) {
  legend_lab[i] <- substitute(italic(b[i]),list(i = i)) %>% as.expression()
}
legend("topleft",legend = legend_lab,lty = 1,col = seq(ncol(nsvar)),cex = 1.2)
legend("top",legend = c("Knots","Bounds"),pch = 16,col = c("black","red"),cex = 1.2,ncol = 2,bty = "n")

## adjudge singularity within a local range
index_row <- x >= 0.2 & x <= 0.6
index_col <- apply(nsvar[index_row,], 2, function(x) length(unique(x)) != 1)
xx <- cbind(1,nsvar[index_row,index_col])
qr(xx)$rank == ncol(xx) # false: complete multicollinearity

############# bs ---------
bsvar <- bs(x,knots = knots,intercept = F,degree = 2,Boundary.knots = bound)
plot(range(x),range(bsvar),type = "n",main = "Basis functions for B-splines",xlab = "x",
     ylab = expression(italic(b[j](x))),font.lab = 3,mgp = c(2,1,0),cex.lab=1.5)
for (i in seq(ncol(bsvar))) {
  lines(x,bsvar[,i],col = i)
}
points(knots,rep(0,length(knots)),pch = 16)
points(bound,c(0,0),pch = 16,col = "red")
legend_lab <- vector(mode = "expression",length = ncol(bsvar))
for (i in 1:ncol(bsvar)) {
  legend_lab[i] <- substitute(italic(b[i]),list(i = i)) %>% as.expression()
}
legend("topleft",legend = legend_lab,lty = 1,col = seq(ncol(bsvar)),cex = 1.2)
legend("top",legend = c("Knots","Bounds"),pch = 16,col = c("black","red"),cex = 1.2,ncol = 2,bty = "n")

## adjudge singularity within a local range
index_row <- x >= 0.2 & x <= 0.6
index_col <- apply(bsvar[index_row,], 2, function(x) length(unique(x)) != 1)
xx <- cbind(1,bsvar[index_row,index_col])
qr(xx)$rank == ncol(xx) # false: complete multicollinearity
# ------
par(opar)
dev.off()
################################################################################
############## plot basis functions for ts                      ################
############# illustrate the problem related to singular matrix ################
################################################################################
tiff(filename = "basis functions for TS.tiff",width = 15/2,height = 7,units = "cm",res = 300,pointsize = 6)
########## ts -----------------
tsvar <- trun_spline(x,knots = knots,cenknot_lower = knots[2],cenknot_upper = knots[3],degree = 2)
plot(range(x),range(tsvar),type = "n",main = "Basis functions for TS",
     xlab = "x",ylab = expression(italic(b[j](x))),font.lab = 3,mgp = c(2,1,0),cex.lab=1.5)
for (i in seq(ncol(tsvar))) {
  lines(x,tsvar[,i],col = i)
}
points(knots,rep(0,length(knots)),pch = 16)
legend_lab <- vector(mode = "expression",length = ncol(tsvar))
for (i in 1:ncol(tsvar)) {
  legend_lab[i] <- substitute(italic(b[i]),list(i = i)) %>% as.expression()
}
legend("topleft",legend = legend_lab,lty = 1,col = seq(ncol(tsvar)),cex = 1.2)
legend("top",legend = "Knots",pch = 16,cex = 1.2,bty = "n")

## adjudge singularity within a local range
index_row <- x >= 0.2 & x <= 0.6
index_col <- apply(tsvar[index_row,], 2, function(x) length(unique(x)) != 1)
xx <- cbind(1,tsvar[index_row,index_col])
qr(xx)$rank == ncol(xx) # false: complete multicollinearity
# -----
dev.off()


######### magnifying ts==========================
tiff(filename = "TS-magnify.tiff",width = 15/2,height = 7,units = "cm",res = 300,pointsize = 6)
tsvar <- trun_spline(x,knots = knots,cenknot_lower = knots[2],cenknot_upper = knots[3],degree = 2)
plot(range(x),range(tsvar),type = "n",ylim = c(0,0.3),
     xlab = "x",ylab = expression(italic(b[j](x))),font.lab = 3,mgp = c(2,1,0),cex.lab=1.5)
for (i in seq(ncol(tsvar))) {
  lines(x,tsvar[,i],col = i)
}
points(knots,rep(0,length(knots)),pch = 16)
legend_lab <- vector(mode = "expression",length = ncol(tsvar))
for (i in 1:ncol(tsvar)) {
  legend_lab[i] <- substitute(italic(b[i]),list(i = i)) %>% as.expression()
}
dev.off()

