########### The functions for TS-based GAM =====================
`%+%` <- function(x,y) paste0(x,y)
`%>%` <- magrittr::`%>%`
trun_spline <- function(x,knots,cenknot_lower = knots[1],cenknot_upper = knots[2],degree = 2){
  a <- sort(knots)
  if(any(a != knots)) stop(" Knots must be increasingly ordered")
  if(!cenknot_lower %in% knots) stop("cenknot_lower must be one of knots")
  if(!cenknot_upper %in% knots) stop("cenknot_upper must be one of knots")
  if(match(cenknot_upper,knots) - match(cenknot_lower,knots) != 1)
    stop("cenknot_upper and cenknot_lower must be two adjacent knots")
  
  nx <- names(x)
  x <- as.vector(x)
  if(degree == 2){
    basis <- cbind(x,x^2)
    for (knot in knots) {
      if(knot <= cenknot_lower){
        a <- (x-knot)^2*(x<=knot)
      } else if(knot >= cenknot_upper) a <- (x-knot)^2*(x>=knot)
      else next
      basis <- cbind(basis,a)
    } 
  }
  else if(degree == 3){
    basis <- cbind(x,x^2,x^3)
    for (knot in knots) {
      if(knot <= cenknot_lower){
        a <- (x-knot)^3*(x<=knot)
      } else if(knot >= cenknot_upper) a <- (x-knot)^3*(x>=knot)
      else next
      basis <- cbind(basis,a)
    } 
  } else stop("degree must be 2 or 3")
  dimnames(basis) <- list(nx, paste("b", seq(ncol(basis)),sep = ""))
  attr(basis,"knots") <- knots
  attr(basis,"cenknot_upper") <- cenknot_upper
  attr(basis,"cenknot_lower") <- cenknot_lower
  attr(basis,"degree") <- degree
  class(basis) <- c("onebasis", "matrix")
  basis
}

rr_ts <- function(basis,coef,vcov,at = NULL, cen = NULL, ci.level = 0.95){
  coef <- as.vector(coef)
  arg <- attributes(basis)[c("degree","cenknot_upper","cenknot_lower","knots")]
  basisat <- do.call("trun_spline",modifyList(arg,list(x = at)))
  basiscen <- do.call("trun_spline",modifyList(arg,list(x = cen)))
  matrixat <- t(t(basisat) - as.vector(basiscen))
  yat <- (matrixat %*% coef) %>% as.vector()
  varat <- diag(matrixat %*% vcov %*% t(matrixat))  
  names(yat) <- names(varat) <- at
  logRR_se = sqrt(varat)
  RR <- exp(yat)
  RR_upper <- exp(yat + qnorm(1-(1-ci.level)/2)*logRR_se)
  RR_lower <- exp(yat - qnorm(1-(1-ci.level)/2)*logRR_se)
  
  res <- list(RR = exp(yat),RR_upper = RR_upper,RR_lower = RR_lower, 
              logRR = yat, logRR_se = logRR_se,cen = cen)
  res
}

fqbic <- function(model) {
  loglik <- sum(dpois(model$y,model$fitted.values,log=TRUE))
  phi <- summary(model)$dispersion
  qaic <- -2*loglik + log(length(model$fitted.values))*summary(model)$df[3]*phi
  return(qaic)
}