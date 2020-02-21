
# library(devtools)
# install_github("timriffe/LifeIneq")
# 
# library(HMDHFDplus)
# LT<-readHMDweb("USA","mltper_1x1",us,pw)
# library(LifeIneq)
# age <- 0:110
# yrs <- c(1950,1975, 2000)
# 
# library(reshape2)
# ax <- acast(LT[LT$Year %in% yrs, ], Age~Year,value.var = "ax")
# dx <- acast(LT[LT$Year %in% yrs, ], Age~Year,value.var = "dx")
# lx <- acast(LT[LT$Year %in% yrs, ], Age~Year,value.var = "lx")
# ex <- acast(LT[LT$Year %in% yrs, ], Age~Year,value.var = "ex")
# 
# p<- c(4,7,20)



#' @title between-within decomposition of lifespan inequality measures
#' @description Partition a lifespan inequality index into additive components of between-group inequality and within-group inequality. Presently implemented for Theil's index, e-edagger, variance, mean log deviation, and the gini coeficient. 
#' 
#' @param age numeric. vector of lower age bounds.
#' @param dx numeric. vector of the lifetable death distribution.
#' @param lx numeric. vector of the lifetable survivorship.
#' @param ex numeric. vector of remaining life expectancy.
#' @param ax numeric. vector of the average time spent in the age interval of those dying within the interval.
#' @import LifeIneq

bw_decomp <- function(age, ax, dx, lx, ex, p, 
                      measure = c("Theil", "edag","variance","MLD")){
  
  # check dims
  K <- length(p)
  stopifnot(all.equal(ncol(ax),
            ncol(dx), 
            ncol(lx), 
            ncol(ex), 
            K))
  N <- length(age)
  stopifnot(all.equal(nrow(ax),
                      nrow(dx), 
                      nrow(lx), 
                      nrow(ex),
                      N))
  
  # validate measure selection
  measure <- match.arg(measure)
  
  # 1) standardize inputs
  p    <- p / sum(p)
  lx   <- lx %*% diag(1 / lx[1, ])
  dx   <- dx %*% diag(1 / colSums(dx))
  
  # 2) get pop avgs and other goods
  
  # weighted dx and lx
  pdxm <- dx %*% diag(p)
  plxm <- lx %*% diag(p)
   
  # pop avg dx and lx
  plx  <- rowSums(plxm)
  pdx  <- rowSums(pdxm)
   
  # need rows to sum to one to weight ax and ex
  plxc <- diag(1 / plx) %*% plxm
  pdxc <- diag(1 / pdx) %*% pdxm
   
  # pop avg ax is dx-weighted ax
  pax  <- rowSums(ax * pdxc)
  # pop avg ex is lx-weighted ex
  pex  <- rowSums(ex * plxc)
   
  # ---------------------------------------------
  # THIS WILL BE SIMPLIFIED WHEN WE HAVE A WRAPPER
  # call up our inequality function
  f    <-  match.fun(paste0("ineq_",measure))
  
  args.need <- names(formals(f))
   
  
   
  # only Gini doesn't use dx
  usedx <- ifelse(measure == "Gini",FALSE, TRUE)
   
  # calculate inequality index
  #for each of the k subgroups and the total
  indices <- rep(0, K)
  if (usedx){
   tot <- f(age = age, 
            dx = pdx, 
            lx = plx, 
            ax = pax, 
            ex = pex)[1]
     for (k in 1:K){
       indices[k] <- f(age = age, 
                       dx = dx[, k], 
                       lx = lx[, k], 
                       ax = ax[, k], 
                       ex = ex[, k])[1]
     }
     
   } else{
     tot <- f(age = age, 
              lx = plx, 
              ax = pax, 
              ex = pex)[1]
     for (k in 1:K){
       indices[k] <- f(age = age, 
                       lx = lx[, k], 
                       ax = ax[, k], 
                       ex = ex[, k])[1]
     }
   }
   # within weighting depends on the measure
   if (measure %in% c("edag","variance","MLD")){ 
     weights <- p
   }
   
  # END part to be simplified w wrapper
  # ---------------------------------------------
  
   # Gini weights same as Theil weights when
   # radices all equal.
   # if (measure %in% "Gini"){ 
   #   w1 <- (lx[1,]^2 * ex[1,])
   #   w2 <- (plx[1]^2 * pex[1])
   #   weights <- (w1 / w2) * p
   # }
   
   if (measure %in% c("Theil","Gini")){
     weights <- p* ex[1, ] / pex[1]
   }
   
   # Combine
   W <- sum(weights * indices)
   
   # Between
   B <- tot - W
   
   # gather minimal goods to return
   out <- list(measure = measure,
               group_ind = indices,
               tot = tot,
               B = B,
               W = W,
               fB = B / tot,
               fW = W / tot)
   return(out)
}

