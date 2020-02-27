
#' @title between-within decomposition of lifespan inequality measures
#' @description Partition a lifespan inequality index into additive components of between-group inequality and within-group inequality. Presently implemented for Theil's index, e-edagger, variance, mean log deviation, and the gini coeficient. 
#' 
#' @param age numeric. vector of lower age bounds.
#' @param dx numeric. vector of the lifetable death distribution.
#' @param lx numeric. vector of the lifetable survivorship.
#' @param ex numeric. vector of remaining life expectancy.
#' @param ax numeric. vector of the average time spent in the age interval of those dying within the interval.
#' @param method character one of \code{"Theil", "edag","variance","MLD","Gini"}
#' @import LifeIneq

bw_decomp <- function(age, ax, dx, lx, ex, prop,
                      method = c("Theil", "edag","variance","MLD","Gini")){
  
  # check dims
  K <- length(prop)
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
  
  # validate method selection
  method <- match.arg(method)
  
  # 1) standardize inputs
  prop <- prop / sum(prop)
  lx   <- lx %*% diag(1 / lx[1, ])
  dx   <- dx %*% diag(1 / colSums(dx))
  
  # 2) get pop avgs and other goods
  
  # weighted dx and lx
  pdxm <- dx %*% diag(prop)
  plxm <- lx %*% diag(prop)
   
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
   

 

  # calculate inequality index
  
  tot <- ineq(age = age, 
              dx = pdx, 
              lx = plx, 
              ax = pax, 
              ex = pex,
              method = method)[1]
  
  # again for each of the k subgroups and the total
  indices <- rep(0, K)
  for (k in 1:K){
    indices[k] <- ineq(age = age, 
                       dx = dx[, k], 
                       lx = lx[, k], 
                       ax = ax[, k], 
                       ex = ex[, k],
                       method = method)[1]
   }

   # within weighting depends on the measure
   if (method %in% c("edag","variance","MLD")){ 
     weights <- prop
   }
   
   if (method %in% c("Theil","Gini")){
     weights <- prop * ex[1, ] / pex[1]
   }
   
   # Combine
   W <- sum(weights * indices)
   
   # Between
   B <- tot - W
   
   # gather minimal goods to return
   out <- list(method = method,
               group_ind = indices,
               tot = tot,
               B = B,
               W = W,
               fB = B / tot,
               fW = W / tot)
   return(out)
}

