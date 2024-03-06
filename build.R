
library(devtools)
document()
check()


library(HMDHFDplus)
library(reshape2)
load_all()
LT<-readHMDweb("USA","mltper_1x1",Sys.getenv("us"),Sys.getenv("pw"))
LT2<-readHMDweb("JPN","mltper_1x1",Sys.getenv("us"),Sys.getenv("pw"))
age <- 0:110
yrs <- c(1950,1975, 2000)


ax <- acast(LT[LT$Year %in% yrs, ], Age~Year,value.var = "ax")
dx <- acast(LT[LT$Year %in% yrs, ], Age~Year,value.var = "dx")
lx <- acast(LT[LT$Year %in% yrs, ], Age~Year,value.var = "lx")
ex <- acast(LT[LT$Year %in% yrs, ], Age~Year,value.var = "ex")
mx1 <- acast(LT[LT$Year %in% yrs, ], Age~Year,value.var = "mx")
mx2 <- acast(LT2[LT2$Year %in% yrs, ], Age~Year,value.var = "mx")
p<- c(4,7,20)
args(bw_decomp)
bw_decomp(age,ax,dx,lx,ex,distribution_type="rl",prop=p, "var")
args(bw_decomp)

bw_decomp(age=age,ax=ax,dx=dx,lx=lx,ex=ex,p=p, method = "gini")
bw_decomp(age=age,ax=ax,dx=dx,lx=lx,ex=ex,p=p, method = "theil")
bw_decomp(age=age,ax=ax,dx=dx,lx=lx,ex=ex,p=p, method = "mld")
bw_decomp(age=age,ax=ax,dx=dx,lx=lx,ex=ex,p=p, method = "edag")


# what about H, cov, sd, aid, mad
# 
# Random curiosity test: age decomp of a fraction-between component. Because yeah
# my_dec_fun <- function(pars, n, age, p, method = "edag", component="fB"){
#   dim(pars) <- c(length(pars)/n,n)
#   LTs <- apply(pars, 2, DemoTools::lt_single_mx, OAG = TRUE)
#   ax  <- lapply(LTs,pull,"nAx") %>% do.call("cbind",.)
#   dx  <- lapply(LTs,pull,"ndx") %>% do.call("cbind",.)
#   lx  <- lapply(LTs,pull,"lx") %>% do.call("cbind",.)
#   ex  <- lapply(LTs,pull,"ex") %>% do.call("cbind",.)
# 
#   bw_decomp(age=age,ax=ax,dx=dx,lx=lx,ex=ex,p=p, method = method)[[component]]
# }
# 
# library(DemoDecomp)
# 
# cc <- horiuchi(my_dec_fun,
#          pars1=c(mx1),
#          pars2=c(mx2),
#          age=0:110,
#          n=3,
#          p=c(1,5,10),
#          method="edag",
#          component="fB",
#          N=10)
# args(horiuchi)
# 
# dim(cc) <- c(111,3)
# 
# matplot(0:110,cc,type='l')
# plot(rowSums(cc))
# sum(cc)
