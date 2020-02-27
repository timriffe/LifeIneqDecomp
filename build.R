
library(devtools)
document()
check()


library(HMDHFDplus)
library(reshape2)
load_all()
LT<-readHMDweb("USA","mltper_1x1",us,pw)
age <- 0:110
yrs <- c(1950,1975, 2000)


ax <- acast(LT[LT$Year %in% yrs, ], Age~Year,value.var = "ax")
dx <- acast(LT[LT$Year %in% yrs, ], Age~Year,value.var = "dx")
lx <- acast(LT[LT$Year %in% yrs, ], Age~Year,value.var = "lx")
ex <- acast(LT[LT$Year %in% yrs, ], Age~Year,value.var = "ex")

p<- c(4,7,20)
args(bw_decomp)
bw_decomp(age=age,ax=ax,dx=dx,lx=lx,ex=ex,p=p, method = "variance")
bw_decomp(age=age,ax=ax,dx=dx,lx=lx,ex=ex,p=p, method = "Gini")
bw_decomp(age=age,ax=ax,dx=dx,lx=lx,ex=ex,p=p, method = "Theil")
bw_decomp(age=age,ax=ax,dx=dx,lx=lx,ex=ex,p=p, method = "MLD")
bw_decomp(age=age,ax=ax,dx=dx,lx=lx,ex=ex,p=p, method = "edag")
