######
## Generate data for precision matrix estimation
######

rm(list=ls())

if (!require(mvtnorm)) install.packages('mvtnorm')
if (!require(moments)) install.packages('moments')
if (!require(Matrix)) install.packages('Matrix')
library(mvtnorm)
library(moments)
library(Matrix)

#Set number of dimension and number of observations
p=10; n=25
#Generate m data sets
m = 5
#########################################################
#Case One: hubs structure
# sigma_inv = diag(1,nrow=p,ncol=p)
# alpha = 0.25      #for magnitude of nonzero partial covariance
# num_hub = p/10    #number of hubs
# for (i in 1:num_hub) {
#   sigma_inv[(1+p/num_hub*(i-1)),(2+p/num_hub*(i-1)):(p/num_hub+p/num_hub*(i-1))] = alpha
#   sigma_inv[(2+p/num_hub*(i-1)):(p/num_hub+p/num_hub*(i-1)),(1+p/num_hub*(i-1))] = alpha
# }
# eigen(sigma_inv)$values
# sigma = solve(sigma_inv); sigma = as.matrix(sigma)
# 
# write.table(sigma_inv,file=paste0("./Data/GHS_sim_p",p,"hubs_sigmainv.csv"),sep=",",row.names=FALSE,col.names=FALSE)
############################################
#Case Two: random structure
set.seed(2016)
#generate off-diagonal elements in \Omega
a = rep(NA,p*(p-1)/2)
for (i in 1:(p*(p-1)/2)) {
  if (runif(1)>0.3) {a[i] = 0}    #elements nonzero with prob=0.01, for p=10
  #if (runif(1)>0.01) {a[i] = 0}    #elements nonzero with prob=0.01, for p=100
  # if (runif(1)>0.002) {a[i] = 0}    #elements nonzero with prob=0.002, for p=200
  else {a[i] = -runif(1,0.2,1)}
}
A = matrix(NA,nrow=p,ncol=p)
A[upper.tri(A)] = a; diag(A) = 1
sigma_inv = forceSymmetric(A); sigma_inv = as.matrix(sigma_inv)
#if \Omega is not positive-definite, generate again till it is
while (min(eigen(sigma_inv)$values) < 0.01) {
  for (i in 1:(p*(p-1)/2)) {
    
    if (runif(1)>0.3) {a[i] = 0}    #elements nonzero with prob=0.01, for p=10
    #if (runif(1)>0.01) {a[i] = 0}    #elements nonzero with prob=0.01, for p=100
    # if (runif(1)>0.002) {a[i] = 0}    #elements nonzero with prob=0.002, for p=200
    else {a[i] = -runif(1,0.2,1)}
  }
  A = matrix(NA,nrow=p,ncol=p)
  A[upper.tri(A)] = a; diag(A) = 1
  sigma_inv = forceSymmetric(A); sigma_inv = as.matrix(sigma_inv)
}
eigen(sigma_inv)$values
sigma = solve(sigma_inv); sigma = as.matrix(sigma)

write.table(sigma_inv,file=paste0("./Data/GHS_sim_p",p,"random_sigmainv.csv"),sep=",",row.names=FALSE,col.names=FALSE)
############################################
omega_elements = t(sigma_inv)[lower.tri(sigma_inv,diag=FALSE)]
sum(omega_elements!=0)

xx_list = list(NA)
for (i in 1:m) {
  set.seed(2050+i)
  xx_list[[i]] = rmvnorm(n=n,mean=rep(0,p),sigma=sigma)
  
  #write.table(xx_list[[i]],file=paste0("./Data/GHS_sim_p",p,"hubs",n,"_data", i, ".csv"),sep=",",row.names=FALSE,col.names=FALSE)
  write.table(xx_list[[i]],file=paste0("./Data/GHS_sim_p",p,"random",n,"_data", i, ".csv"),sep=",",row.names=FALSE,col.names=FALSE)
  }

save(xx_list,sigma_inv,sigma,omega_elements,n,p,
	#file=paste0("./Data/GHS_sim_p",p,"hubs",n,"_data.RData")
	file=paste0("./Data/GHS_sim_p",p,"random",n,"_data.RData")
	)
