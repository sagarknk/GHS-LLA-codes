######
## Estimate precision matrix by frequentist GLS
######

n = 25
p = 10
m = 5    #number of replicates

prec_struc = 2

if (prec_struc == 1){
  struc_label = "hubs"
}else{
  struc_label = "random"
}

load(file=paste0("./Data/GHS_sim_p",p,struc_label,n,"_data.RData"))
	
if (!require(glasso)) install.packages('glasso')
if (!require(cvTools)) install.packages('cvTools')
library(glasso)
library(cvTools)

###########################################
#Cross validation for Graphical LASSO
log_likelihood <- function(precision, emp_cov) {
  p      <- nrow(precision)
  logdet <- determinant(precision, logarithm=T)$modulus
  loglik <- 0.5 * (logdet - sum(emp_cov * precision) - p*log(2*pi))
  return(as.numeric(loglik))
}

glasso_cv <- function(ts, k=5, rholist, verbose=T, penalize.diag) {
  library(glasso)
  library(cvTools)
  n     <- nrow(ts)
  folds <- cvFolds(n, k, type="consecutive")
  
  loglikes <- sapply(1:k, function(ki) {
#     if (verbose) cat("Fold ", ki, "\n")
    S_train <- cov(ts[folds$which!=ki,])
    S_test  <- cov(ts[folds$which==ki,])
    if (penalize.diag==T) {
      GLP     <- glassopath(S_train, rholist=rho_seq, trace=0, penalize.diagonal=TRUE, maxit=20) }
    else if (penalize.diag==F) {
      GLP     <- glassopath(S_train, rholist=rho_seq, trace=0, penalize.diagonal=FALSE, maxit=20) }
    loglike <- apply(GLP$wi, 3, function(P_train) log_likelihood(P_train, S_test))
    loglike
  })
  
  ind     <- which.max(rowMeans(loglikes))
  rhomax  <- rholist[ind]
  S       <- cov(ts)
  a       <- glasso(S, rhomax)
  a$rhomax <- rhomax
    
  return(a)
}


# sequence of rho (tuning parameter) for CV:
rho_seq = exp(seq(from=-3,to=0,by=0.1))
GL_est_T_list = list(NA); GL_est_F_list = list(NA)
GL_cov_T_list = list(NA); GL_cov_F_list = list(NA)
rhomax_GLT = rep(NA,m); rhomax_GLF = rep(NA,m)
GLT_time = rep(NA,m); GLF_time = rep(NA,m)
for (i in 1:m) {
  print(i)
  xx = xx_list[[i]]
  
  t1 = Sys.time()
  gl_cv_T = glasso_cv(ts=xx, k=5, verbose=T, rholist=rho_seq, penalize.diag=T)
  GLT_time[i] = Sys.time()-t1
  t2 = Sys.time()
  gl_cv_F = glasso_cv(ts=xx, k=5, verbose=T, rholist=rho_seq, penalize.diag=F)
  GLF_time[i] = Sys.time()-t2
  
  GL_est_T_list[[i]] = gl_cv_T$wi
  GL_cov_T_list[[i]] = gl_cv_T$w
  rhomax_GLT[i] = gl_cv_T$rhomax
  GL_est_F_list[[i]] = gl_cv_F$wi
  GL_cov_F_list[[i]] = gl_cv_F$w
  rhomax_GLF[i] = gl_cv_F$rhomax
}

#Calculate Stein's loss, Frobenius norm
#sensitivity, specificity, precision, accuracy of the estimated precision
ind_nonzero = which(omega_elements!=0)
GL_T_loss = rep(NA,m); GL_T_fro = rep(NA,m)
sen_GLT = rep(NA,m); spe_GLT = rep(NA,m); fpr_GLT = rep(NA,m)
prec_GLT = rep(NA,m); acc_GLT = rep(NA,m)

for (i in 1:m) {
  Sigma_GL_T = GL_cov_T_list[[i]];   Omega_GL_T = GL_est_T_list[[i]]
  
  GL_T_loss[i] = log(det(Sigma_GL_T%*%sigma_inv))+sum(diag(Omega_GL_T%*%sigma))-p
  GL_T_fro[i] = norm(Omega_GL_T-sigma_inv,type="f")
  
  TP = sum((abs(t(Omega_GL_T)[lower.tri(Omega_GL_T)][ind_nonzero])!=0))
  FP = sum((abs(t(Omega_GL_T)[lower.tri(Omega_GL_T)][-ind_nonzero])!=0))
  TN = sum((abs(t(Omega_GL_T)[lower.tri(Omega_GL_T)][-ind_nonzero])==0))
  FN = sum((abs(t(Omega_GL_T)[lower.tri(Omega_GL_T)][ind_nonzero])==0))
  sen_GLT[i] = TP/sum(omega_elements!=0)
  spe_GLT[i] = TN/sum(omega_elements==0)
  prec_GLT[i] = TP/(TP+FP)
  acc_GLT[i] = (TP+TN)/(TP+TN+FP+FN)
  fpr_GLT[i] = 1-spe_GLT[i]
}

means = cbind(mean(GL_T_loss),mean(GL_T_fro),mean(sen_GLT),mean(fpr_GLT))
cat('GL1 mean: loss, Fnorm, TPR, FPR',round(means,3))
cat('\n')
sds = cbind(sd(GL_T_loss),sd(GL_T_fro),sd(sen_GLT),sd(fpr_GLT))
cat('GL1 sd: loss, Fnorm, TPR, FPR',round(sds,3))
cat('\n')
GL_F_loss = rep(NA,m); GL_F_fro = rep(NA,m)
sen_GLF = rep(NA,m); spe_GLF = rep(NA,m); fpr_GLF = rep(NA,m)
prec_GLF = rep(NA,m); acc_GLF = rep(NA,m)

for (i in 1:m) {
  Sigma_GL_F = GL_cov_F_list[[i]];   Omega_GL_F = GL_est_F_list[[i]]
  
  GL_F_loss[i] = log(det(Sigma_GL_F%*%sigma_inv))+sum(diag(Omega_GL_F%*%sigma))-p
  GL_F_fro[i] = norm(Omega_GL_F-sigma_inv,type="f")
  
  TP = sum((abs(t(Omega_GL_F)[lower.tri(Omega_GL_F)][ind_nonzero])!=0))
  FP = sum((abs(t(Omega_GL_F)[lower.tri(Omega_GL_F)][-ind_nonzero])!=0))
  TN = sum((abs(t(Omega_GL_F)[lower.tri(Omega_GL_F)][-ind_nonzero])==0))
  FN = sum((abs(t(Omega_GL_F)[lower.tri(Omega_GL_F)][ind_nonzero])==0))
  sen_GLF[i] = TP/sum(omega_elements!=0)
  spe_GLF[i] = TN/sum(omega_elements==0)
  prec_GLF[i] = TP/(TP+FP)
  acc_GLF[i] = (TP+TN)/(TP+TN+FP+FN)
  fpr_GLF[i] = 1-spe_GLF[i]
}

means = cbind(mean(GL_F_loss),mean(GL_F_fro),mean(sen_GLF),mean(fpr_GLF))
cat('GL2 mean: loss, Fnorm, TPR, FPR',round(means,3))
cat('\n')
sds = cbind(sd(GL_F_loss),sd(GL_F_fro),sd(sen_GLF),sd(fpr_GLF))
cat('GL2 sd: loss, Fnorm, TPR, FPR',round(sds,3))
cat('\n')


###########################################################

MCC_1 = rep(NA, m)

for (i in 1:m) {
  #Sigma_GL_T = GL_cov_T_list[[i]];  
  Omega_GL_T = GL_est_T_list[[i]]
  
  #GL_T_loss[i] = log(det(Sigma_GL_T%*%sigma_inv))+sum(diag(Omega_GL_T%*%sigma))-p
  #GL_T_fro[i] = norm(Omega_GL_T-sigma_inv,type="f")
  
  TP = as.numeric(sum((abs(t(Omega_GL_T)[lower.tri(Omega_GL_T)][ind_nonzero])!=0)))
  FP = as.numeric(sum((abs(t(Omega_GL_T)[lower.tri(Omega_GL_T)][-ind_nonzero])!=0)))
  TN = as.numeric(sum((abs(t(Omega_GL_T)[lower.tri(Omega_GL_T)][-ind_nonzero])==0)))
  FN = as.numeric(sum((abs(t(Omega_GL_T)[lower.tri(Omega_GL_T)][ind_nonzero])==0)))
  
  MCC_1[i] = ((TP*TN)-(FP*FN))/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  
  # sen_GLT[i] = TP/sum(omega_elements!=0)
  # spe_GLT[i] = TN/sum(omega_elements==0)
  # prec_GLT[i] = TP/(TP+FP)
  # acc_GLT[i] = (TP+TN)/(TP+TN+FP+FN)
  # fpr_GLT[i] = 1-spe_GLT[i]
}

###########################################################

MCC_2 = rep(NA, m)

for (i in 1:m) {
  #Sigma_GL_F = GL_cov_F_list[[i]];   
  Omega_GL_F = GL_est_F_list[[i]]
  
  #GL_F_loss[i] = log(det(Sigma_GL_F%*%sigma_inv))+sum(diag(Omega_GL_F%*%sigma))-p
  #GL_F_fro[i] = norm(Omega_GL_F-sigma_inv,type="f")
  
  TP = as.numeric(sum((abs(t(Omega_GL_F)[lower.tri(Omega_GL_F)][ind_nonzero])!=0)))
  FP = as.numeric(sum((abs(t(Omega_GL_F)[lower.tri(Omega_GL_F)][-ind_nonzero])!=0)))
  TN = as.numeric(sum((abs(t(Omega_GL_F)[lower.tri(Omega_GL_F)][-ind_nonzero])==0)))
  FN = as.numeric(sum((abs(t(Omega_GL_F)[lower.tri(Omega_GL_F)][ind_nonzero])==0)))
  
  MCC_2[i] = ((TP*TN)-(FP*FN))/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  # sen_GLF[i] = TP/sum(omega_elements!=0)
  # spe_GLF[i] = TN/sum(omega_elements==0)
  # prec_GLF[i] = TP/(TP+FP)
  # acc_GLF[i] = (TP+TN)/(TP+TN+FP+FN)
  # fpr_GLF[i] = 1-spe_GLF[i]
}
cat('MCC mean and SD for GL1 is:', round(mean(MCC_1),3), round(sd(MCC_1),3))
cat('\n')
cat('MCC mean and SD for GL2 is:', round(mean(MCC_2),3), round(sd(MCC_2),3))
cat('\n')
cat(round(mean(GLT_time),3))
cat('\n')
cat(round(mean(GLF_time),3))
