source("./R/Basics.R")
source("./R/Modeling.R")
source("./R/Simulation.R")
source("./R/poly_mirt_wrappers.R")
library(Rcpp)
sourceCpp("./Rcpp/poly_mirt_rcpp.cpp")

#####################################
########## Example 1: MGRM ##########
#####################################

# 1. Generate data
## 1.1 Setup
J=300 # number of items
sample_coef=25
N=J*sample_coef # sample size
K=5 # number of factors
C=5 # number of options
rho=c(0.1,0.3) # rho=0 for orthogonal case

## 1.2 Parameters and data
A_true <- generate_A(J,K,zero_prob=0.7)  # true loading matrix
tau_true <- generate_tau(J,C=C)
Data <- generate_data(N,J,rho=rho,A_true,theta_distr="normal",
                      model="GRM",C=C,tau_true=tau_true)
Y <- Data$Y

# 2. Three-stage estimation
## 2.1 The first stage (rotation)
Params <- GRM_HiMMIRT(Y,K,C,alpha=1,rotation=T,max_iter=200,tol=1e-3)
### tol=1e-3 as it's only a warm start

## 2.2 The second stage (lasso penalization)
lambda_max <- ifelse(sample_coef==10,J*2,J*3)
result_list <- GRM_autolasso_HiMIRT(Y,K,C,Params=Params,try_out=15,
                                    lambda_range=c(0,lambda_max),
                                    orthogonal=F) # fit with strong rule

df <- as.data.frame(t(sapply(1:(length(result_list)-1),
                             function(i) calSelect(result_list[[i]],
                                                   Data, N, K))))
BIC_num <- which.min(df$BIC)
Lambda_BIC <- ifelse(result_list[[BIC_num]]$A==0,0,1). # design matrix

result_BIC <- GRM_CHiMIRT(Y,Params,Lambda_BIC,K,
                          orthogonal=F,
                          max_iter=200,tol=1e-5)

## 2.3 The third stage (cardinality refinement)
result_seq <- GRM_card_HiMIRT(Y,K,C,Params=result_BIC,
                              random_start=F,
                              kappaC=sum(Lambda_BIC!=0),
                              orthogonal=F, tol=1e-5,
                              alpha=1,max_iter=200)

repeat {
  small_num <- sum(abs(result_seq$A)<.3 & result_seq$A!=0)
  kappaC_new <- sum(result_seq$A!=0)-small_num
  if (small_num==0 || kappaC_new<K) break
  result_seq <- GRM_card_HiMIRT(Y,K,C,Params=result_seq,
                                random_start=F,
                                kappaC=kappaC_new,
                                orthogonal=F, tol=1e-5,
                                alpha=1,max_iter=200)
}

# 3. Check performance
calResults(result_seq,Data,N,J,K,C,model="GRM")
calResults(result_BIC,Data,N,J,K,C,model="GRM")
### result_BIC usually has lower specificity
### one can also check rotation's performance by controlling tol=1e-5. It's
### unstable when K is large.

######################################
########## Example 2: MGPCM ##########
######################################

# 1. Generate data
## 1.1 Setup
J=300
sample_coef=25
N=J*sample_coef
K=5
C=5
rho=c(0.1,0.3)

## 1.2 Parameters and data
A_true <- generate_A(J,K,zero_prob=0.7)  # true loading matrix
beta_true <- generate_beta(J,C=C)
Data <- generate_data(N,J,rho=rho,A_true,theta_distr="normal",
                      model="GPCM",C=C,beta_true=beta_true)
Y <- Data$Y

# 2. Three-stage estimation
## 2.1 The first stage (rotation)
Params <- GPCM_HiMMIRT(Y,K,C,rotation=T,
                       max_iter=200,tol=1e-3)
### tol=1e-3 as it's only a warm start

## 2.2 The second stage (lasso penalization)
lambda_max <- ifelse(sample_coef==10,J*5,J*12)
result_list <- GPCM_autolasso_HiMIRT(Y,K,C,Params=Params,try_out=15,
                                     lambda_range=c(0,lambda_max),
                                     orthogonal=F) # fit with strong rule

df <- as.data.frame(t(sapply(1:(length(result_list)-1),
                             function(i) calSelect(result_list[[i]],
                                                   Data, N, K))))
BIC_num <- which.min(df$BIC)
Lambda_BIC <- ifelse(result_list[[BIC_num]]$A==0,0,1) # design matrix

result_BIC <- GPCM_CHiMIRT(Y,Params,Lambda_BIC,K,
                           orthogonal=F,
                           max_iter=200,tol=1e-5)
### set tol=1e-3 * K / (sample_coef * J) here for maximal estimation accuracy

## 2.3 The third stage (cardinality refinement)
result_seq <- GPCM_card_HiMIRT(Y,K,C,Params=result_BIC,
                               random_start=F,
                               kappaC=sum(Lambda_BIC!=0),
                               orthogonal=F, tol=1e-5,
                               max_iter=200)

repeat {
  small_num <- sum(abs(result_seq$A)<.3 & result_seq$A!=0)
  kappaC_new <- sum(result_seq$A!=0)-small_num
  if (small_num==0 || kappaC_new<K) break
  result_seq <- GPCM_card_HiMIRT(Y,K,C,Params=result_seq,
                                random_start=F,
                                kappaC=kappaC_new,
                                orthogonal=F, tol=1e-5,
                                max_iter=200)
}

# 3. Check performance
calResults(result_seq,Data,N,J,K,C,model="GPCM")
calResults(result_BIC,Data,N,J,K,C,model="GPCM")
### result_BIC usually has lower specificity
### one can also check rotation's performance by controlling tol=1e-5. It's
### unstable when K is large.
