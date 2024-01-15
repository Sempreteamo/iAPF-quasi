library(mvtnorm)
library(MASS)
library(FKF)
library(profvis)
library(mvnfast)

##parameters
start = proc.time()
set.seed(1)
Num <- 200 #total number of particles
N <- vector()
N[1] <- Num
Time = 100
Lag = 8 #lag can be any integers <= Time which is divided by Time
alpha = 0.42
d = 35
k <- 5
tau <- 0.5
kappa = 0.5
A <- matrix(nrow = d, ncol = d)
for (i in 1:d){
  for (j in 1:d){
    A[i,j] = alpha^(abs(i-j) + 1)
  }
}
B = C = D = diag(1, nrow = d, ncol = d)
X_true <-  matrix(0, nrow = Time, ncol = d )
obs <-  matrix(0, nrow = Time, ncol = d )
dt <- ct <- matrix(0, d, 1)
Tt <- A
P0 <- Zt <- Ht <- Gt <- diag(1, d, d)
a0 <- rep(0, d)
index <- 1
psi_final <- matrix(NA, nrow = Time, ncol = 2*d)
MH_distance <- vector()
KS_distance <- vector()

Obs <- function(){
  #set.seed(123)
  X_true[1,] <- rnorm(d)     
  for(t in 2:Time){ 
    #set.seed(123)#observations
    X_true[t,] <- rnorm(d) + A%*%X_true[t-1,]  #t(rnorm(d) + A%*%x)
  }
  #set.seed(123)
  return(matrix(rnorm(Time*d, X_true, 1), ncol = d))
}
obs <- Obs()

X <- array(NA, dim = c(Time, Num, d))
X_ <- array(NA, dim = c(Time, Num, d))
w <- matrix(NA, Time, Num)
Z = 0
Z_apf <- vector()

fkf.obj <- fkf(a0, P0, dt, ct, Tt, Zt, Ht, Gt, yt = t(obs))
fks.obj <- fks(fkf.obj)
fkf.obj_Z <- fkf(a0, P0, dt, ct, Tt, Zt, Ht, Gt, yt = t(obs))$logLik #???

Smoothing <- function(specific_time){
  fks_mean <- fks.obj$ahatt[,specific_time]
  fks_cov <- fks.obj$Vt[,,specific_time]
  return(list(fks_mean, fks_cov))
}

##Start our iAPF-quasi algorithm

source('iapf_only_controlled_new_ed2_nd.R')

####Initialization####
iAPF1 <- function(time_step, L) {
  #initialization for time t = 1...L to obtain psi1, w1, X1
  L <- L
  output <- Init(time_step, L)
  X <- output[[1]] 
  w <- output[[2]]
  psi <- output[[3]]
  
  #choose psi from 1 to 3L/4
  psi_final <- psi[complete.cases(psi), ][1:ceiling(3*L/4),]
  
  #Iteration
  #obtain 
  if(L < Time){
    for(time_step in seq(2*L,Time,L)){
      
      #I didn't include any resampling in this step
      #Run iAPF with the initial distribution we defined 
      
      output <- init_APF(time_step, w, X, L)
      X_apf <- output[[1]]
      w_apf <- output[[2]]
      Z_apf <- output[[3]]
      
      output2 <- psi_APF(time_step, X_apf, Z_apf, w, X, L)
      
      #smoothing particles
      X <- output2[[1]]
      w <- output2[[2]]
      psi <- output2[[3]]
      gap_matrix <- matrix(NA, nrow = L/2, ncol = 2*d)
      
      #take psi from kL+L/4+1 to kL+3L/4+1
      psi_final <- rbind(psi_final, gap_matrix, psi[complete.cases(psi), ][ceiling(L/4 + 1):ceiling(3*L/4),])
    }  
  }
  return(psi_final)
}


psi_final <- iAPF1(Lag, Lag)
####the second layer of iAPF from t=1 to L to obtain psi####

iAPF2 <- function(psi_final, time_step, L1, L) {
#t = 1, . . . , L/2 to obtain X2, w2 and psi2
  psi_final <- psi_final
  L1 <- L1
  L <- L
  
output <- Init(time_step, L1) #pass
X <- output[[1]] 
w <- output[[2]]
psi <- output[[3]]    

####Algorithm####
if(L < Time){
  #start from L/2 to 3L/2 to obtain psi2
  count = 0
  for(time_step in seq(ceiling(3*L/2),Time,L)){
    #I didn't include any resampling in this step
    #Run iAPF with the initial distribution we defined 
    
    output <- init_APF(time_step, w, X, L)
    X_apf <- output[[1]]
    w_apf <- output[[2]]
    Z_apf <- output[[3]]
    
    output2 <- psi_APF(time_step, X_apf, Z_apf, w, X, L)
    
    #smoothing particles
    X <- output2[[1]]
    w <- output2[[2]]
    psi <- output2[[3]]
    
    psi_final[ceiling(3*L/4+1+count*L):ceiling(3*L/4+count*L+ L/2),] <- psi[complete.cases(psi), ][ceiling(L/4 + 1):ceiling(3*L/4),]
    
    count = count + 1
  }  
}

if(time_step != Time){
  time_step <- Time
  
  output <- init_APF(time_step, w, X, L1)
  X_apf <- output[[1]]
  w_apf <- output[[2]]
  Z_apf <- output[[3]]
  
  output2 <- psi_APF(time_step, X_apf, Z_apf, w, X, L1)
  
  #smoothing particles
  X <- output2[[1]]
  w <- output2[[2]]
  psi <- output2[[3]]
  psi_final <- rbind(psi_final, psi[complete.cases(psi), ][ceiling(L1/2 + 1):L1,])
}
  return(psi_final)
}

psi_final <- iAPF2(psi_final, Lag/2, Lag/2, Lag)

####psi-APF####
#purely smoothing particles and normalizing constant
smoothing_APF <- function(psi_pa, N, Time){ 
  #l >= 2
  X_apf <- array(NA, dim = c(Time, N[1], d))
  w_apf <- matrix(NA, Time, N[1])
  Z_apf <- 0
  
  re = 0
  
  #when n = kL, we use the new mu.tilda to initialize particles
  
  X_apf[1,,] <- mu_aux(psi_pa, 1, N, 1)
  for(i in 1:Num){
    w_apf[1,i] <- g_aux_smoo(obs[1], X_apf[1,i,],1, psi_pa, Time)
  }
  
  for(t in 2:Time){
    
    if(ESS(t,w_apf, is.log = TRUE) <= kappa*Num){
      re = re + 1
      mx <- max(w_apf[t-1,])
      w_ <- exp(w_apf[t-1,1:Num]-mx)/sum(exp(w_apf[t-1, 1:Num] - mx))
      Z_apf = Z_apf + log(mean(exp(w_apf[t-1,]-mx))) + mx
      #set.seed(123)
      mix <- sample(1:Num,Num, replace = TRUE, prob = w_)
      X_apf <- X_apf[, mix,, drop = FALSE]
      
      for(i in 1:Num){
        #filtering particles
        X_apf[t,i,] <- f_aux(X_apf[t-1, i,], psi_pa, t)
        w_apf[t,i] <- g_aux_smoo(obs[t,], X_apf[t,i,], t, psi_pa, Time) 
        
      }
    }else{
      
      for(i in 1:Num){
        #filtering particles
        X_apf[t,i,] <- f_aux(X_apf[t-1,i,], psi_pa, t) 
        w_apf[t,i] <- w_apf[t-1,i] + g_aux_smoo(obs[t,], X_apf[t,i,], t, psi_pa, Time)
      }
    }
    
  }
  
  mx <- max(w_apf[Time, 1:Num])
  Z_apf <- Z_apf + log(mean(exp(w_apf[Time, 1:Num]-mx))) + mx
  
  return(list(X_apf, w_apf, Z_apf, re))
}

output2 <- smoothing_APF(psi_final, N, Time)
X <- output2[[1]]
w <- output2[[2]]
Z <- output2[[3]]
re <- output2[[4]] #above are not the same

mx <- max(w[Time,]) 
w_ <- exp(w[Time,]-mx)/sum(exp(w[Time,] - mx))

#Mahalanobis distance
for(specific in 1:Time){
  weighted_mean <- colSums(w_*X[specific,,])
  output_s <- Smoothing(specific)
  fks_mean <- output_s[[1]]
  fks_cov <- output_s[[2]]
  MH_distance[specific] <- mahalanobis(weighted_mean, fks_mean, fks_cov)/d
}

plot(MH_distance)
normalizing_c <- exp(Z-fkf.obj_Z)


#Calculation of the mean and variance of empirical distribution
#Here 'weighted_mean' is the empirical mean, 'variance' is the empirical variance
Empirical <- function(X, w_, specific_time){
  weighted_mean <- sum(w_*X[specific_time,])
  s_mean <- sum(w_*X[specific_time,]^2)
  variance <- s_mean - weighted_mean^2
  return(list(weighted_mean, variance))
}

#output_e <- Empirical(X, w_, specific)
#weighted_mean <- output_e[[1]]
#variance <- output_e[[2]]
#empirical distribution function
#the equation of dKS(F, G) below is in the notes

dKS <- function(X, w_, specific_time) {
  output_s <- Smoothing(specific_time)
  fks_mean <- output_s[[1]]
  fks_var <- output_s[[2]]
  
  #the replacement of the weights and particles matrices
  X_ <- X
  w_c <- w_
  
  d <- vector()
  order <- order(X[specific_time,,])
  X_[specific_time,,] <- sort(X[specific_time,,])
  w_c <- w_[order]
  cumsum_w <- cumsum(w_c)
  f <- pnorm(X_[specific_time,,], mean = fks_mean, sd = sqrt(fks_var)) 
  d <- abs(f - cumsum_w)
  return(max(d))
}.

#KS_distance <- dKS(X, w_, specific)
normalizing_c <- exp(Z-fkf.obj_Z)
#KS_distance_start <- dKS(X, w_, specific1)

for(specific in 1:Time){
  KS_distance[specific] <- dKS(X, w_, specific)
}

