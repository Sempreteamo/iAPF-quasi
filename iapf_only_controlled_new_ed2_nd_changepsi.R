#ess resampling
ESS <- function(t,w, is.log=FALSE){
  if(is.log) {
    mx <- max(w[t-1,])
    s <- sum(exp(w[t-1,]-mx))
    ess <- 1/sum((exp(w[t-1,]-mx)/s)^2)
  }else{
    s <- sum(w[t-1,])
    ess <- 1/sum((w[t-1,]/s)^2) 
  }
  return(ess)  
}

residual <- function(t, w){
  mx <- max(w[t-1,])
  w_ <- exp(w[t-1,] - mx)/sum(exp(w[t-1,] - mx))
  
  Ntm <- as.integer(Num*w_)
  
  mix <- unlist(lapply(1:Num, function(i) {rep(i, Ntm[i])}))
  mr <- Num - sum(Ntm)
  
  w_hat <- w_ - Ntm/Num
  w_hat <- w_hat*Num/mr
  
  mix <- c(sample(1:Num, mr, replace = TRUE, prob = w_hat), mix)
  
  return(mix)
  
} 

Num_apf <- function(Z, l, k){
  return(sd(exp(Z[max(l-k,1):l]-max(Z[max(l-k,1):l])))/mean(exp(Z[max(l-k,1):l]-max(Z[max(l-k,1):l]))))
}

#trans prob
#N(; Ax, B)
f <- function(x){
  return (rnorm(d) + A%*%x)   #trans prob
}

#obs prob  
#N(; Cx, D)
g <- function(y, x){  
  return (log((2*pi)^(-d/2)) + (-1/2)*t(y-x)%*%(y-x)) #obs prob  C%*%x = x
}

#twisted mu
mu_aux <- function(psi_pa, l, N, t){  
  return(rmvn(N[l], diag(((psi_pa[t, (d+1):(d+d)])^(-1)+1)^(-1), nrow=d,ncol=d)%*%
                (diag((psi_pa[t, (d+1):(d+d)])^(-1), nrow=d,ncol=d)%*%psi_pa[t,1:d]), 
              diag(((psi_pa[t, (d+1):(d+d)])^(-1)+1)^(-1), nrow=d,ncol=d)))
}

#twisted g
g_aux <- function(y, x, t, psi_pa, n, L){  
  if(t == (n-L+1)){
    return(g(y, x) + psi_tilda(x, psi_pa, t, n) + log((2*pi)^(-d/2)*det(diag(psi_pa[t, (d+1):(d+d)]+1, nrow=d,ncol=d))^
                                                        (-1/2)*exp((-1/2)*t(-psi_pa[t, 1:d])%*%
                                                                     diag((psi_pa[t, (d+1):(d+d)]+1)^(-1), nrow=d,ncol=d)%*%
                                                                     (-psi_pa[t, 1:d]))) - psi_t(x, psi_pa, t, n))  #initialisation of g = t=1 or t=L?
  }else{
    return(g(y, x) + psi_tilda(x, psi_pa, t, n) - psi_t(x, psi_pa, t, n))  #g_2:T 
  }
}

g_aux_smoo <- function(y, x, t, psi_pa, n){  
  if(t == 1){
    return(g(y, x) + psi_tilda(x, psi_pa, t, n) + log((2*pi)^(-d/2)*det(diag(psi_pa[t, (d+1):(d+d)]+1, nrow=d,ncol=d))^
                                                        (-1/2)*exp((-1/2)*t(-psi_pa[t, 1:d])%*%
                                                                     diag((psi_pa[t, (d+1):(d+d)]+1)^(-1), nrow=d,ncol=d)%*%
                                                                     (-psi_pa[t, 1:d]))) - psi_t(x, psi_pa, t, n))  #initialisation of g = t=1 or t=L?
  }else{
    return(g(y, x) + psi_tilda(x, psi_pa, t, n) - psi_t(x, psi_pa, t, n))  #g_2:T 
  }
}

#twisted f
f_aux <- function(x, psi_pa, t){
  return(rmvn(1, diag(((psi_pa[t, (d+1):(d+d)])^(-1)+1)^(-1), nrow=d,ncol=d)%*%
                (A%*%x + diag(psi_pa[t, (d+1):(d+d)]^(-1), nrow=d,ncol=d)%*%psi_pa[t,1:d]), 
              diag(((psi_pa[t, (d+1):(d+d)])^(-1)+1)^(-1), nrow=d,ncol=d)))  #f_2:T 
}

#psi.tilda_t = f (xt, ψt+1); psi.tilda_n = 1
psi_tilda <- function(x, psi_pa, t, n){  #from 0 to T. 0,T = 1 
  if(t == n){
    psi_tilda <- 0
  }else{   #psi_pa_t = psi_t
    psi_tilda <- log((2*pi)^(-d/2)) + log(det(diag(psi_pa[t+1, (d+1):(d+d)]+1, nrow=d, ncol=d))^(-1/2)) +
      (-1/2)*t(A%*%x - psi_pa[t+1, 1:d])%*%diag((psi_pa[t+1, (d+1):(d+d)]+1)^(-1), nrow=d,ncol=d)%*%
      (A%*%x-psi_pa[t+1, 1:d]) 
  }
  return(psi_tilda)
}

#ψt(xt) = N (xt; m_t, Σ_t), m_t, Σ_t obtained in Psi function
psi_t <- function(x, psi_pa, t, n){ 
  if(t == (n + 1)){
    psi_t <- 0
  }else{
    psi_t <- log((2*pi)^(-d/2)) + log(det(diag(psi_pa[t, (d+1):(d+d)], nrow=d,ncol=d))^(-1/2)) +						
      (-1/2)*t(x-psi_pa[t, 1:d])%*%diag((psi_pa[t, (d+1):(d+d)])^(-1), nrow=d,ncol=d)%*%						
      (x-psi_pa[t, 1:d])
  }
  return(psi_t)
}

#the distribution here and the distribution below should be double checked

#the initial distribution mu at time kL
# ∑W_kL(?)^i*f(X_{k−1}^i,(k−1)L, ·), i in 1:N
#I used the weights at time (k-1)L because these are the closest weights, 
#sample particles X_{kL-L+1} using the particles at time X_(k-1)L

change_mu <- function(n, w, X, L){
  sam <- matrix(0,Num,d)
  mx <- max(w[n-L,])
  w_ <- exp(w[n-L,]-mx)/sum(exp(w[n-L,] - mx))
  s <- sample(1:Num, size = Num, replace = TRUE, prob = w_) 
  mus <- X[n-L,,]
  for(i in 1:Num){
    sam[i,] <- rnorm(d) + A%*%mus[s[i],]
  }
  return(list(sam, w_))
}

#the initial distribution mu.tilda at time kL
#the particles and weights I use are exactly the same as above, and I use
#f.tilda to generate particles instead of f
change_mupsi <- function(n, w, X, psi_pa, t, N, l, L){
  sam <- matrix(0,Num,d)
  w_adj <- vector()
  
  #here we need to adjust the weights 
  
  for(i in 1:Num){
    w_adj[i] <- w[n-L,i]*exp(1/2*(t(A%*%X[n-L,i,]) + t(psi_pa[t,1:d])%*%
                                    diag(psi_pa[t, (d+1):(d+d)]^(-1), nrow=d,ncol=d))%*%diag((psi_pa[t, (d+1):(d+d)]^(-1) + 1)^(-1), nrow=d,ncol=d)%*%
                               (A%*%X[n-L,i,] + diag(psi_pa[t, (d+1):(d+d)]^(-1), nrow=d,ncol=d)%*%psi_pa[t,1:d]) - 1/2*(t(A%*%X[n-L,i,])%*%A%*%X[n-L,i,] +
                                                                                                                           t(psi_pa[t,1:d])%*%diag(psi_pa[t, (d+1):(d+d)]^(-1), nrow=d,ncol=d)%*%psi_pa[t,1:d]))
  }
  
  mx <- max(w_adj)
  w_ <- exp(w_adj-mx)/sum(exp(w_adj - mx))
  s <- sample(1:Num, size = Num, replace = TRUE, prob = w_) 
  mus <- X[n-L,,]
  for(i in 1:N[l]){
    sam[i,] <- f_aux(mus[s[i],], psi_pa, t)
    
  }
  return(list(sam, w_))
}

####iapf
####Initialization####
Init <- function(time_step, L){
  L <- L
  n <- time_step
  output <- init_APF(n, w = 0, X = 0, L)
  X_apf <- output[[1]]
  w_apf <- output[[2]]
  Z_apf <- output[[3]]
  
  output2 <- psi_APF(n, X_apf, Z_apf, w = 0, X = 0, L)
  X <- output2[[1]]
  w <- output2[[2]]
  psi <- output2[[3]]
  Z <- output2[[4]]
  
  return(list(X, w, psi, Z))
}

#Below is the iAPF algorithm
#First we do the initialization
####init_APF####
init_APF <- function(n, w, X, L){  #pure filtering
  l = 1
  Z_apf[1] = 0
  X_init <- array(NA, dim = c(Time, Num, d))
  w_init <- matrix(NA, Time, N[l])
  
  #when n = L, we use the mu as the initialization distribution mu; 
  #when n = kL, we use the new distribution
  #I tried to modify the previous paths from 1 to n-L+1
  #X_init[1:(n-L+1),,] <- rnorm(N[l]*d)  
  #for(i in 1:N[l]){
   # w_init[1:(n-L+1),i] <- g(obs[n-L+1,], X_init[n-L+1,i,])  
  #}
  
  if(n == L){
    X_init[1:(n-L+1),,] <- rnorm(N[l]*d)  
    for(i in 1:N[l]){
      w_init[1:(n-L+1),i] <- g(obs[n-L+1,], X_init[n-L+1,i,])  
    }
  }else{
    output <- change_mu(n, w, X, L)
    X_init[1:(n-L+1),,] <- output[[1]]
    w_ <- 1
    
    for (i in 1:Num){
      w_init[1:(n-L+1), i] <- log(sum(w_*dmvn(X_init[n-L+1,,], A%*%X[n-L,i,], B)))
    }
  }
  
  for(t in (n-L+2):n){
    
    if(ESS(t, w_init, is.log=TRUE) <= kappa*N[l]){
      mx <- max(w_init[t-1,])
      w_ <- exp(w_init[t-1,] - mx)/sum(exp(w_init[t-1,] - mx))
      Z_apf[1] = Z_apf[1] + log(mean(exp(w_init[t-1,] - mx))) + mx
      mix <- sample(1:N[l], N[l], replace = TRUE, prob = w_)
      #mix <- residual(t, w_init)
      # at the initialization stage, we want filtering particles for psi
      
      for(i in 1:N[l]){
        X_init[t,i,] <- f(X_init[t-1,mix[i],]) 
        w_init[t,i] <- g(obs[t,], X_init[t,i,])  
      }
      
    }else{
      for(i in 1:N[l]){
        X_init[t,i,] <- f(X_init[t-1,i,]) 
        w_init[t,i] <- w_init[t-1,i] + g(obs[t,], X_init[t,i,])  
      }
    }
  }
  
  mx <- max(w_init[n, 1:N[l]])
  Z_apf[1] <- Z_apf[1] + log(mean(exp(w_init[n, 1:N[l]]-mx))) + mx
  return(list(X_init = X_init, w_init = w_init, Z_apf))
}

#Within the iAPF scheme, when the number of iterations l >= 2, run APF
####APF####
APF <- function(n, w, X, psi_pa, l, Z_apf, N, L){ #purely filtering particles
  #l >= 2
  X_apf <- array(NA, dim = c(Time, N[l], d))
  w_apf <- matrix(NA, Time, N[l])
  Z_apf[l] <- 0
  
  #when n = kL, we use the new mu.tilda to initialize particles
  #X_apf[1:(n-L+1),1:N[l],] <- mu_aux(psi_pa, l, N, n-L+1)
  #for(i in 1:N[l]){
   # w_apf[1:(n-L+1),i] <- g_aux(obs[n-L+1,], X_apf[n-L+1,i,],n-L+1, psi_pa, n, L) 
  #}
  if(n == L){
    X_apf[1:(n-L+1),1:N[l],] <- mu_aux(psi_pa, l, N, n-L+1)
    for(i in 1:N[l]){
      w_apf[1:(n-L+1),i] <- g_aux(obs[n-L+1,], X_apf[n-L+1,i,],n-L+1, psi_pa, n, L) 
    }
  }else{
    output <- change_mupsi(n, w, X, psi_pa, n-L+1, N, l, L)
    X_apf[1:(n-L+1),1:N[l],] <- output[[1]]
    w_ <- 1
    
    value <- t(A%*%t(X[n-L,,]))
    cov_matrix <-  diag(((psi_pa[n-L+1, (d+1):(d+d)])^(-1)+1)^(-1), nrow=d,ncol=d)
    mean_d1 <- diag(((psi_pa[n-L+1, (d+1):(d+d)])^(-1)+1)^(-1), nrow=d,ncol=d)
    mean_d2 <- diag(psi_pa[n-L+1, (d+1):(d+d)]^(-1), nrow=d,ncol=d)%*%psi_pa[n-L+1,1:d]
    
    for (i in 1:N[l]){
      w_apf[1:(n-L+1), i] <- log(sum(w_*dmvn(X_apf[n-L+1,,], mean_d1%*%(A%*%t(X[n-L,,])[,i] + mean_d2), cov_matrix)))
    }
  }
  
  for(t in (n-L+2):n){
    
    if(ESS(t,w_apf, is.log = TRUE) <= kappa*N[l]){
      
      mx <- max(w_apf[t-1,])
      w_ <- exp(w_apf[t-1,1:N[l]]-mx)/sum(exp(w_apf[t-1, 1:N[l]] - mx))
      Z_apf[l] = Z_apf[l] + log(mean(exp(w_apf[t-1,]-mx))) + mx
      mix <- sample(1:N[l],N[l], replace = TRUE, prob = w_)
      #mix <- residual(t, w_apf)
      
      for(i in 1:N[l]){
        #filtering particles
        X_apf[t,i,] <- f_aux(X_apf[t-1, mix[i],], psi_pa, t)
        w_apf[t,i] <- g_aux(obs[t,], X_apf[t,i,], t, psi_pa, n, L) 
      }
    }else{
      
      for(i in 1:N[l]){
        #filtering particles
        X_apf[t,i,] <- f_aux(X_apf[t-1,i,], psi_pa, t) 
        w_apf[t,i] <- w_apf[t-1,i] + g_aux(obs[t,], X_apf[t,i,], t, psi_pa, n, L)
      }
    }
    
  }
  mx <- max(w_apf[n, 1:N[l]])
  Z_apf[l] <- Z_apf[l] + log(mean(exp(w_apf[n, 1:N[l]]-mx))) + mx
  
  return(list(X_apf, w_apf, Z_apf))
}

####Psi####
Psi <- function(l, n, X_apf, N, L){
  psi <- matrix(NA, nrow = Time, ncol = N[l])
  psi_pa <- matrix(NA, nrow = Time, ncol = 2*d)
  
  #calculate psi
  for(t in n:(n-L+1)){
    if(t == n){
      for(i in 1:N[l]){
        psi[t,i] <- (1 / ((2 * pi)^(d / 2))) * 
          exp(-0.5 * t(X_apf[t,i,] - obs[t,]) %*% (X_apf[t,i,] - obs[t,]))
      }
      
      
    }else{
      for(i in 1:N[l]){
        
        psi[t,i] <- exp(g(obs[t,],X_apf[t,i,]))*dmvn(as.vector(A%*%X_apf[t,i,]), 
                                                     psi_pa[t+1, 1:d], diag(psi_pa[t+1, (d+1):(d+d)]+1, nrow=d,ncol=d))
      }
    }
    
    
    #calculate psi_t
    fn <- function(x, X_apf, psi){
      #lambda <- vector()
      #for(i in 1:N[l]){
      # lambda <-  2*sum((1 / ((2 * pi)^(d / 2) * sqrt(prod((x[(d+1):(d+d)])))) * 
      #                    exp(-0.5 * t(X_apf[t,i,] - x[1:d]) %*% 
      #                         diag(x[(d+1):(d+d)]^(-1), nrow=d,ncol=d) %*% (X_apf[t,i,] - x[1:d]))%*%
      #                  psi[t,1:N[l]]))/sum(psi[t,1:N[l]]^2) #2* or not 2*?
      #}
      lambda <-  2*sum(dmvn(X_apf[t,1:N[l],],x[1:d],
                            diag(exp(x[(d+1):(d+d)]), nrow=d,ncol=d))%*%psi[t,1:N[l]])/sum(psi[t,1:N[l]]^2) #2* or not 2*?
      return(sum((psi[t,1:N[l]] - (1/lambda)*dmvn(X_apf[t,1:N[l],],
                                                  x[1:d], diag(exp(x[(d+1):(d+d)]), nrow=d,ncol=d)))^2))
    }
    
    #get the distribution of psi_t
    if(t == n){
      psi_pa[t,] <- optim(par = c(colMeans(X_apf[t,1:N[l],]), rep(1, d)),
                          fn = fn, X_apf = X_apf, psi = psi)$par
    }else{
      
      psi_pa[t,] <- optim(par = c(X_apf[t,which.max(psi[t,1:N[l]]),], rep(1, d)), 
                          fn = fn, X_apf = X_apf, psi = psi)$par
      #optim(par = c(X_apf[t,which.max(psi[t,1:N[l]]),], rep(1, d)), 
      #                fn = fn, X_apf = X_apf, psi = psi, method='L-BFGS-B',
      #              lower=c(rep(-Inf, d),rep(0.3, d)), upper=rep(Inf, 2*d))$par
      #c(X_apf[t,which.max(psi[t,1:N[l]]),], rep(1, d))
      
    }
    psi_pa[t, (d+1):(d+d)] <- exp(psi_pa[t, (d+1):(d+d)])
    
    #print(psi_pa[t, 1])
    #print(obs[t])
  }
  return(psi_pa)
}

#The main part of iAPF to iterate and terminate
####psi_APF####
psi_APF <- function(n, X_apf, Z_apf, w, X, L){
  l = 1
  
  while(TRUE){
    
    output <- list()
    
    if(l != 1){
      #generate filtering particles X_apf for psi the next iteration
      #APF outputs filtering X_apf for the next psi, and smoothing X_apf_s
      #for the final calculation
      
      output <- APF(n, w, X, psi_pa, l, Z_apf, N, L)
      X_apf <- output[[1]]
      w_apf <- output[[2]]
      Z_apf <- output[[3]]
    }
    
    #to speed up the algorithm, I just fix the number of iterations to be k.
    #Here k = 5
    
    if(l <= k ){
      
      #receive filtering particles X_apf for psi
      psi_pa <- Psi(l, n, X_apf, N, L) 
      
      if(l > k & N[max(l-k,1)] == N[l] & is.unsorted(Z_apf[max(l-k,1):l])){  
        N[l+1] <- 2*N[l]
        
      }else{
        N[l+1] <- N[l]
      }
      
      l <- l+1
    }else break
  }
  
  #output psi
  return(list(X_apf, w_apf, psi_pa, Z_apf[l]))
}
