A=Sys.time()
library(mvtnorm)
library(Matrix)
library(psych)
Index<-as.integer(Sys.getenv("PBS_ARRAYID"))
set.seed(123+Index)
Durchlaeufe=40
# Necessary Functions
maxna <- function(x){
  return(max(x,na.rm = TRUE))
}
meanna <- function(x){
  return(mean(x, na.rm = TRUE))
}
MWurzel<-function(X){
  SVD <- svd(X)
  WurzelX <- SVD$u%*%(tcrossprod(sqrt(diag(SVD$d)),(SVD$v)))
  return(WurzelX)
}
# misk: Number of dependent replicates in group i at time s of subject k
misk_fun <- function(data, i,s,k){
  return(length(data[[i]][[s]][[k]]))
}
# mis: Number of observed replicates within group i at time s
mis_fun <- function(lambda,misk,ni,i,s){
  sum <- 0 
  for(k in 1:ni[i]){
    sum <- sum + lambda[[i]][k,s]*misk[[i]][k,s]
  }
  return(sum)
}
F_is <- function(data,lambda,lambdais,misk, mis,x,i,s,w,ni){
  sum <- 0
  if(w == 1){
    for(k in 1:ni[i]){
      wisk <- 1/(lambdais[i,s]*misk[[i]][k,s])
      for(u in 1:misk[[i]][k,s]){
        sum <- sum + wisk *lambda[[i]][k,s]*(x>data[[i]][[s]][[k]][u]) + wisk * 1/2*lambda[[i]][k,s]*(x == data[[i]][[s]][[k]][u])
      }
    }
    return(sum)
  }
  if(w == 2){
    wisk <- 1/mis[i,s]
    for(k in 1:ni[i]){
      for(u in 1:misk[[i]][k,s]){
        sum <- sum + wisk *lambda[[i]][k,s]*(x>data[[i]][[s]][[k]][u]) + wisk * 1/2*lambda[[i]][k,s]*(x == data[[i]][[s]][[k]][u])
      }
    }
    return(sum)
  }
}
# Empirical Distribution Function Paavo
F_isP <- function(data,lambda,lambdais,misk, mis,x,i,s,w,ni){
  Ergebnis <- rep(0,ni[i])
  if(w == 1){
    wisklambda <- 1/(lambdais[i,s]*misk[[i]][,s])*lambda[[i]][,s]
    for(k in 1:ni[i]){
      Ergebnis[k]= wisklambda[k]*sum((x>data[[i]][[s]][[k]])+1/2*(x == data[[i]][[s]][[k]]))
      
    }
    
  }
  if(w == 2){
    wisklambda <- 1/mis[i,s]*lambda[[i]][,s]
    for(k in 1:ni[i]){
      Ergebnis[k]= wisklambda[k]*sum((x>data[[i]][[s]][[k]]) +  1/2*(x == data[[i]][[s]][[k]]))
    }
  }
  return(sum(Ergebnis))}
# Relative effect estimator
p_is <- function(data, lambda,lambdais, i,s,w,ni,misk,mis,a,d){
  sum <- 0 
  if(w == 1){
    for(j in 1:a){
      for(t in 1:d){
        for(k in 1:ni[i]){
          wisk <- 1/(lambdais[i,s]*misk[[i]][k,s])
          for(u in 1:misk[[i]][k,s]){
            sum <- sum + lambda[[i]][k,s] * wisk * F_isP(data, lambda,lambdais, misk, mis, data[[i]][[s]][[k]][u], j,t,w,ni)
            
          }
        }
      }
    }
    return(sum/(a*d))
  }
  if(w == 2){
    wisk <- 1/mis[i,s]
    for(j in 1:a){
      for(t in 1:d){
        for(k in 1:ni[i]){
          for(u in 1:misk[[i]][k,s]){
            sum <- sum + lambda[[i]][k,s] * wisk * F_isP(data, lambda,lambdais, misk, mis, data[[i]][[s]][[k]][u], j,t,w,ni)
          }
        }
      }
    }
    return(sum/(a*d))
  }
}
# Pairwise relative effect estimator
p_is_ht <- function(data, lambda, lambdais, i,s,h,t,w,ni,misk,mis,a,d){
  if(w == 1){
    sum <- 0
    for(k in 1:ni[h]){
      wisk <- 1/(lambdais[h,t]*misk[[h]][k,t])
      for(u in 1:misk[[h]][k,t]){
        sum <- sum + lambda[[h]][k,t]*wisk*F_isP(data, lambda,lambdais, misk, mis, data[[h]][[t]][[k]][u], i,s,w,ni)
      }
    }
    return(sum)
  }
  if(w == 2){
    sum <- 0
    wisk <- 1/mis[h,t]
    for(k in 1:ni[h]){
      for(u in 1:misk[[h]][k,t]){
        sum <- sum + lambda[[h]][k,t]*wisk*F_isP(data, lambda,lambdais, misk, mis, data[[h]][[t]][[k]][u], i,s,w,ni)
      }
    }
    return(sum)
  }
}
# Psi
Psi <- function(data, lambda, lambdais, misk, mis, i,s,h,k,w,ni,a,d){
  if(h != i){
    sum <- 0
    for(t in 1:d){
      if(w == 1) whtk <- 1/(lambdais[h,t]*misk[[h]][k,t])
      if(w == 2) whtk <- 1/(mis[h,t])
      for(u in 1:misk[[h]][k,t]){
        sum <- sum + lambda[[h]][k,t]*whtk*F_isP(data, lambda,lambdais, misk, mis, data[[h]][[t]][[k]][u], i,s,w,ni)
      }
    }
    return(-ni[h]/(a*d)*sum)
  }
  if(h == i){
    sum <- 0
    if(w == 1) wisk <- 1/(lambdais[i,s]*misk[[i]][k,s])
    if(w == 2) wisk <- 1/(mis[i,s])
    for(j in 1:a){
      if(j==i) sum <- sum
      if(j!=i){
        for(t in 1:d){
          for(u in 1:misk[[i]][k,s]){
            sum <- sum + lambda[[i]][k,s] * wisk * F_isP(data, lambda,lambdais, misk, mis, data[[i]][[s]][[k]][u], j,t,w,ni)
          }
        }
      }
    }
    for(t in 1:d){
      for(u in 1:misk[[i]][k,s]){
        sum <- sum + lambda[[i]][k,s]*wisk * F_isP(data, lambda, lambdais, misk, mis, data[[i]][[s]][[k]][u], i,t,w,ni)
      }
    }
    for(t in 1:d){
      if(w == 1) witk <- 1/(lambdais[i,t]*misk[[i]][k,t])
      if(w == 2) witk <- 1/(mis[i,t])
      for(u in 1:misk[[i]][k,t]){
        sum <- sum - lambda[[i]][k,t]*witk*F_isP(data, lambda, lambdais, misk, mis, data[[i]][[t]][[k]][u], i,s,w,ni)
      }
    }
    return(ni[h]/(a*d)*sum)
  }
}
# Beta
Beta <- function(data, lambda, lambdais, misk, mis, i,s,h,k,w,ni,a,d){
  if(h != i){
    sum <- 0
    for(t in 1:d){
      if(w == 1) whtk <- 1/(lambdais[h,t] * misk[[h]][k,t])
      if(w == 2) whtk <- 1/mis[h,t]
      
      sum <- sum + lambda[[h]][k,t]* misk[[h]][k,t]*whtk*p_is_ht(data, lambda,lambdais, i,s,h,t,w,ni,misk,mis,a,d)
      
    }
    return(-ni[h]/(a*d)*sum)
  }
  if(h == i){
    sum <- 0
    if(w == 1) wisk <- 1/(lambdais[i,s]*misk[[i]][k,s])
    if(w == 2) wisk <- 1/mis[i,s]
    for(j in 1:a){
      if(j == i) sum <- sum 
      if(j != i){
        for(t in 1:d){
          sum <- sum + misk[[i]][k,s]*lambda[[i]][k,s]*wisk*p_is_ht(data, lambda,lambdais, j,t,i,s,w,ni,misk,mis,a,d)
        }
      }
    }
    for(t in 1:d){
      sum <- sum + misk[[i]][k,s]*lambda[[i]][k,s]*wisk*p_is_ht(data, lambda,lambdais, i,t,i,s,w,ni,misk,mis,a,d)
    }
    for(t in 1:d){
      if(w == 1) witk <- 1/(lambdais[i,t]*misk[[i]][k,t])
      if(w == 2) witk <- 1/mis[i,t]
      sum <- sum - misk[[i]][k,t]*lambda[[i]][k,t]*witk*p_is_ht(data, lambda,lambdais, i,s,i,t,w,ni,misk,mis,a,d)
    }
    return(ni[h]/(a*d)*sum)
  }
}
#### Function to calculate test  ####
getTest <- function(a,d,ni,N,C,data, lambda,w){
  ### Parameters ###
  # Number of groups
  #a <-length(data)
  # Number of repeated measures
  #d <- unique(lengths(data))
  # ni: Sample size within group i 
  #ni <- unlist(lapply(lapply(data, lengths),unique))
  # N: total sample size
  #N <- sum(ni)
  
  # lambdais: Number of observations in group i at time s
  lambdais <- matrix(unlist(lapply(lambda, colSums)), byrow = TRUE, 
                     ncol = d)
  if(sum(lambdais==0) > 0 | sum(lambdais == 1) > 0){
    res_ATS <- NA
    res_MCTP <- NA
    phat <- rep(NA, a*d)
      }
  else{
  # misk: Number of dependent replicates in group i at time s of subject k
  misk <- replicate(n=a, matrix(nrow = max(ni), ncol = d), simplify = FALSE)
  for(i in 1:a){
    for(s in 1:d){
      for(k in 1:ni[i]){
        misk[[i]][k,s] <- misk_fun(data,i,s,k)
      }
    }
  }
  
  # mis: Number of observed replicates within group i at time s
  mis <- matrix(rep(0,a*d), ncol = d)
  for(i in 1:a){
    for(s in 1:d){
      mis[i,s] <- mis_fun(lambda, misk, ni,i,s)
    }
  }
  # Contrast Matrix
  C <- diag(a*d)-1/(a*d)
  
  # Calculation of Covariance Matrix
  Vhat <- matrix(0,nrow=a*d, ncol = a*d)
  omega <- matrix(0,nrow = a, ncol = a*d)
  for(h in 1:a){
    Psi_mat <- matrix(0,nrow=a*d, ncol = ni[h])
    Beta_mat <- matrix(0,nrow=a*d, ncol = ni[h])
    for(k in 1:ni[h]){
      Psi_help <- rep(0,a*d)
      Beta_help <- rep(0,a*d)
      for(i in 1:a){
        for(s in 1:d){
          Psi_help[(i-1)*d+s] <-  Psi(data, lambda, lambdais, misk, mis, i,s,h,k,w,ni,a,d)
          Beta_help[(i-1)*d+s] <-  Beta(data, lambda, lambdais, misk, mis, i,s,h,k,w,ni,a,d)
        }
      }
      Psi_mat[,k] <- Psi_help
      Beta_mat[,k] <- Beta_help
    }
    Vhat <- Vhat + N/ni[h]*(Psi_mat-Beta_mat)%*%t(Psi_mat-Beta_mat)/(ni[h]-1)
    Phi_mat <- C%*%Psi_mat
    B_mat <- C%*%Beta_mat
    #omega[h,] <- diag((Phi_mat-B_mat)%*%t(Phi_mat-B_mat)/(ni[h]-1))
    omega[h,] <- apply(Phi_mat, 1, var)
  }
  
  # Relative Effect Estimator
  phat <- rep(0,a*d)
  for(i in 1:a){
    for(s in 1:d){
      phat[(i-1)*d+s] <-  p_is(data, lambda, lambdais,i,s,w,ni,misk,mis,a,d)
    }
  }
  
  nc <- nrow(C)
  Cp <- C%*%phat
  CVC <- C%*%Vhat%*%t(C)
  
  # Check if data is generative
  if(sum(Vhat==0) > 0 | sum(is.na(CVC))>0 |isCorrelation(cov2cor(CVC)) == F) {
    res_ATS <- NA
    res_MCTP <- NA
  }
  if(sum(Vhat==0) ==0 & sum(is.na(CVC))==0 & isCorrelation(cov2cor(CVC)) == TRUE) {
    # ATS
    trTV <- sum(diag(CVC))
    trTV2 <- trTV^2
    trTVTV <- sum(diag(CVC%*%CVC))
    dfF <- trTV2/trTVTV
    ATS <- N*t(Cp)%*%Cp/trTV*dfF
    crit_ATS <- qchisq(0.95,dfF)
    res_ATS <- (ATS>crit_ATS)
    
    # MCTP 
    Test <- sqrt(N)*Cp/sqrt(c(diag(CVC)))
    T0 <- max(abs(Test))
    nu <- rep(0, nrow(C))
    for(l in 1:nrow(C)){
      
        sum1 <- sum(omega[,l]/ni)
        sum2 <- sum(omega[,l]^2/(ni^2*(ni-1)))
      
      nu[l] <- sum1^2/sum2
    }
    mnu <- max(1,round(min(nu)))
    pval <- 1- pmvt(lower = -T0,upper = T0, delta=rep(0,nrow(C)), corr = cov2cor(CVC), df = mnu)[1]
    res_MCTP <- (pval < 0.05)
  }}
  return(list(res_ATS, res_MCTP, phat))
}
# Simulation Function
mySimu <- function(n1n2, R, MISK, RHO, SIGMA, w, nsim){
  a <- 2
  d <- 3
  # Sample Size within groups
  # Balanced settings
  if(n1n2 == "15-15") ni <- c(15,15)
  if(n1n2 == "30-30") ni <- c(30,30)
  # # Unbalanced settings
  if(n1n2 == "20-10") ni <- c(20,10)
  if(n1n2 == "40-20") ni <- c(40,20)
  N=sum(ni)
  # Missing values
  if(R == "0.0") r <- c(0,0)
  if(R == "0.1") r <- c(0.1,0.1)
  if(R == "0.3") r <- c(0.3,0.3)
  if(R == "0.0-0.2") r <- c(0.0,0.2)
  
  
  T1ATS <- rep(0,nsim)
  T1MCTP <- rep(0,nsim)
  phat <- matrix(rep(0, nsim*a*d), ncol = nsim)
  
  for(isim in 1:nsim){
    # Number of dependent replicates of subject k at time s in group i
    if(MISK == "1") misk_sim <- rep(1,d)
    if(MISK == "2") misk_sim <- rep(2,d)
    if(MISK == "BINOM-5-0.6+1"){
      misk_sim <- replicate(n=a, matrix(nrow = max(ni), ncol = d), simplify = FALSE)
      for(i in 1:a){
        for(s in 1:d){
          for(k in 1:ni[i]){
            misk_sim[[i]][k,s] <- rbinom(1,5,0.6)+1
          }
        }
      }
    }
    if(MISK == "BINOM-10-0.4+1"){
      misk_sim <- replicate(n=a, matrix(nrow = max(ni), ncol = d), simplify = FALSE)
      for(i in 1:a){
        for(s in 1:d){
          for(k in 1:ni[i]){
            misk_sim[[i]][k,s] <- rbinom(1,10,0.4)+1
          }
        }
      }
    }
    # Correlation between dependent replicates
    if(RHO == "0.0") rho <- 0
    if(RHO == "0.3") rho <- 0.3
    if(RHO == "0.9") rho <- 0.9
    if(RHO == "BINOM-10-0.6/10"){
      rho <- replicate(n=a, matrix(nrow = max(ni), ncol = d), simplify = FALSE)
      for(i in 1:a){
        for(s in 1:d){
          for(k in 1:ni[i]){
            rho[[i]][k,s] <- rbinom(1,10,0.6)/10
          }
        }
      }
    }
    # Covariance matrices
    if(SIGMA == "Homo"){
      Sigma <- matrix(c(1,0.2, 0.2, 0.2, 1, 0.2, 0.2, 0.2, 1), byrow = T, ncol = d)
    }
    if(SIGMA == "Het"){
      Sigma <- matrix(c(1, 0.1, 0.2, 0.1, 1.2, 0.3, 0.2, 0.3, 1.5), byrow = T, ncol = d)
    }
    
    ## Initialize data strucutre ##
    # Initialize list for group 1
    dat1 <- list(list(c(0)),
                 list(c(0)),
                 list(c(0)))
    for(i in 1:c(ni[1]-1)){
      dat1[[1]] <- append(dat1[[1]],0)
      dat1[[2]] <- append(dat1[[2]],0)
      dat1[[3]] <- append(dat1[[3]],0)
    }
    # Initialize list for group 2
    dat2 <- list(list(c(0)),
                 list(c(0)),
                 list(c(0)))
    for(i in 1:c(ni[2]-1)){
      dat2[[1]] <- append(dat2[[1]],0)
      dat2[[2]] <- append(dat2[[2]],0)
      dat2[[3]] <- append(dat2[[3]],0)
    }
    if((MISK == "1")){
      V <- Sigma
      VWurzel <- MWurzel(V)
      x <-t(VWurzel%*% matrix(rnorm(ni[1]*sum(misk_sim)),sum(misk_sim),ni[1]))
      y <-t(VWurzel%*% matrix(rnorm(ni[2]*sum(misk_sim)),sum(misk_sim),ni[2]))
      
      # Fill list for group 1
      for(k in 1:ni[1]){
        help <- 1             
        for(s in 1:d){
          for(u in 1:misk_sim[s]){
            dat1[[s]][[k]][u] <- x[k,help]
            help <- help + 1
          }
        }
      }
      
      # Fill list for group 2
      for(k in 1:ni[2]){
        help <- 1             
        for(s in 1:d){
          for(u in 1:misk_sim[s]){
            dat2[[s]][[k]][u] <- y[k,help]
            help <- help + 1
          }
        }
      }
      # Combine both lists 
      data <- list(dat1,dat2)
    }
    if((MISK =="2")){
      if(RHO == "0.0" | RHO == "0.3" | RHO == "0.9"){
        cmisk <- c(0,cumsum(misk_sim))
        # Correlation in cluster
        rho <- rep(rho, sum(misk_sim)^2)
        V <- matrix(rho, ncol = sum(misk_sim))
        # Blow up covariance matrix
        for(i in 1:d){
          for(j in 1:d){
            if(i != j){
              for(k in 1:misk_sim[j]){
                for(l in 1:misk_sim[i]){
                  V[cmisk[i]+l,cmisk[j]+k] <- Sigma[i,j]
                }
              }
            }
          }
        }
        diag(V) <- rep(diag(Sigma), misk_sim)
        VWurzel <- MWurzel(V)
        x <-t(VWurzel%*% matrix(rnorm(ni[1]*sum(misk_sim)),sum(misk_sim),ni[1]))
        y <-t(VWurzel%*% matrix(rnorm(ni[2]*sum(misk_sim)),sum(misk_sim),ni[2]))
        # Fill list for group 1
        for(k in 1:ni[1]){
          help <- 1             
          for(s in 1:d){
            for(u in 1:misk_sim[s]){
              dat1[[s]][[k]][u] <- x[k,help]
              help <- help + 1
            }
          }
        }
        # Fill list for group 2
        for(k in 1:ni[2]){
          help <- 1             
          for(s in 1:d){
            for(u in 1:misk_sim[s]){
              dat2[[s]][[k]][u] <- y[k,help]
              help <- help + 1
            }
          }
        }
        # Combine both lists 
        data <- list(dat1,dat2)
      }
      if(RHO == "BINOM-10-0.6/10"){
        # Combine both lists
        data <- list(dat1,dat2)
        for(ii in 1:a){
          for(kk in 1:ni[ii]){
            rho_k <- rho[[ii]][kk,]
            #rho_k <- rep(rho_k, sum(misk_sim)^2)
            V <- as.matrix(bdiag(matrix(rep(rho_k[1], misk_sim[1]^2), ncol = misk_sim[1]),
                                 matrix(rep(rho_k[2], misk_sim[2]^2), ncol = misk_sim[2]),
                                 matrix(rep(rho_k[3], misk_sim[3]^2), ncol = misk_sim[3])))
            
            cmisk <- c(0,cumsum(misk_sim))
            # Blow up covariance matrix
            for(i in 1:d){
              for(j in 1:d){
                if(i != j){
                  for(k in 1:misk_sim[j]){
                    for(l in 1:misk_sim[i]){
                      V[cmisk[i]+l,cmisk[j]+k] <- Sigma[i,j]
                    }
                  }
                }
              }
            }
            diag(V) <- rep(diag(Sigma), misk_sim)
            VWurzel <- MWurzel(V)
            x <-t(VWurzel%*% matrix(rnorm(1*sum(misk_sim)),sum(misk_sim),1))
            #y <-t(VWurzel%*% matrix(rnorm(1*sum(misk_sim)),sum(misk_sim),1))
            data[[ii]][[1]][[kk]] <- x[1:2]
            data[[ii]][[2]][[kk]] <- x[3:4]
            data[[ii]][[3]][[kk]] <- x[5:6]
          }
        }
      }
    }
    if(MISK == "BINOM-5-0.6+1" | MISK == "BINOM-10-0.4+1"){
      if(RHO == "0.0" | RHO == "0.3" | RHO == "0.9"){
        # Obtain maximum of dependent replicates
        misk_1k_max <- apply(misk_sim[[1]], 2, maxna)
        misk_2k_max <- apply(misk_sim[[2]], 2, maxna)
        
        misk_k_mat <- rbind(misk_1k_max, misk_2k_max)
        misk_k <- apply(misk_k_mat, 2, max)
        # Cummulative number of dependent replicates
        cmisk <- c(0,cumsum(misk_k))
        # Correlation in cluster
        rho <- rep(rho, sum(misk_k)^2)
        V <- matrix(rho, ncol = sum(misk_k))
        # Blow up covariance matrix
        for(i in 1:d){
          for(j in 1:d){
            if(i != j){
              for(k in 1:misk_k[j]){
                for(l in 1:misk_k[i]){
                  V[cmisk[i]+l,cmisk[j]+k] <- Sigma[i,j]
                }
              }
            }
          }
        }
        diag(V) <- rep(diag(Sigma), misk_k)
        VWurzel <- MWurzel(V)
        x <-t(VWurzel%*% matrix(rnorm(ni[1]*sum(misk_k)),sum(misk_k),ni[1]))
        y <-t(VWurzel%*% matrix(rnorm(ni[2]*sum(misk_k)),sum(misk_k),ni[2]))
        
        # Fill list for group 1
        for(k in 1:ni[1]){
          help <- 1             
          for(s in 1:d){
            for(u in 1:misk_k[s]){
              dat1[[s]][[k]][u] <- x[k,help]
              help <- help + 1
            }
          }
        }
        
        # Fill list for group 2
        for(k in 1:ni[2]){
          help <- 1             
          for(s in 1:d){
            for(u in 1:misk_k[s]){
              dat2[[s]][[k]][u] <- y[k,help]
              help <- help + 1
            }
          }
        }
        # Combine both lists 
        data <- list(dat1,dat2)
        # Delete observations according to misk_sim
        for(i in 1:a){
          for(s in 1:d){
            for(k in 1:ni[i]){
              data[[i]][[s]][[k]] <- data[[i]][[s]][[k]][1:misk_sim[[i]][k,s]]
            }
          }
        }
      }
      if(RHO == "BINOM-10-0.6/10"){
        # Obtain maximum of dependent replicates
        misk_1k_max <- apply(misk_sim[[1]], 2, maxna)
        misk_2k_max <- apply(misk_sim[[2]], 2, maxna)
        
        misk_k_mat <- rbind(misk_1k_max, misk_2k_max)
        misk_k <- apply(misk_k_mat, 2, max)
        # Cumulative number of dependent replicates
        cmisk <- c(0,cumsum(misk_k))
        # Combine both lists
        data <- list(dat1,dat2)
        for(ii in 1:a){
          for(kk in 1:ni[ii]){
            rho_k <- rho[[ii]][kk,]
            #rho_k <- rep(rho_k, sum(misk_sim)^2)
            V <- as.matrix(bdiag(matrix(rep(rho_k[1], misk_k[1]^2), ncol = misk_k[1]),
                                 matrix(rep(rho_k[2], misk_k[2]^2), ncol = misk_k[2]),
                                 matrix(rep(rho_k[3], misk_k[3]^2), ncol = misk_k[3])))
            
            # Blow up covariance matrix
            for(i in 1:d){
              for(j in 1:d){
                if(i != j){
                  for(k in 1:misk_k[j]){
                    for(l in 1:misk_k[i]){
                      V[cmisk[i]+l,cmisk[j]+k] <- Sigma[i,j]
                    }
                  }
                }
              }
            }
            diag(V) <- rep(diag(Sigma), misk_k)
            VWurzel <- MWurzel(V)
            x <-t(VWurzel%*% matrix(rnorm(1*sum(misk_k)),sum(misk_k),1))
            data[[ii]][[1]][[kk]] <- x[1:misk_k[1]]
            data[[ii]][[2]][[kk]] <- x[misk_k[1]+1:misk_k[2]]
            data[[ii]][[3]][[kk]] <- x[misk_k[2]+1:misk_k[3]]
            # Delete observations according to misk_sim
            for(i in 1:a){
              for(s in 1:d){
                for(k in 1:ni[i]){
                  data[[i]][[s]][[k]] <- data[[i]][[s]][[k]][1:misk_sim[[i]][k,s]]
                }
              }
            }
          }
        }
      }
    }
    # Missing Values
    lambda <- list(
      matrix(rbinom(ni[1]*d,1,1-r[1]), ncol = d),
      matrix(rbinom(ni[2]*d,1,1-r[2]), ncol = d)
    )
    test <- getTest(a,d,ni,N,C,data, lambda, w)
    T1ATS[isim] <- test[[1]]
    T1MCTP[isim] <- test[[2]]
    phat[,isim] <- test[[3]]-0.5
  }
  # Bias and MSE
  MSE <- sum(apply((phat)^2,1, meanna))/(a*d)
  Bias <- sum(apply(phat ,1, meanna))/(a*d)
  
  
  
  result <- data.frame( n1n2 = n1n2, r = R, MISK = MISK, RHO = RHO, S = SIGMA, w = w, 
                       T1ATS = mean(T1ATS, na.rm = TRUE), T1MCTP = mean(T1MCTP, na.rm = TRUE), MSE, Bias)
  
  write.table(result,paste("results",Index,".txt",sep=""), eol="\r\n", row.names = F, col.names = F, quote = F, append = T)
}
# Parameters for Simulation 
n1n2 <- c("15-15", "30-30", "20-10", "40-20")
R <- c("0.0", "0.1", "0.3","0.0-0.2")
MISK <- c("1", "2", "BINOM-5-0.6+1", "BINOM-10-0.4+1")
RHO <- c("0.0", "0.3", "0.9","BINOM-10-0.6/10")
SIGMA <- c("Homo", "Het")
W <- c(1,2)
for(h in 1:length(n1n2)){
  for(hh in 1:length(R)){
   for(hhh in 1:1){
    for(hhhh in 1:1){
        for(hhhhh in 1:length(SIGMA)){
          for(hhhhhh in 1:length(W)){
            mySimu(n1n2 = n1n2[h], R = R[hh], MISK = MISK[hhh],
                   RHO = RHO[hhhh], SIGMA = SIGMA[hhhhh],
                   w = W[hhhhhh],nsim = Durchlaeufe)
          }
        }
      }
    }
  }
}


for(h in 1:length(n1n2)){
  for(hh in 1:length(R)){
    for(hhh in 2:length(MISK)){
      for(hhhh in 1:length(RHO)){
        for(hhhhh in 1:length(SIGMA)){
          for(hhhhhh in 1:length(W)){
            mySimu(n1n2 = n1n2[h], R = R[hh], MISK = MISK[hhh],
                   RHO = RHO[hhhh], SIGMA = SIGMA[hhhhh],
                   w = W[hhhhhh],nsim = Durchlaeufe)
          }
        }
      }
    }
  }
}
B=Sys.time()
B-A
