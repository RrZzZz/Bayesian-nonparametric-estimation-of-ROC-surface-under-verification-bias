DP_ROC <- function(niter, S, L.org, lambda){
  # niter: number of iterations for MCMC run
  # S: observations 
  # L.org: labels corresponding to S. 1,2,3 represents different classes while 4 means unknown.
  # lambda: disease prevalence rates
  
  # assign initial values to labels using naive bayes
  L=L.org
  
  den1 <- density(S[L==1])
  den2 <- density(S[L==2])
  den3 <- density(S[L==3])
  
  for(i in 1:length(L)){
    if(L[i]==4){
      L[i] <- which.max(c(den1$y[which.min(abs(den1$x-S[i]))], den2$y[which.min(abs(den2$x-S[i]))], den1$y[which.min(abs(den2$x-S[i]))]))
    }
  }
  
  # initial values for parameters
  mu <- S
  sigma <- sqrt(1/rgamma(length(S), 0.1, 0.1))
  n <- c(sum(L==1), sum(L==2), sum(L==3))
  m <- c(mean(S[L==1]), mean(S[L==2]), mean(S[L==3]))
  M <- c(1, 1, 1)
  a <- c(1, 1, 1)
  s <- c(0.1, 0.1, 0.1)
  beta <- c(0.1, 0.1, 0.1)
  
  # grid to store function values
  grid = seq(min(S)-sd(S), max(S)+sd(S), length.out=200)
  val <- list()
  
  # store function values on the grid points
  F1 <- matrix(0, niter, length(grid))
  F2 <- matrix(0, niter, length(grid))
  F3 <- matrix(0, niter, length(grid))
  
  f.un <- list()
  
  for(iter in 1:niter){
    # assign labels to unverified points
    for(i in 1:length(S)){
      if(L.org[i]==4){
        mu.new <- mu[-i]
        L.new <- L[-i]
        sigma.new <- sigma[-i]
        n.new <- c(sum(L.new==1), sum(L.new==2), sum(L.new==3))
        fval <- c()
        for(k in 1:3){
          v=rbeta(1, M[k], n.new[k])
          N=100
          if(runif(1)<v){
            q.gen <- rdirichlet(1, rep(M[k]/N, N)) 
            sigma.gen <- sqrt(1/rgamma(N, s[k], beta[k]))
            mu.gen <- c()
            for(i in 1:length(q.gen)){
              mu.gen[i] <- rnorm(1, m[k], sigma.gen[i]/a[k])
            }
            f <- function(t){
              sum(dnorm(t, mu.gen, sigma.gen)*q.gen)
            }
          }else{
            w.gen <- rdirichlet(1, rep(1, n.new[k]))
            f <- function(t){
              sum(dnorm(t,mu.new[which(L.new==k)], sigma.new[which(L.new==k)])*w.gen)
            }
          }
          fval[k] <- f(S[i])
          
        }
        L[i] <- which(as.vector(rmultinom(1, 1, lambda*fval) )!=0)
        
      }
      n <- c(sum(L==1), sum(L==2), sum(L==3))
      v <- L[i]
      q <- c()
      
      # update values of mu[i] and sigma[i]
      for(j in 1:length(S)){
        if(L[j]==v){
          q[j]=1/sigma[j]*exp(-(S[i]-mu[j])^2/(2*sigma[j]^2))
        }else{
          q[j]=0
        }
      }
      q[i]=0
      q[length(S)+1] <- (M[v]*sqrt(a[v])*gamma(s[v]+0.5)*beta[v]^s[v])/(sqrt(1+a[v])*gamma(s[v])*(beta[v]+a[v]*(S[i]-m[v])^2/(2*(1+a[v])))^(s[v]+0.5))
      
      q <- q/sum(q)
      r <- which(rmultinom(1,1,prob=q)==1)
      if(r==length(S)+1){
        mu[i] <- rnorm(1, (S[i]+a[v]*m[v])/(1+a[v]), sigma[i]/sqrt((1+a[v])) )
        sigma[i] <- sqrt(1/(rgamma(1, s[v]+0.5, beta[v]+a[v]*(S[i]-m[v])^2/(2*(1+a[v])))))
      }else{
        mu[i] <- mu[r]
        sigma[i] <- sigma[r]
      }
      
      
    }
    
    # sample the distributions for three categories
    for(k in 1:3){
      
      v=rbeta(1, M[k], n[k])
      N=100
      if(runif(1)<v){
        q.gen <- rdirichlet(1, rep(M[k]/N, N)) #Let M=10*n
        sigma.gen <- sqrt(1/rgamma(N, s[k], beta[k]))
        mu.gen <- c()
        for(i in 1:length(q.gen)){
          mu.gen[i] <- rnorm(1, m[k], sigma.gen/a[k])
        }
        F <- function(t){
          sum(pnorm(t, mu.gen, sigma.gen)*q.gen)
        }
      }else{
        w.gen <- rdirichlet(1, rep(1, n[k]))
        F <- function(t){
          sum(pnorm(t,mu[which(L==k)], sigma[which(L==k)])*w.gen)
        }
      }
      val[[k]] <- unlist(lapply(grid, F))
      
    }
    
    
    F1[iter, ] <- val[[1]]
    F2[iter, ] <- val[[2]]
    F3[iter, ] <- val[[3]]
  }
  
  
  # take the means for every points of the distributions
  F1.mean <- colMeans(F1[(niter*0.05):niter, ])
  F2.mean <- colMeans(F2[(niter*0.05):niter, ])
  F3.mean <- colMeans(F3[(niter*0.05):niter, ])
  
  F1.mean <- c(0, F1.mean, 1)
  F2.mean <- c(0, F2.mean, 1)
  F3.mean <- c(0, F3.mean, 1)
  
  ROC <- function(c1, c2){
    tcf2 <- F2.mean[c2]-F2.mean[c1]
  }
  
  tcf1 <- F1.mean
  tcf3 <- rev(1-F3.mean)
  tcf2 <- outer(c(1:202), c(202:1), FUN = "ROC")
  
  tcf2[tcf2<0]=0
  
  # plot ROC surface
  zz <- tcf2*20
  
  library(rgl)
  zlim        <- range(zz,na.rm=T)
  zlen        <- zlim[2] - zlim[1] + 1
  color.range <- rev(rainbow(zlen))       # height color lookup table
  colors      <- color.range[zz-zlim[1]+1] # assign colors to heights for each point
  persp3d(main="true ROC surface",tcf1, tcf3, tcf2, col=colors)
  
  surface3d(tcf1, tcf3, tcf2, back = "lines")
  surface3d(tcf1, tcf3, tcf2, front = "lines")

  # calculate the VUS
  grid.length1 <- tcf1[-1]-tcf1[-length(tcf1)]
  grid.length3 <- tcf3[-1]-tcf3[-length(tcf3)]
  grid.length1 <- as.matrix(grid.length1)
  grid.length3 <- as.matrix(grid.length3)
  
  gridsq <- grid.length1%*%t(grid.length3)
  gridval <- 1/4*(tcf2[1:201, 1:201]+tcf2[1:201, 2:202]+tcf2[2:202, 1:201]+tcf2[2:202, 2:202])
  vus <- sum(gridsq*gridval)

  return(list(vus=vus, F1=F1.mean, F2=F2.mean, F3=F3.mean))
}
