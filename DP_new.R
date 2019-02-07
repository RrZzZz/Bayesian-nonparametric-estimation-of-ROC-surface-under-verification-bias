library(MCMCpack)
library(mvtnorm)
library(dplyr)
# read data in
data <- read.csv("//wolftech.ad.ncsu.edu/cos/stat/Redirect/rzhu4/Downloads/hcc-dataset/hcc-data.csv", header = F)
L.org <- read.table("//wolftech.ad.ncsu.edu/cos/stat/Redirect/rzhu4/Downloads/hcc-dataset/hcc-stages.txt")
L.org <- L.org$V1
L.org[1]=2 # For some reason the first label is not read correctly so we manually adjust the label
table(L.org)

i=36 # Serum albumin
data[which(data[,i]=="?"),i]=NA
data[,i]=as.numeric(as.character(data[,i]))
summary(data[,i])
L.org[order(data[,i])]
L=case_when(
  L.org==0 | L.org==1 ~0,
  L.org==2 | L.org==3 ~1,
  L.org==4 ~2
) # combine some levels together to form 3 classes

# Exploratory analysis
L[order(data[,i])]
mean(data[which(L==0), i], na.rm = T)
mean(data[which(L==1), i], na.rm = T)
mean(data[which(L==2), i], na.rm = T)
S=data[,i]

gender=data[,1]

L[order(S)]
L=L[-c(which(is.na(S)))]
gender=gender[-c(which(is.na(S)))]
S=S[-c(which(is.na(S)))]



# Since the missing proportion is very small, we will use the sample proportion as an estimate of the population disease prevalence rates
lambda <- as.numeric(table(L)/sum(table(L)))

L=L+1
L[which(is.na(L))] <- 4
L.org=L

den1 <- density(S[L==1])
den2 <- density(S[L==2])
den3 <- density(S[L==3])

#initial values
for(i in 1:length(L)){
  if(L[i]==4){
    L[i] <- sample(c(1:3),size=1, prob=c(den1$y[which.min(abs(den1$x-S[i]))], den2$y[which.min(abs(den2$x-S[i]))], den1$y[which.min(abs(den2$x-S[i]))])*lambda)
  }
}

mu <- rbind(S, rep(0, length(S)))
sigma <- sqrt(1/rgamma(length(S), 0.1, 0.1))

n <- c(sum(L==1), sum(L==2), sum(L==3))

m <- rbind(c(mean(S[which(L==1&gender==0)]), mean(S[which(L==2&gender==0)]), mean(S[which(L==3&gender==0)])), c(mean(S[which(L==1&gender==1)]), mean(S[which(L==2&gender==1)]), mean(S[which(L==3&gender==1)])))
m[2,]=m[2,]-m[1,]
M <- c(1, 1, 1)
a=list()
a[[1]]=diag(c(0.1,0.1))
a[[2]]=diag(c(0.1,0.1))
a[[3]]=diag(c(0.1,0.1))
s <- c(0.1, 0.1, 0.1)
beta <- c(0.1, 0.1, 0.1)

niter=100000

grid = seq(min(S)-sd(S), max(S)+sd(S), length.out=200)

val0 <- list()
val1 <- list()

# Store posterior distributions
F10 <- matrix(0, niter, length(grid))
F20 <- matrix(0, niter, length(grid))
F30 <- matrix(0, niter, length(grid))
F11 <- matrix(0, niter, length(grid))
F21 <- matrix(0, niter, length(grid))
F31 <- matrix(0, niter, length(grid))

f.un <- list()

ptm <- proc.time()

# MCMC iterations
for(iter in 1:niter){
  
  for(i in 1:length(S)){
    if(L.org[i]==4){
      mu.new <- mu[,-i]
      L.new <- L[-i]
      sigma.new <- sigma[-i]
      n.new <- c(sum(L.new==1), sum(L.new==2), sum(L.new==3))
      fval <- c()
      for(k in 1:3){
        v=rbeta(1, M[k], n.new[k])
        N=100
        if(runif(1)<v){
          q.gen <- rdirichlet(1, rep(M[k]/N, N)) 
          sigmasq.gen <- 1/rgamma(N, s[k], beta[k])
          mu.gen <- matrix(0, 2, length(q.gen))
          for(i in 1:length(q.gen)){
            mu.gen[,i] <- rmvnorm(1, m[,k], sigmasq.gen[i]*solve(a[[k]]))
          }
          f <- function(t, z){
            sum(dnorm(t, matrix(c(1,z),1,2)%*%mu.gen, sigmasq.gen)*q.gen)
          }
        }else{
          w.gen <- rdirichlet(1, rep(1, n.new[k]))
          f <- function(t, z){
            sum(dnorm(t, matrix(c(1,z),1,2)%*%mu.new[,which(L.new==k)], sigma.new[which(L.new==k)])*w.gen)
          }
        }
        fval[k] <- f(S[i], gender[i])
        
      }
      L[i] <- which(as.vector(rmultinom(1, 1, lambda*fval) )!=0)
      
    }
    n <- c(sum(L==1), sum(L==2), sum(L==3))
    v <- L[i]
    q <- c()
    
    for(j in 1:length(S)){
      if(L[j]==v){
        q[j]=1/sigma[j]*exp(-(S[i]-mu[1,j]-mu[2,j]*gender[i])^2/(2*sigma[j]^2))
      }else{
        q[j]=0
      }
    }
    q[i]=0
    if(gender[i]==0){
      sigma.star=a[[v]]+matrix(c(1,0,0,0),2,2)
      mu.star=solve(sigma.star)%*%(a[[v]]%*%m[,v]+S[i]*c(1, 0))
    }else{
      sigma.star=a[[v]]+matrix(c(1,1,1,1),2,2)
      mu.star=solve(sigma.star)%*%(a[[v]]%*%m[,v]+S[i]*c(1, 1))
    }
  
    q[length(S)+1] <- M[[v]]*beta[v]*gamma(s[v]+0.5)*sqrt(det(a[[v]]))/((2*pi)^0.5*gamma(s[v])*sqrt(det(sigma.star)) )*
      (beta[v]+0.5*(S[i]^2+t(m[,v])%*%a[[v]]%*%m[,v]-t(mu.star)%*%sigma.star%*%mu.star)^-(s[v]+0.5) )
    q <- q/sum(q)
    r <- which(rmultinom(1,1,prob=q)==1)
    if(r==length(S)+1){
      mu[,i] <- rmvnorm(1, mu.star, solve(sigma.star)*sigma[i]^2 )
      sigma[i] <- sqrt(1/(rgamma(1, s[v]+0.5, beta[v]+0.5*(S[i]^2+t(m[,v])%*%a[[v]]%*%m[,v]-t(mu.star)%*%sigma.star%*%mu.star)  ) ))
    }else{
      mu[,i] <- mu[,r]
      sigma[i] <- sigma[r]
    }
    
    
  }
  
  
  for(k in 1:3){
    
    v=rbeta(1, M[k], n[k])
    N=100
    if(runif(1)<v){
      q.gen <- rdirichlet(1, rep(M[k]/N, N)) 
      sigmasq.gen <- 1/rgamma(N, s[k], beta[k])
      mu.gen <- matrix(0, 2, length(q.gen))
      for(i in 1:length(q.gen)){
        mu.gen[,i] <- rmvnorm(1, m[,k], sigmasq.gen[i]*solve(a[[k]]))
      }
      F <- function(t, z){
        sum(pnorm(t, matrix(c(1,z),1,2)%*%mu.gen, sigmasq.gen)*q.gen)
      }
      
    }else{
      w.gen <- rdirichlet(1, rep(1, n[k]))
      F <- function(t, z){
        sum(pnorm(t, matrix(c(1,z),1,2)%*%mu[,which(L==k)], sigma[which(L==k)])*w.gen)
      }
    }
    F.z0 <- function(t){
      F(t,0)
    }
    F.z1 <- function(t){
      F(t,1)
    }
    
    
    val0[[k]] <- unlist(lapply(grid, F.z0))
    val1[[k]] <- unlist(lapply(grid, F.z1))
  }
  
  
  F10[iter, ] <- val0[[1]]
  F20[iter, ] <- val0[[2]]
  F30[iter, ] <- val0[[3]]
  
  F11[iter, ] <- val1[[1]]
  F21[iter, ] <- val1[[2]]
  F31[iter, ] <- val1[[3]]
}

# Stop the clock
proc.time() - ptm

# > proc.time() - ptm
# user  system elapsed 
# 7068.75    4.64 7082.61  2 hours approximately

# distributions for gender=0
F10.mean <- colMeans(F10[20001:niter, ])
F20.mean <- colMeans(F20[20001:niter, ])
F30.mean <- colMeans(F30[20001:niter, ])

# distributions for gender=1
F11.mean <- colMeans(F11[20001:niter, ])
F21.mean <- colMeans(F21[20001:niter, ])
F31.mean <- colMeans(F31[20001:niter, ])

# ratio=sum(gender==1)/length(gender)
# 
# distributions overall
# F1.mean <- ratio*F11.mean+(1-ratio)*F10.mean
# F2.mean <- ratio*F21.mean+(1-ratio)*F20.mean
# F3.mean <- ratio*F31.mean+(1-ratio)*F30.mean
# 
# F1.mean <- c(0, F1.mean, 1)
# F2.mean <- c(0, F2.mean, 1)
# F3.mean <- c(0, F3.mean, 1)

F1.mean <- c(0, F10.mean, 1)
F2.mean <- c(0, F20.mean, 1)
F3.mean <- c(0, F30.mean, 1)
# 
# F1.mean <- c(0, F11.mean, 1)
# F2.mean <- c(0, F21.mean, 1)
# F3.mean <- c(0, F31.mean, 1)

ROC <- function(c1, c2){
  tcf2 <- F2.mean[c2]-F2.mean[c1]
}

tcf1 <- F3.mean
tcf3 <- rev(1-F1.mean)
tcf2 <- outer(c(1:202), c(202:1), FUN = "ROC")

tcf2[tcf2<0]=0

zz <- tcf2*20

library(rgl)
zlim        <- range(zz,na.rm=T)
zlen        <- zlim[2] - zlim[1] + 1
color.range <- rev(rainbow(zlen))       # height color lookup table
colors      <- color.range[zz-zlim[1]+1] # assign colors to heights for each point
persp3d(main="DP estimator",tcf1, tcf3, tcf2, col=colors)

surface3d(tcf1, tcf3, tcf2, back = "lines")
surface3d(tcf1, tcf3, tcf2, front = "lines")


grid.length1 <- tcf1[-1]-tcf1[-length(tcf1)]
grid.length3 <- tcf3[-1]-tcf3[-length(tcf3)]
grid.length1 <- as.matrix(grid.length1)
grid.length3 <- as.matrix(grid.length3)

gridsq <- grid.length1%*%t(grid.length3)



gridval <- 1/4*(tcf2[1:201, 1:201]+tcf2[1:201, 2:202]+tcf2[2:202, 1:201]+tcf2[2:202, 2:202])
vus <- sum(gridsq*gridval)

setwd("//wolftech.ad.ncsu.edu/cos/stat/Redirect/rzhu4/Downloads/research 2/new real data")
save.image(file="gender0.Rdata")


