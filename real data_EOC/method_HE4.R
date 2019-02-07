library(MCMCpack)

# read data in
setwd("//wolftech.ad.ncsu.edu/cos/stat/Redirect/rzhu4/Downloads/research 2/real data")
dat <- read.csv("data.csv")[-1,]

x <- dat$Early

L <- c()
L[which(x==0)]=2
L[which(x==1)]=3

y <- dat$X.1

L[which(y=="GP Control")]=1
L[which(y=="BD Control")]=1
L[which(y=="SS Control")]=1


CA125 <- dat$Harvard.1
HE4 <- dat$FHCRC

# delete missing values
index <- c(which(CA125==""), which(HE4==""))
index <- unique(index)
CA125 <- CA125[-index]
HE4 <- HE4[-index]
L <- L[-index]

CA125 <- as.numeric(paste(CA125))
HE4 <- as.numeric(paste(HE4))

# Exploratory analysis
L[order(CA125)]
L[order(HE4)]

data <- cbind(L, CA125, HE4)
data <- data[-which(CA125>1000),]

L <- data[,1]
CA125 <- log(data[,2])
HE4 <- log(data[,3])
lambda <- c(10,1,1)/12

S <- HE4
# S <- CA125

L[which(is.na(L))] <- 4
L.org=L

# initial values 

den1 <- density(S[L==1])
den2 <- density(S[L==2])
den3 <- density(S[L==3])

for(i in 1:length(L)){
  if(L[i]==4){
    L[i] <- sample(c(1:3),size=1, prob=c(den1$y[which.min(abs(den1$x-S[i]))], den2$y[which.min(abs(den2$x-S[i]))], den1$y[which.min(abs(den2$x-S[i]))])*lambda)
  }
}

mu <- S
sigma <- sqrt(1/rgamma(length(S), 0.1, 0.1))
#initial value for sigma???

n <- c(sum(L==1), sum(L==2), sum(L==3))

m <- c(mean(L[which(L==1)]), mean(L[which(L==2)]), mean(L[which(L==3)]))
M <- c(1, 1, 1)
a <- c(1, 1, 1)
s <- c(0.1, 0.1, 0.1)
beta <- c(0.1, 0.1, 0.1)

niter=100000

grid = seq(min(S)-sd(S), max(S)+sd(S), length.out=200)

val <- list()

F1 <- matrix(0, niter, length(grid))
F2 <- matrix(0, niter, length(grid))
F3 <- matrix(0, niter, length(grid))

f.un <- list()

# MCMC
for(iter in 1:niter){
  
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


# F1.mean <- colMeans(F1[500:(iter-1), ])
# F2.mean <- colMeans(F2[500:(iter-1), ])
# F3.mean <- colMeans(F3[500:(iter-1), ])

# distributions
F1.mean <- colMeans(F1[5000:niter, ])
F2.mean <- colMeans(F2[5000:niter, ])
F3.mean <- colMeans(F3[5000:niter, ])


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

zz <- tcf2*20

library(rgl)
zlim        <- range(zz,na.rm=T)
zlen        <- zlim[2] - zlim[1] + 1
color.range <- rev(rainbow(zlen))       # height color lookup table
colors      <- color.range[zz-zlim[1]+1] # assign colors to heights for each point
persp3d(main="trinormal estimator",tcf1, tcf3, tcf2, col=colors)

surface3d(tcf1, tcf3, tcf2, back = "lines")
surface3d(tcf1, tcf3, tcf2, front = "lines")

###############################################

# k=1
#   xgrid= grid
#   plot(xgrid, dnorm(xgrid, mu0[k], sigma0[k]), type="l", ylim=c(0, 0.3), lwd=3)
# 
#   for(i in 1:niter){
#     lines(xgrid, result[[i]][[k]], col=i)
#   }
#   # lines(xgrid, dnorm(xgrid, mu0[k], sigma0[k]), type="l", lwd=3)
# 
# 
#   plot(xgrid, dnorm(xgrid, mu0[k], sigma0[k]), type="l", ylim=c(0, 0.3), lwd=3)
# 
#   for(i in 1:200){
#     lines(xgrid, result[[i]][[k]], col=i)
#   }
# 
# 
# hist(z)


###############################################

grid.length1 <- tcf1[-1]-tcf1[-length(tcf1)]
grid.length3 <- tcf3[-1]-tcf3[-length(tcf3)]
grid.length1 <- as.matrix(grid.length1)
grid.length3 <- as.matrix(grid.length3)

gridsq <- grid.length1%*%t(grid.length3)



gridval <- 1/4*(tcf2[1:201, 1:201]+tcf2[1:201, 2:202]+tcf2[2:202, 1:201]+tcf2[2:202, 2:202])
vus <- sum(gridsq*gridval)

save.image(file="result_HE4.Rdata")


