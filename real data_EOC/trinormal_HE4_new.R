library(matrixStats)
library(truncnorm)
library(MCMCpack)
library(bcROCsurface)

# functions to calculate ROC surface
genabcd <- function(truevalue){
  a=1/truevalue[2]
  b=(truevalue[1])/truevalue[2]
  c=1/truevalue[4]
  d=(truevalue[3])/truevalue[4]
  c(a,b,c,d)
}

VUS2 <- function(val){
  a <- val[1]
  b <- val[2]
  c <- val[3]
  d <- val[4]
  f <- function(s){
    pnorm(a*s-b)*pnorm(-c*s+d)*dnorm(s)
  }
  integrate(f, -Inf, Inf)$value
}

# trinormality method without verification bias
logrank <- function(L, n, intval, iter.num=50000){
  n0 <- n[1]
  n1 <- n[2]
  n2 <- n[3]
  mu1 <- intval[1]
  mu2 <- intval[3]
  sigma1 <- intval[2]
  sigma2 <- intval[4]
  N <- length(L)
  sample <- matrix(0, iter.num, 4)
  sample1 <- matrix(0, iter.num, 4)
  
  for(i in 1:iter.num){
    y0 <- rnorm(n0,0,1)
    y1 <- rnorm(n1, mu1, sigma1)
    y2 <- rnorm(n2, mu2, sigma2)
    Q <- c(y0,y1,y2)
    Q.s <- sort(Q)
    Q.s <- c(-Inf, Q.s, Inf)
    for(j in 2:(N+1)){
      if(L[j-1]==0){Q.s[j] <- rtruncnorm(1, Q.s[j-1], Q.s[j+1], 0, 1)}
      if(L[j-1]==1){Q.s[j] <- rtruncnorm(1, Q.s[j-1], Q.s[j+1], mu1, sigma1)}
      if(L[j-1]==2){Q.s[j] <- rtruncnorm(1, Q.s[j-1], Q.s[j+1], mu2, sigma2)}
    }
    y0 <- Q.s[which(L==0)+1]
    y1 <- Q.s[which(L==1)+1]
    y2 <- Q.s[which(L==2)+1]
    sigma1 <- sqrt(1/rgamma(1, 0.5*(n1-1), rate = var(y1)*(n1-1)/2))
    mu1 <- rtruncnorm(1, b=0, mean=mean(y1), sd=sigma1/sqrt(n1))
    sigma2 <- sqrt(1/rgamma(1, 0.5*(n2-1), rate = var(y2)*(n2-1)/2))
    mu2 <- rtruncnorm(1, a=0, mean=mean(y2), sd=sigma2/sqrt(n2))
    
    sample[i,] <- c(mu1, sigma1, mu2, sigma2) 
    sample1[i, ] <- genabcd(c(mu1, sigma1, mu2, sigma2))
    
    if(i %% 500==0){
      par(mfrow = c(4, 1), mar = c(2, 3, 2, 2) + 1)
      for(k in 1:4){
        plot(sample1[1:i, k], type = "l",
             xlab = "iteration", ylab = paste("param", k))
      }
    }
  }
  postmean <- colMeans(sample[(iter.num/10):iter.num,])
  postmedian <- colMedians(sample[(iter.num/10):iter.num,])
  postmean1 <- colMeans(sample1[(iter.num/10):iter.num,])
  postmedian1 <- colMedians(sample1[(iter.num/10):iter.num,])
  return(list(sample=sample, postmean=postmean, postmedian=postmedian, 
              sample1=sample1, postmean1=postmean1, postmedian1=postmedian1))
  
}

# trinormality method under verification bias
logrankvb <- function(z, L, n, iter.num){
  # mu1 <- intval[1]
  # mu2 <- intval[3]
  # sigma1 <- intval[2]
  # sigma2 <- intval[4]
  N <- length(L)
  
  sample <- matrix(0, iter.num, 4)
  sample1 <- matrix(0, iter.num, 4)
  
  paramd <- n
  #prior information is three group distributed evenly and the scale of the param is related to how confident we are about this prior information
  
  p <- rdirichlet(1,paramd) 
  L.sim <- L
  for(j in 1:N){
    if(L[j]==3){
      L.sim[j] <- which(rmultinom(1,1,p)==1)-1
    }
  }
  
  x0 <- z[which(L.sim==0)]
  x1 <- z[which(L.sim==1)]
  x2 <- z[which(L.sim==2)]
  
  mu0 <- mean(x0)
  sigma0 <- sd(x0)
  y0 <- (x0-mu0)/sigma0
  y1 <- (x1-mu0)/sigma0
  y2 <- (x2-mu0)/sigma0
  y20 <- y2
  
  mu1 <- mean(y1)
  mu2 <- mean(y20)
  sigma1 <- sd(y1)
  sigma2 <- sd(y20)
  Q <- c(y0,y1,y2)
  Q.s <- sort(Q)
  Q.s <- c(-Inf, Q.s, Inf)
  
  for(i in 1:iter.num){
    n.ll <- c(sum(L.sim==0), sum(L.sim==1), sum(L.sim==2))
    
    for(j in 2:(N+1)){
      if(L.sim[j-1]==0){Q.s[j] <- rtruncnorm(1, Q.s[j-1], Q.s[j+1], 0, 1)}
      if(L.sim[j-1]==1){Q.s[j] <- rtruncnorm(1, Q.s[j-1], Q.s[j+1], mu1, sigma1)}
      if(L.sim[j-1]==2){Q.s[j] <- rtruncnorm(1, Q.s[j-1], Q.s[j+1], mu2, sigma2)}
    }
    y0 <- Q.s[which(L.sim==0)+1]
    y1 <- Q.s[which(L.sim==1)+1]
    y2 <- Q.s[which(L.sim==2)+1]
    
    sigma1 <- sqrt(1/rgamma(1, 0.5*(n.ll[2]-1), rate = var(y1)*(n.ll[2]-1)/2))
    mu1 <- rtruncnorm(1, b=0, mean=mean(y1), sd=sigma1/sqrt(n.ll[2]))
    sigma2 <- sqrt(1/rgamma(1, 0.5*(n.ll[3]-1), rate = var(y2)*(n.ll[3]-1)/2))
    mu2 <- rtruncnorm(1, a=0, mean=mean(y2), sd=sigma2/sqrt(n.ll[3]))
    
    sample[i,] <- c(mu1, sigma1, mu2, sigma2)
    sample1[i,] <- genabcd(c(mu1, sigma1, mu2, sigma2)) 
    
    p <- rdirichlet(1, paramd+n.ll)
    
    L.sim <- L
    Q.s <- Q.s[-1]
    for(j in 1:N){
      if(L[j]==3){
        pmulti <- rep(0,3)
        dn <- c(dnorm(Q.s[j]),dnorm(Q.s[j], mu1, sigma1),dnorm(Q.s[j], mu2, sigma2))
        pmulti[1] <- (p[1]*dn[1])/(sum(p*dn))
        pmulti[2] <- (p[2]*dn[2])/(sum(p*dn))
        pmulti[3] <- (p[3]*dn[3])/(sum(p*dn))
        L.sim[j] <- which(rmultinom(1,1,pmulti)==1)-1
      }
    }
    Q.s <- c(-Inf, Q.s)
    
    if(i %% 500==0){
      par(mfrow = c(4, 1), mar = c(2, 3, 2, 2) + 1)
      for(k in 1:4){
        plot(sample1[1:i, k], type = "l",
             xlab = "iteration", ylab = paste("param", k))
      }
    }
  }
  postmean <- colMeans(sample[(iter.num/10):iter.num,])
  postmedian <- colMedians(sample[(iter.num/10):iter.num,])
  postmean1 <- colMeans(sample1[(iter.num/10):iter.num,])
  postmedian1 <- colMedians(sample1[(iter.num/10):iter.num,])
  return(list(sample=sample, postmean=postmean, postmedian=postmedian, 
              sample1=sample1, postmean1=postmean1, postmedian1=postmedian1))
}

# read data

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

# data <- cbind(L, CA125, CA153)

index <- c(which(CA125==""), which(HE4==""))
index <- unique(index)
CA125 <- CA125[-index]
HE4 <- HE4[-index]
L <- L[-index]

CA125 <- as.numeric(paste(CA125))
HE4 <- as.numeric(paste(HE4))

L[order(CA125)]
L[order(HE4)]

data <- cbind(L, CA125, HE4)
data <- data[-which(CA125>1000),]

L <- data[,1]
CA125 <- log(data[,2])
HE4 <- log(data[,3])
S <- HE4

L <- L[order(S)]

for(i in 1:length(L)){
  if(is.na(L[i])){L[i]=4}
  if(L[i]==2){L[i]=0}
  if(L[i]==3){L[i]=2}
  if(L[i]==4){L[i]=3}
}

L.new <- data[,1]
for(i in 1:length(L.new)){
  if(is.na(L.new[i])){
    L.new[i] <- sample(c(1:3), size=1, prob=rep(1,3)/3)
  }
}


n0 <- sum(L.new==2)
n1 <- sum(L.new==1)
n2 <- sum(L.new==3)


result <- logrankvb(S, L, c(n0, n1, n2), iter.num=100000)


write.csv(result$sample, file='HE4_musigma.csv')
write.csv(result$sample1, file='HE4_abcd.csv')

genabcd(colMeans(result$sample[20000:100000,]))
colMeans(result$sample1[20000:100000,])
VUS2(colMeans(result$sample1[20000:100000,]))


# > genabcd(colMeans(result$sample[50000:300000,]))
# [1]  1.1300068 -1.2170748  0.8977535  0.7551002
# > colMeans(result$sample1[50000:300000,])
# [1]  1.1754825 -1.2178029  0.9289647  0.7733957
# > VUS2(colMeans(result$sample1[50000:300000,]))
# [1] 0.5157985