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

V <- as.numeric(!is.na(L))
library(bcROCsurface)
# Preparing the missing disease status
Dna <- preDATA(L, S)
Dvec.na <- Dna$Dvec
Dfact.na <- Dna$D

# FI estimator
rho.out <- rhoMLogit(Dfact.na ~ CA125+HE4, test = TRUE)
ROCs("fi", T = S, Dvec = Dvec.na, V = V, rhoEst = rho.out, ncp = 50)
VUS.FI <- vus(method="fi", T=S, Dvec=Dvec.na, V=V, rhoEst = rho.out)

# MSI estimator
ROCs("msi", T = S, Dvec = Dvec.na, V = V, rhoEst = rho.out, ncp = 50)
VUS.MSI <- vus("msi", T=S, Dvec=Dvec.na, V=V, rhoEst = rho.out)

#IPW estimator
pi.out <- psglm(V ~ CA125+HE4, test = TRUE)
ROCs("ipw", T = S, Dvec = Dvec.na, V = V, piEst = pi.out, ncp=50)
VUS.IPW <- vus("ipw", T = S, Dvec = Dvec.na, V = V, piEst = pi.out)

# SPE estimator
ROCs("spe", T = S, Dvec = Dvec.na, V = V, rhoEst = rho.out, piEst = pi.out, ncp=50)
VUS.SPE <- vus("spe", T = S, Dvec = Dvec.na, V = V, rhoEst = rho.out, piEst = pi.out)

