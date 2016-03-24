## UPDATED March 24 2016 - moved code from "EMM2.R" 
## This file contains much of the same information from EMM2, but cleaner
## Matrix model for growth rate het research
## Dissertation chapter 1


#parameters s, gbar, sigma, phi, F, P
#sigma < gbar
#gbar is a probability - so it's positive and < 1
#0 < phi < 1
#S and P also on (0,1)
#F is fertility - positive

#clear objects from environment for simulation output (large)
rm(list = ls())

library(colorRamps)
library(RColorBrewer)
library(colorspace)
library(graphics)
library(popbio)
library(scatterplot3d)
library(ggplot2)
library(plot3D)

#function for estimating Bienvenus's T
Bien.T <- function(L, v, w, r){
  time <- L * drop(v%*%w)/(drop(v%*%r%*%w))
  return(time)
}

#Fix P and F
P <- 0.7
F <- 1

#create vectors to store output
gb <- vector()
sig <- vector()
eig <- vector()
eig2 <- vector()
#poc- vector to store parent-offspring correlation phi
poc <- vector()
jsur <- vector()
r <- vector()
genT <- vector()
nrepd <- vector()
dpr <- vector()
BT <- vector()
fd <- vector()
#fund1 <- vector()
#fund2 <- vector()
#fund3 <- vector()
#fund4 <- vector()
gbr <- vector()

gbar <- 0
ii <- 0

#Simulation 

gb <- 0.01 * (1:99)
S <- 0
for (h in 8:9){
  S <- (h - 2) * 0.1
  phi <- 0
  for (k in 0:20){
    phi <- (k - 10) * 0.1
    for (i in 1:99){
      gbar <- gb[i]
      for (j in 0:((100 * min (gbar, 1- gbar)-1))){
        ii <- ii + 1
        sigma <- 0.01 * j
        A <-matrix(ncol=4, nrow=4, byrow=TRUE, data=c(
          S*(1-(gbar-sigma)), (1 + phi) * F/2, 0, (1-phi) * F/2,
          S * (gbar-sigma), P, 0, 0,
          0, (1 - phi) *F/2, S*(1 - (gbar + sigma)), (1 + phi) * F/2,
          0, 0, S * (gbar + sigma), P) )
        
        gbr[ii] <- gbar
        jsur[ii] <- S
        sig[ii] <- sigma
        eig[ii] <- Re(eigen(A)[[1]][[1]])
        eig2[ii] <- Re(eigen(A)[[1]][[2]])
        r[ii] <- Re(log(eig[ii]))
        #r[ii] <- Re(log(Re(eigen(A)[[1]][[1]])))
        genT[ii] <- generation.time(A, r = c(1,3), c = c(2,4))
        nrepd[ii] <- net.reproductive.rate(A, r = c(1,3), c = c(2,4))
        dpr[ii] <- damping.ratio(A)
        poc[ii] <- phi
        #
        split.mat <- splitA(A, r = c(1,3), c = c(2,4))
        fert.mat <- split.mat$F
        ssd <- eigen(A)$vectors[,1]/sum(eigen(A)$vectors[,1])
        rvv <- (eigen(t(A))$vectors[,1]/eigen(t(A))$vectors[1,1])
        bienvenu <- Bien.T(eig[ii], rvv, ssd, fert.mat)
        BT[ii] <- bienvenu
        fund <- fundamental.matrix(A, split.mat$T)
        #fd[ii]<- fund
        #spits out first element so this is only for slow juveniles staying slow
        #I'm not sure how to store the output
        #fund1[ii] <- fund$N
        #fund2[ii] <- fund$var
        #fund3[ii] <- fund$meaneta
        #fund4[ii] <- fund$vareta
      }
    }
  }
}

eig.dat <- data.frame(gbar = gbr, js = jsur, sigma = sig, eigen = eig, eigen2 = eig2, pc = poc, instr = r, R0 = nrepd, damp = dpr, time = genT, BVT = BT)


##New data frame without T2/Bienvenu's T 
eig.new <- eig.dat[,1:10]
names(eig.new) <- c("gbar", "jsur", "sigma", "lam","eigen2", "phi", "r", "R0","DampR", "time")
head(eig.new)
## Simulation output with s  = 0.6 and 0.7
eigen <- write.csv(eig.new, file = "eigen.csv")

#Read in stored file with s = 0.3, 0.4, and 0. 5
eig.2 <- read.csv("eignew2.csv", head=TRUE)

##Have eig.new and eig.2, former has s = 0.6,0.7 and latter is 0.3-0.5

##Figures - labeled as intended for manuscript (Fig 1 is life cycle)

#Figure 2 - lambda as a function of sigma across different g at 2 phis
#2 figures (or more, can increase to 4 panels) in one, phi = 1 and phi = 0
datA <- subset(eig.new, phi=="1" & jsur == "0.7" & gbar == "0.5")
datB <- subset(eig.new, phi=="0" & jsur == "0.7" & gbar == "0.5")
plotA <- ggplot(datA, aes(datA$sigma, datA$lam))
plotA + geom_line(colour = "red", linetype="solid", size = 1) + xlim(0,0.5) +
  ylim(1,1.4) + labs(x = expression(sigma), y = expression(lambda)) 

ggplot(datA, aes(x=sigma, y =lam)) + geom_line(aes(color="Phi")) +
  geom_line(data=datB, aes(colour="blue")) + labs(x = expression(sigma),
                                    y = expression(lamba))

#Figure 3 - R0 as a function of sigma (similar style as 3)

#Figure 4 - T as a function of sigma (similar style as 3)

#Figure 5 - Damping ratio panel (similar style as 3)

#Figure 6 - mat plot, SSD (NOTE- might want to use just this and not Fig 5)

