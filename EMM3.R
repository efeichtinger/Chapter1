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
for (h in 10:11){
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
#eigen9 <- write.csv(eig.new, file= "eigen9.csv")
#eigen <- write.csv(eig.new, file = "eigen.csv")

#Read in all stored files 
eig.2 <- read.csv("eigen2.csv", head=TRUE)
str(eig.2)
#tmp <- read.csv("eigen2.csv", nrows = 150500)
eig.3 <- read.csv("eigen.csv", head=TRUE)
str(eig.3)
eig.4 <- read.csv("eigenapril.csv", head=TRUE)
str(eig.4)
eig.5 <- read.csv("eigen9.csv", head=TRUE)
str(eig.5)


eig.all <- rbind(eig.4,eig.2,eig.3,eig.5)
str(eig.all)
eig.all["r"] <-log(eig.all$lam)
names(eig.all)[names(eig.all)=='gbar'] <- 'gamma'

##Figures - labeled as intended for manuscript (Fig 1 is life cycle)

#Figure 2 - lambda as a function of sigma across different g at 2 phis
#2 figures (or more, can increase to 4 panels) in one, phi = 1 and phi = 0
datA <- subset(eig.all, phi=="1" & jsur == "0.7" & gamma == "0.5")
datB <- subset(eig.all, phi=="0" & jsur == "0.7" & gamma == "0.5")
datAB <- rbind(datA, datB)

#Just one line first for practice 
plotA <- ggplot(datAB, aes(sigma, lam))
plotA + geom_line(colour = "red", linetype="solid", size = 1) + xlim(0,0.5) +
  ylim(1,1.4) + labs(x = expression(sigma), y = expression(lambda)) 

#Start here I guess - multiple phi's on one plot with 1 g and 1 s
ggplot(datA, aes(x=sigma, y =lam)) + 
  geom_line(aes(color="Phi = 1")) +
  geom_line(data=datB, aes(colour="Phi = 0")) + labs(x = expression(sigma),
        y = expression(lambda), color="Legend") 

#Could do with multiple g's (and s?) for 4 different phi's in a panel 
datC <- subset(eig.all, phi==1 & gamma == 0.1 & jsur == 0.5)
datD <- subset(eig.all, phi==1 & gamma == 0.3 & jsur == 0.5)
datE <- subset(eig.all, phi==1 & gamma == 0.5 & jsur == 0.5)
datF <- subset(eig.all, phi==1 & gamma == 0.7 & jsur == 0.5)
datG <- subset(eig.all, phi==1 & gamma == 0.9 & jsur == 0.5)

#phi 1
p1 <- ggplot(datC, aes(x=sigma, y=lam)) +
geom_line(aes(color="0.1")) + 
geom_line(data=datD, aes(color="0.3")) +
geom_line(data=datE, aes(color="0.5")) +
geom_line(data=datF, aes(color="0.7")) +
geom_line(data=datG, aes(color="0.9")) +
  labs(color = "g", x = expression(sigma), y = expression(lambda))

#phi 0
datH <- subset(eig.all, phi==0& gamma == 0.1 & jsur == 0.5)
datI <- subset(eig.all, phi==0 & gamma == 0.3 & jsur == 0.5)
datJ <- subset(eig.all, phi==0 & gamma == 0.5 & jsur == 0.5)
datK <- subset(eig.all, phi==0 & gamma == 0.7 & jsur == 0.5)
datL <- subset(eig.all, phi==0 & gamma == 0.9 & jsur == 0.5)

datM <- subset(eig.all, phi == 0.5 & gamma==0.1 & jsur==0.5)
datN <- subset(eig.all, phi == 0.5 & gamma==0.3 & jsur==0.5)
datO <- subset(eig.all, phi ==0.5 & gamma==0.5 & jsur==0.5)
datP <- subset(eig.all, phi == 0.5 & gamma==0.7 & jsur==0.5)
datQ <- subset(eig.all, phi ==0.5 & gamma==0.9 & jsur==0.5)

datR <- subset(eig.all, phi == -0.3 & gamma==0.1 & jsur==0.5)
datS <- subset(eig.all, phi == -0.3 & gamma==0.3 & jsur==0.5)
datT <- subset(eig.all, phi == -0.3 & gamma==0.5 & jsur==0.5)
datU <- subset(eig.all, phi == -0.3 & gamma==0.7 & jsur==0.5)
datV <- subset(eig.all, phi == -0.3 & gamma==0.9 & jsur==0.5)


p2 <- ggplot(datH, aes(x=sigma, y=lam)) +
  geom_line(aes(color="0.1")) + 
  geom_line(data=datI, aes(color="0.3")) +
  geom_line(data=datJ, aes(color="0.5")) +
  geom_line(data=datK, aes(color="0.7")) +
  geom_line(data=datL, aes(color="0.9")) +
  labs(color = "g", x = expression(sigma), y = expression(lambda))

dat.pan <- rbind(datC,datD,datE,datF,datG,datH,datI,datJ,datK,datL)
dat.pan2 <- rbind(datM,datN,datO,datP,datQ,datR,datS,datT,datU,datV)
dat.pan3 <- rbind(datE,datJ,datO,datT)
dat.all <- rbind(dat.pan,dat.pan2)


#Create something for labels first 
value1 <- c(-0.3,0,0.5,1)
value2 <- c(0.1,0.3,0.5,0.7,0.9)
labsx <- list(bquote(g==.(value1)),bquote(phi==.(value2)))
labsx <- list(bquote(phi==.(value1)),bquote(g==.(value2)))

#mf labeller

#Figure 2 
#Want to get phi as a symbol and gbar as a g with a line over it 
p3 <- ggplot(dat.all, aes(sigma, lam)) + geom_line()
p3 + facet_grid(gamma~ phi, labeller=label_parsed) + 
labs(x=expression(sigma), y=expression(lambda))


#Figure 3 - R0 as a function of sigma (similar style as 3)
#change to R0 with subscript
p4 <- ggplot(dat.all, aes(sigma, R0)) + geom_line()
p4 + facet_grid(gamma ~ phi) +
labs(x=expression(sigma), y="Net Reproductive Rate")

#Figure 4 - T as a function of sigma (similar style as 3)
##FIX
p5 <- ggplot(dat.all, aes(sigma, time)) + geom_line()
p5 + facet_grid(gamma ~ phi) +
  labs(x=expression(sigma),y="Generation Time")

#Figure 5 - Damping ratio panel (similar style as 3)
p6 <- ggplot(dat.all, aes(sigma, DampR)) + geom_line()
p6 + facet_grid(gamma ~ phi) +
  labs(x=expression(sigma), y="Damping Ratio") +
  #changes text size in panels 
theme(strip.text.x = element_text(size = 13)) +
  theme(strip.text.y = element_text(size = 13))


#Figure 6 - mat plot, SSD (NOTE- might want to use just this and not Fig 5)
#Monotypic populations (this can result from phi = 1, but don't present it this way
## s = 0.5, g = 0.5, sigma = 0.05, P = 0.7
## Slow
sl <- c(0.275,1,0.225,0.7)
slow <- matrix(sl, nrow=2, ncol =2, byrow =TRUE)
fa <- c(0.225, 1, 0.275, 0.7)
fast <- matrix(fa, nrow=2, ncol=2, byrow=TRUE)


colnames(slow) <- c("slow juvenile", "slow adult")
colnames(fast) <- c("fast juvenile", "fast adult")
rownames(slow) <- c("slow juvenile", "slow adult")
rownames(fast) <- c("fast juvenile", "fast adult")

tf <- splitA(slow, r = 1, c=2)
Tmat1 <- tf$T
f1 <- fundamental.matrix(Tmat1)

tf2 <- splitA(fast, r =1, c=2)
Tmat2 <- tf2$T
f2 <- fundamental.matrix(Tmat2)

matplot2(pop.projection(slow, c(1,1), 100)$stage.vectors, col= 10:11, 
         lwd = 3, proportions = TRUE, legend= "topright")
matplot2(pop.projection(fast, c(1,1), 100)$stage.vectors, col = 10:11,
         lwd = 3, proportions = TRUE, legend= "topright")

matplot2(pop.projection(slow, c(1,1), 100)$stage.vectors, col= 10:11, 
         lwd = 3, proportions = FALSE, legend= "topright")
matplot2(pop.projection(fast, c(1,1), 100)$stage.vectors, col = 10:11,
         lwd = 3, proportions = FALSE, legend= "topright")

