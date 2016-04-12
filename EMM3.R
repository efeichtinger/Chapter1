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
library(grid)

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
eig.all["r"] <-log(eig.all$lam)
names(eig.all)[names(eig.all)=='gbar'] <- 'gamma'
str(eig.all)

##Figures - labeled as intended for manuscript (Fig 1 is life cycle)

#Figure 2 - lambda as a function of sigma across different g at 2 phis
#2 figures (or more, can increase to 4 panels) in one, phi = 1 and phi = 0
datA <- subset(eig.all, phi=="1" & jsur == "0.7" & gamma == "0.5")
datB <- subset(eig.all, phi=="0" & jsur == "0.7" & gamma == "0.5")
datAB <- rbind(datA, datB)
datC <- subset(eig.all, phi==0.9 & gamma == 0.1 & jsur == 0.5)
datD <- subset(eig.all, phi==0.9 & gamma == 0.3 & jsur == 0.5)
datE <- subset(eig.all, phi==0.9 & gamma == 0.5 & jsur == 0.5)
datF <- subset(eig.all, phi==0.9 & gamma == 0.7 & jsur == 0.5)
datG <- subset(eig.all, phi==0.9 & gamma == 0.9 & jsur == 0.5)

#phi 1 Code does not work
#p1 <- ggplot(datC, aes(x=sigma, y=lam)) +
#geom_line(aes(color="0.1")) + 
#geom_line(data=datD, aes(color="0.3")) +
#geom_line(data=datE, aes(color="0.5")) +
#geom_line(data=datF, aes(color="0.7")) +
#geom_line(data=datG, aes(color="0.9")) +
  #labs(color = "g", x = expression(sigma), y = expression(lambda))

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


#p2 <- ggplot(datH, aes(x=sigma, y=lam)) +
  #geom_line(aes(color="0.1")) + 
  #geom_line(data=datI, aes(color="0.3")) +
  #geom_line(data=datJ, aes(color="0.5")) +
  #geom_line(data=datK, aes(color="0.7")) +
  #geom_line(data=datL, aes(color="0.9")) +
  #labs(color = "g", x = expression(sigma), y = expression(lambda))

dat.pan <- rbind(datC,datD,datE,datF,datG,datH,datI,datJ,datK,datL)
dat.pan2 <- rbind(datM,datN,datO,datP,datQ,datR,datS,datT,datU,datV)
dat.pan3 <- rbind(datE,datJ,datO,datT)
dat.all <- rbind(dat.pan,dat.pan2)


#Create something for labels first 
#value1 <- c(-0.3,0,0.5,1)
#value2 <- c(0.1,0.3,0.5,0.7,0.9)
#labsx <- list(bquote(gamma==.(value1)),bquote(phi==.(value2)))
#labsx <- list(bquote(phi==.(value1)),bquote(g==.(value2)))

#Label function -  taken from stackoverflow
#http://stackoverflow.com/questions/14181234/facet-labels-involving-a-greek-symbol
my.label <- function (expr1 = gamma == .(x), expr2 = phi == .(x)) 
{
  quoted1<- substitute(expr1)
  quoted2 <- substitute(expr2)
  function(variable, value) {
    value <- as.character(value)
    if(variable == 'gamma')
      lapply(value, function(x)
        eval(substitute(bquote(expr1, list(x = x)),list(expr1 = quoted1))))
    else
      lapply(value, function(x) 
        eval(substitute(bquote(expr2, list(x = x)),list(expr2 = quoted2))))
  }
}

#Figure 2 
#Want to get phi as a symbol and gbar as a g with a line over it

p3 <- ggplot(dat.all, aes(sigma, lam)) + geom_line()
p3 + facet_grid(gamma~ phi, labeller=my.label()) + 
labs(x=expression(sigma), y=expression(lambda)) +
  theme(strip.text.x = element_text(size = 13)) +
  theme(strip.text.y = element_text(size = 13)) + 
  theme(axis.text.x = element_text(size = 11, angle = 45)) +
  theme(axis.text.y = element_text(size = 11)) +
  theme(axis.title.x = element_text(size = 15, face = 'bold')) +
  theme(axis.title.y = element_text(size = 15, face = 'bold'))

p3b <- ggplot(dat.all, aes(sigma, r)) + geom_line()
p3b + facet_grid(gamma~ phi, labeller=my.label()) + 
  labs(x=expression(sigma), y="r") +
  theme(strip.text.x = element_text(size = 13)) +
  theme(strip.text.y = element_text(size = 13)) + 
  theme(axis.text.x = element_text(size = 11, angle = 45)) +
  theme(axis.text.y = element_text(size = 11)) +
  theme(axis.title.x = element_text(size = 15, face = 'bold')) +
  theme(axis.title.y = element_text(size = 15, face = 'bold'))



#Figure 3 - R0 as a function of sigma (similar style as 3)
#change to R0 with subscript
p4 <- ggplot(dat.all, aes(sigma, R0)) + geom_line()
p4 + facet_grid(gamma ~ phi,labeller=my.label()) +
labs(x=expression(sigma), y=expression('R'[0])) +
  theme(strip.text.x = element_text(size = 13)) +
  theme(strip.text.y = element_text(size = 13)) +
  theme(axis.text.x = element_text(size = 11, angle = 45)) +
  theme(axis.text.y = element_text(size = 11)) +
  theme(axis.title.x = element_text(size = 15, face = 'bold')) +
  theme(axis.title.y = element_text(size = 15, face = 'bold'))


#Figure 4 - T as a function of sigma (similar style as 3)
##FIX
p5 <- ggplot(dat.all, aes(sigma, time)) + geom_line()
p5 + facet_grid(gamma ~ phi,labeller=my.label()) +
  labs(x=expression(sigma),y=expression('T'[1])) +
  theme(strip.text.x = element_text(size = 13)) +
  theme(strip.text.y = element_text(size = 13)) +
  theme(axis.text.x = element_text(size = 11, angle = 45)) +
  theme(axis.text.y = element_text(size = 11)) +
  theme(axis.title.x = element_text(size = 15, face = 'bold')) +
  theme(axis.title.y = element_text(size = 15, face = 'bold'))


#Figure 5 - Damping ratio panel (similar style as 3)
p6 <- ggplot(dat.all, aes(sigma, DampR)) + geom_line()
p6 + facet_grid(gamma ~ phi,labeller=my.label()) +
  labs(x=expression(sigma), y="Damping Ratio") +
  #changes text size in panels 
  theme(strip.text.x = element_text(size = 13)) +
  theme(strip.text.y = element_text(size = 13)) +
  theme(axis.text.x = element_text(size = 11, angle = 45)) +
  theme(axis.text.y = element_text(size = 11)) +
  theme(axis.title.x = element_text(size = 15, face = 'bold')) +
  theme(axis.title.y = element_text(size = 15, face = 'bold'))



#Figure 6 - mat plot, SSD (NOTE- might want to use just this and not Fig 5)
#Monotypic populations (this can result from phi = 1, but don't present it this way
## s = 0.5, g = 0.5, sigma = 0.05, P = 0.7
## Slow
sl <- c(0.275,1,0.225,0.7)
slow <- matrix(sl, nrow=2, ncol =2, byrow =TRUE)
fa <- c(0.225, 1, 0.275, 0.7)
fast <- matrix(fa, nrow=2, ncol=2, byrow=TRUE)


colnames(slow) <- c("Slow juvenile", "Slow adult")
colnames(fast) <- c("Fast juvenile", "Fast adult")
rownames(slow) <- c("Slow juvenile", "Slow adult")
rownames(fast) <- c("Fast juvenile", "Fast adult")

tf <- splitA(slow, r = 1, c=2)
Tmat1 <- tf$T
f1 <- fundamental.matrix(Tmat1)

tf2 <- splitA(fast, r =1, c=2)
Tmat2 <- tf2$T
f2 <- fundamental.matrix(Tmat2)

matplot2(pop.projection(slow, c(1,1), 100)$stage.vectors, col= 10:11, 
         lwd = 3, proportions = TRUE, legend= "right")
matplot2(pop.projection(fast, c(1,1), 100)$stage.vectors, col = 10:11,
         lwd = 3, proportions = TRUE, legend= "right")

matplot2(pop.projection(slow, c(1,1), 100)$stage.vectors, col= 10:11, 
         lwd = 3, proportions = FALSE, legend= "topleft")
matplot2(pop.projection(fast, c(1,1), 100)$stage.vectors, col = 10:11,
         lwd = 3, proportions = FALSE, legend= "topleft")

### Population dynamics 

eig.slow <- eigen.analysis(slow)
eig.fast <- eigen.analysis(fast)

R0.slow <- net.reproductive.rate(slow, r = 1, c = 2)
R0.fast <- net.reproductive.rate(fast, r = 1, c = 2)

t.slow <- generation.time(slow, r = 1, c = 2)
t.fast <- generation.time(fast, r = 1, c = 2)

eig.slow
eig.fast
R0.slow
R0.fast
t.slow
t.fast

## monotypic populations at different g's 
#g = 0.1, all other parameters the same 
s2 <- matrix(c(0.475,1,0.025,0.7), nrow = 2, ncol = 2, byrow=TRUE)
f2 <- matrix(c(0.425,1,0.075,0.7), nrow = 2, ncol = 2, byrow=TRUE)
colnames(s2) <- c("Slow juvenile", "Slow adult")
colnames(f2) <- c("Fast juvenile", "Fast adult")
rownames(s2) <- c("Slow juvenile", "Slow adult")
rownames(f2) <- c("Fast juvenile", "Fast adult")
s2
f2

#g = 0.2
s3 <- matrix(c(0.425,1,0.075,0.7), nrow = 2, ncol = 2, byrow=TRUE)
f3 <- matrix(c(0.375,1,0.125,0.7), nrow = 2, ncol = 2, byrow=TRUE)
colnames(s3) <- c("Slow juvenile", "Slow adult")
colnames(f3) <- c("Fast juvenile", "Fast adult")
rownames(s3) <- c("Slow juvenile", "Slow adult")
rownames(f3) <- c("Fast juvenile", "Fast adult")

#g = 0.3
s4 <- matrix(c(0.375,1,0.125,0.7), nrow = 2, ncol = 2, byrow=TRUE)
f4 <- matrix(c(0.325,1,0.175,0.7, nrow = 2, ncol = 2, byrow=TRUE))
colnames(s4) <- c("Slow juvenile", "Slow adult")
colnames(f4) <- c("Fast juvenile", "Fast adult")
rownames(s4) <- c("Slow juvenile", "Slow adult")
rownames(f4) <- c("Fast juvenile", "Fast adult")

#g = 0.4
s5 <- matrix(c(0.325,1,0.175,0.7), nrow = 2, ncol = 2, byrow=TRUE)
f5 <- matrix(c(0.275,1,0.225,0.7), nrow = 2, ncol = 2, byrow=TRUE)
colnames(s5) <- c("Slow juvenile", "Slow adult")
colnames(f5) <- c("Fast juvenile", "Fast adult")
rownames(s5) <- c("Slow juvenile", "Slow adult")
rownames(f5) <- c("Fast juvenile", "Fast adult")

#g = 0.6 
s6 <- matrix(c(0.225,1,0.275,0.7), nrow = 2, ncol = 2, byrow=TRUE)
f6 <- matrix(c(0.175,1,0.325,0.7), nrow = 2, ncol = 2, byrow=TRUE)
colnames(s6) <- c("Slow juvenile", "Slow adult")
colnames(f6) <- c("Fast juvenile", "Fast adult")
rownames(s6) <- c("Slow juvenile", "Slow adult")
rownames(f6) <- c("Fast juvenile", "Fast adult")

#g = 0.7 
s7 <- matrix(c(0.175,1,0.325,0.7), nrow = 2, ncol = 2, byrow=TRUE)
f7 <- matrix(c(0.125,1,0.375,0.7), nrow = 2, ncol = 2, byrow=TRUE)
colnames(s7) <- c("Slow juvenile", "Slow adult")
colnames(f7) <- c("Fast juvenile", "Fast adult")
rownames(s7) <- c("Slow juvenile", "Slow adult")
rownames(f7) <- c("Fast juvenile", "Fast adult")

#g = 0.8 
s8 <- matrix(c(0.125,1,0.375,0.7), nrow = 2, ncol = 2, byrow=TRUE)
f8 <- matrix(c(0.075,1,0.425,0.7), nrow = 2, ncol = 2, byrow=TRUE)
colnames(s8) <- c("Slow juvenile", "Slow adult")
colnames(f8) <- c("Fast juvenile", "Fast adult")
rownames(s8) <- c("Slow juvenile", "Slow adult")
rownames(f8) <- c("Fast juvenile", "Fast adult")

#g = 0.9
s9 <- matrix(c(0.075,1,0.425,0.7), nrow = 2, ncol = 2, byrow=TRUE)
f9 <- matrix(c(0.025,1,0.475,0.7), nrow = 2, ncol = 2, byrow=TRUE)
colnames(s9) <- c("Slow juvenile", "Slow adult")
colnames(f9) <- c("Fast juvenile", "Fast adult")
rownames(s9) <- c("Slow juvenile", "Slow adult")
rownames(f9) <- c("Fast juvenile", "Fast adult")


## Pop dynamics s2-s9 and f2-f9
#g = 0.1
eig.s2 <- eigen.analysis(s2)
eig.f2 <- eigen.analysis(f2)
R0.s2 <- net.reproductive.rate(slow, r = 1, c = 2)
R0.f2 <- net.reproductive.rate(fast, r = 1, c = 2)
t.s2 <- generation.time(slow, r = 1, c = 2)
t.f2 <- generation.time(fast, r = 1, c = 2)

#g = 0.2
eig.s3 <- eigen.analysis(s2)
eig.f3 <- eigen.analysis(f2)
R0.s3 <- net.reproductive.rate(slow, r = 1, c = 2)
R0.f3 <- net.reproductive.rate(fast, r = 1, c = 2)
t.s3 <- generation.time(slow, r = 1, c = 2)
t.f3 <- generation.time(fast, r = 1, c = 2)

#g = 0.3 
eig.s4 <- eigen.analysis(s2)
eig.f4 <- eigen.analysis(f2)
R0.s4 <- net.reproductive.rate(slow, r = 1, c = 2)
R0.f4 <- net.reproductive.rate(fast, r = 1, c = 2)
t.s4 <- generation.time(slow, r = 1, c = 2)
t.f4 <- generation.time(fast, r = 1, c = 2)

#g = 0.4
eig.s5 <- eigen.analysis(s2)
eig.f5 <- eigen.analysis(f2)
R0.s5 <- net.reproductive.rate(slow, r = 1, c = 2)
R0.f5 <- net.reproductive.rate(fast, r = 1, c = 2)
t.s5 <- generation.time(slow, r = 1, c = 2)
t.f5 <- generation.time(fast, r = 1, c = 2)

#g = 0.6 
eig.s6 <- eigen.analysis(s2)
eig.f6 <- eigen.analysis(f2)
R0.s6 <- net.reproductive.rate(slow, r = 1, c = 2)
R0.f6 <- net.reproductive.rate(fast, r = 1, c = 2)
t.s6 <- generation.time(slow, r = 1, c = 2)
t.f6 <- generation.time(fast, r = 1, c = 2)

#g = 0.7
eig.s7 <- eigen.analysis(s2)
eig.f7 <- eigen.analysis(f2)
R0.s7 <- net.reproductive.rate(slow, r = 1, c = 2)
R0.f7 <- net.reproductive.rate(fast, r = 1, c = 2)
t.s7 <- generation.time(slow, r = 1, c = 2)
t.f7 <- generation.time(fast, r = 1, c = 2)

#g = 0.8
eig.s8 <- eigen.analysis(s2)
eig.f8 <- eigen.analysis(f2)
R0.s8 <- net.reproductive.rate(slow, r = 1, c = 2)
R0.f8 <- net.reproductive.rate(fast, r = 1, c = 2)
t.s8 <- generation.time(slow, r = 1, c = 2)
t.f8 <- generation.time(fast, r = 1, c = 2)

#g = 0.9
eig.s9 <- eigen.analysis(s2)
eig.f9 <- eigen.analysis(f2)
R0.s9 <- net.reproductive.rate(slow, r = 1, c = 2)
R0.f9 <- net.reproductive.rate(fast, r = 1, c = 2)
t.s9 <- generation.time(slow, r = 1, c = 2)
t.f9 <- generation.time(fast, r = 1, c = 2)


#### Need to make a data frame of g, lambda, R0, T from the monotypic pops
#################################

#### 7 April 2016 - practice code 
#2 figures (or more, can increase to 4 panels) in one, phi = 1 and phi = 0

#Just one line first for practice 
plotA <- ggplot(datAB, aes(sigma, lam))
plotA + geom_line(colour = "red", linetype="solid", size = 1) + xlim(0,0.5) +
  ylim(1,1.4) + labs(x = expression(sigma), y = expression(lambda)) 

#Start here I guess - multiple phi's on one plot with 1 g and 1 s
ggplot(datA, aes(x=sigma, y =lam)) + 
  geom_line(aes(color="Phi = 1")) +
  geom_line(data=datB, aes(colour="Phi = 0")) + labs(x = expression(sigma),
                                                     y = expression(lambda), color="Legend") 

## Juvenile survival 
dJ1 <- subset(eig.all, phi==1 & gamma == 0.6 & jsur == 0.1)
dJ2 <- subset(eig.all, phi==1 & gamma == 0.6 & jsur == 0.2)
dJ3 <- subset(eig.all, phi==1 & gamma == 0.6 & jsur == 0.3)
dJ4 <- subset(eig.all, phi==1 & gamma == 0.6 & jsur == 0.4)
dJ5 <- subset(eig.all, phi==1 & gamma == 0.6 & jsur == 0.5)
dJ6 <- subset(eig.all, phi==1 & gamma == 0.6 & jsur == 0.6)
dJ7 <- subset(eig.all, phi==1 & gamma == 0.6 & jsur == 0.7)
dJ8 <- subset(eig.all, phi==1 & gamma == 0.6 & jsur == 0.8)
dJ9 <- subset(eig.all, phi==1 & gamma == 0.6 & jsur == 0.9)

d <- rbind(dJ1, dJ2, dJ3, dJ4, dJ5, dJ6, dJ7, dJ8, dJ9)


pJ <- ggplot(d, aes(sigma, lam)) + geom_line()
pJ + facet_grid(jsur~ phi, labeller=label_both) + 
  labs(x=expression(sigma), y=expression(lambda)) +
  theme(strip.text.x = element_text(size = 13)) +
  theme(strip.text.y = element_text(size = 13)) + 
  theme(axis.text.x = element_text(size = 11, angle = 45)) +
  theme(axis.text.y = element_text(size = 11)) +
  theme(axis.title.x = element_text(size = 15, face = 'bold')) +
  theme(axis.title.y = element_text(size = 15, face = 'bold'))

dK1 <- subset(eig.all, phi==0 & gamma == 0.6 & jsur == 0.1)
dK2 <- subset(eig.all, phi==0 & gamma == 0.6 & jsur == 0.2)
dK3 <- subset(eig.all, phi==0 & gamma == 0.6 & jsur == 0.3)
dK4 <- subset(eig.all, phi==0 & gamma == 0.6 & jsur == 0.4)
dK5 <- subset(eig.all, phi==0 & gamma == 0.6 & jsur == 0.5)
dK6 <- subset(eig.all, phi==0 & gamma == 0.6 & jsur == 0.6)
dK7 <- subset(eig.all, phi==0 & gamma == 0.6 & jsur == 0.7)
dK8 <- subset(eig.all, phi==0 & gamma == 0.6 & jsur == 0.8)
dK9 <- subset(eig.all, phi==0 & gamma == 0.6 & jsur == 0.9)

dd <- rbind(dK1,dK2,dK3,dK4,dK5,dK6,dK7,dK8,dK9)

pd <- ggplot(dd, aes(sigma, lam)) + geom_line()
pd + facet_grid(jsur~ phi, labeller=label_both) + 
  labs(x=expression(sigma), y=expression(lambda)) +
  theme(strip.text.x = element_text(size = 13)) +
  theme(strip.text.y = element_text(size = 13)) + 
  theme(axis.text.x = element_text(size = 11, angle = 45)) +
  theme(axis.text.y = element_text(size = 11)) +
  theme(axis.title.x = element_text(size = 15, face = 'bold')) +
  theme(axis.title.y = element_text(size = 15, face = 'bold'))
