### Attempt to add the mean age at first reproduction
### Code taken from EMM3.R
### That simulation works so I didn't want to mess with it 

#parameters s, gbar, sigma, phi, F, P
#sigma < gbar
#gbar is a probability - so it's positive and < 1
#0 < phi < 1
#S and P also on (0,1)
#F is fertility - positive

#clear objects from environment for simulation output (large)
rm(list = ls())

library(popbio)
library(ggplot2)

#Define Identity matrix for use in finding mean age at 1st reproduction
I <- matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),
             nrow = 4, ncol = 4, byrow=TRUE)

#Define column vector for use in finding mean age at 1st reproduction
e <- matrix(c(1,1,1,1), nrow = 4, ncol = 1)


### Functions

# function for estimating Bienvenus's T
Bien.T <- function(L, v, w, r){
  time <- L * drop(v%*%w)/(drop(v%*%r%*%w))
  return(time)
}


#This function is Eq (5.51) in Caswell, all entries all matrices
#input for m is M' - Caswell Eq. 5.55
#input for i is Identity matrix 
#input for t is T' - Caswell Eq. 5.55
B.prime <- function(m,i,t){
  b <- m %*% (solve(i - t))
  return(b) 
}


#Function computes T(c), inputs and output are matrices
#the input for b is b2
#the input for t is the transition matrix 
#output is a matrix, the T(c) matrix 
tc.mat <- function(b,t){
  tc <- solve(diag(b)) %*% t %*% diag(b)
  return(tc)
}

#Function computes the mean age at first reprodcuction
#input at i is the identity matrix 
#input at e is a column vector of 1's with dim 1 X s 
#input t is tc.mat 
mean.age <- function(e,i,t){
  ma <- ((t(e)) %*% (solve(i-t)))
  return(ma)
}

#####################################################

# Simulation 
# May 2016

#Fix P and F
P <- 0.7
F <- 1

#create vectors to store output
gb <- vector()
sig <- vector()
eig <- vector()
eig2 <- vector()
poc <- vector()
jsur <- vector()
r <- vector()
genT <- vector()
nrepd <- vector()
dpr <- vector()
gbr <- vector()
age <- vector()
MES <- vector()
MEF <- vector()
gbar <- 0
ii <- 0

#Simulation 

gb <- 0.01 * (1:99)
S <- 5
for (h in 6:6){
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
        
        #transition matrix for 2 type model
        tr <-matrix(ncol=4, nrow=4, byrow=TRUE, data=c(
          S*(1-(gbar-sigma)), 0, 0, 0,
          S * (gbar-sigma), P, 0, 0,
          0, 0, S*(1 - (gbar + sigma)), 0,
          0, 0, S * (gbar + sigma), P))
        
        Tp <- matrix(ncol=4, nrow=4, byrow=TRUE, data = c(
          S*(1-(gbar-sigma)), 0, 0, 0,
          S*(gbar-sigma), 0, 0, 0,
          0, 0, S*(1-(gbar + sigma)), 0,
          0, 0, S*(gbar + sigma), 0))
        
        Mp <- matrix(ncol=4, nrow=2, byrow=TRUE, 
                          data = c(1-S, 0, 1-S, 0,
                                   0, 1, 0, 1))
        
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
        B <- B.prime(Mp, I, tr)
        b2 <- B[2,1:4]
        Tc <- tc.mat(b2,tr)
        Tc[2,2] <- 0
        age <- mean.age(e, I, Tc)
        #the output of age is a row vector, I want the 1st and 3rd entries 
        #in different objects so it can collect 
        MES[ii] <- age[1,1]
        MEF[ii] <- age[1,3]
        
        
      }
    }
  }
}

eig.dat <- data.frame(gbar = gbr, js = jsur, sigma = sig, eigen = eig, eigen2 = eig2, pc = poc, instr = r, R0 = nrepd, damp = dpr, time = genT, ages =MES, agef =MEF)

#new.data <- write.csv(eig.dat, file= "data.set.csv")

new.data <- read.csv("data.set.csv", header=TRUE)

#Subset for graphing
d1 <- subset(new.data, pc == -0.3 & gbar == 0.1)
d2 <- subset(new.data, pc == -0.3 & gbar == 0.3)
d3 <- subset(new.data, pc == -0.3 & gbar == 0.5)
d4 <- subset(new.data, pc == -0.3 & gbar == 0.7)
d5 <- subset(new.data, pc == -0.3 & gbar == 0.9)
d6 <- subset(new.data, pc == 0 & gbar == 0.1)
d7 <- subset(new.data, pc == 0 & gbar == 0.3)
d8 <- subset(new.data, pc == 0 & gbar == 0.5)
d9 <- subset(new.data, pc == 0 & gbar == 0.7)
d10 <- subset(new.data, pc == 0 & gbar == 0.9)
d11 <- subset(new.data, pc == 0.5 & gbar == 0.1)
d12 <- subset(new.data, pc == 0.5 & gbar == 0.3)
d13 <- subset(new.data, pc == 0.5 & gbar == 0.5)
d14 <- subset(new.data, pc == 0.5 & gbar == 0.7)
d15 <- subset(new.data, pc == 0.5 & gbar == 0.9)
d16 <- subset(new.data, pc == 1 & gbar == 0.1)
d17 <- subset(new.data, pc == 1 & gbar == 0.3)
d18 <- subset(new.data, pc == 1 & gbar == 0.5)
d19 <- subset(new.data, pc == 1 & gbar == 0.7)
d20 <- subset(new.data, pc == 1 & gbar == 0.9)

dat.all <- rbind(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16,d17,
          d18,d19,d20)

p <- ggplot(dat.all, aes(sigma, agef)) + geom_line()
p + facet_grid(gbar ~ pc, labeller=label_parsed) +
  labs(x=expression(sigma), y="age") 

        
