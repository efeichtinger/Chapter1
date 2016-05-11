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
library(plyr)


##########################################
# 5/10 Change to a function so the dimensions can be changed
# Define Identity matrix for use in finding mean age at 1st reproduction
# hard code
#I <- matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),
             #nrow = 4, ncol = 4, byrow=TRUE)

# Function for flexibility 
Idmat <- function(dmns){
  i <- diag(1, dmns)
  return(i)
}

#e <- matrix(c(1,1,1,1), nrow = 4, ncol = 1)

# function for e flexibility in size 
# e is a column vector with dimensions the number of stages in model
# input of dmns is numeric
# output a column vector 
evec <- function(dmns){
  e <- matrix(rep(1,dmns,nrow=dmns,ncol=1))
  return(e)
}


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

I <- Idmat(4)
e <- evec(4)

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
S <- 6
for (h in 7:8){
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
          S*gbar, 0, 0, 0,
          0, 0, S*(1-(gbar + sigma)), 0,
          0, 0, S*gbar, 0))
        
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
        B <- B.prime(Mp, I, Tp)
        b2 <- B[2,1:4]
        Tc <- tc.mat(b2, Tp)
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
names(eig.dat) <- c("gbar", "jsur", "sigma", "lam","eigen2", "phi", "r", "R0","DampR", "time", "mages", "magef")

##################################
#Mean of the means across all parameter values (phi doesn't matter here)
# May 3 2016, S = 0.4 
#new.data <- write.csv(eig.dat, file= "May3b.csv")
new.data <- read.csv("May3b.csv", header=TRUE)
names(new.data) <- c("X", "gbar","jsur","sigma","lam","eigen2","phi",
            "r","R0","damp","time","ages","agef")

# May 10 2016, S = 0.5 & 0.6 
data.56 <- read.csv("May10.csv", header = TRUE)
names(data.56) <- c("X", "gbar","jsur","sigma","lam","eigen2","phi",
                    "r","R0","damp","time","ages","agef")

d41 <- subset(new.data, gbar == 0.1, select = gbar:agef)
d42 <- subset(new.data, gbar == 0.3, select = gbar:agef)
d43 <- subset(new.data, gbar == 0.5, select = gbar:agef)
d44 <- subset(new.data, gbar == 0.7, select = gbar:agef)
d45 <- subset(new.data, gbar == 0.9, select = gbar:agef)


d51 <- subset(data.56, gbar == 0.1 & jsur == 0.5, select = gbar:agef)
d52 <- subset(data.56, gbar == 0.3 & jsur == 0.5, select = gbar:agef)
d53 <- subset(data.56, gbar == 0.5 & jsur == 0.5, select = gbar:agef)
d54 <- subset(data.56, gbar == 0.7 & jsur == 0.5, select = gbar:agef)
d55 <- subset(data.56, gbar == 0.9 & jsur == 0.5, select = gbar:agef)

#d5all <- rbind(d51,d52,d53,d54,d55)

d61 <- subset(data.56, gbar == 0.1 & jsur == 0.6, select = gbar:agef)
d62 <- subset(data.56, gbar == 0.3 & jsur == 0.6, select = gbar:agef)
d63 <- subset(data.56, gbar == 0.5 & jsur == 0.6, select = gbar:agef)
d64 <- subset(data.56, gbar == 0.7 & jsur == 0.6, select = gbar:agef)
d65 <- subset(data.56, gbar == 0.9 & jsur == 0.6, select = gbar:agef)

#d6all <- rbind(d61,d62,d63,d64,d65)

dall <- rbind(d41,d42,d43,d44,d45,d51,d52,d53,d54,d55,d61,d62,d63,d64,d65)

#S = 0.5
px <- ggplot(dall, aes(x=sigma, y=ages)) + geom_line()
px + facet_grid(jsur~gbar, labeller = label_both) +
  geom_line(data=dall, aes(x=sigma, y=agef)) +
  labs(x=expression(sigma), y="Mean Age at 1st Reproduction")



####################################################

#subset by S 
sur5 <- subset(data.56, jsur == 0.5)
sur6 <- subset(data.56, jsur == 0.6)

# S = 0.5
mins5 <- min(sur5$mages)
minf5 <- min(sur5$magef)
maxs5 <- max(sur5$mages)
maxf5 <- max(sur5$magef)
ms5 <- mean(sur5$mages)
mf5 <-mean(sur5$magef)
sds5 <-sd(sur5$mages)
sdf5 <- sd(sur5$magef)

# S = 0.6
mins6 <- min(sur6$mages)
minf6 <- min(sur6$magef)
maxs6 <- max(sur6$mages)
maxf6 <- max(sur6$magef)
ms6 <- mean(sur6$mages)
mf6 <-mean(sur6$magef)
sds6 <-sd(sur6$mages)
sdf6 <- sd(sur6$magef)


#slows 0.5
sl5 <- cbind(ms5, sds5, mins5, maxs5, 0.5)
sl5 <- as.data.frame(sl5)
sl5["Type"] <- "Slow"
colnames(sl5)[1:5] <- c("mean", "sd", "min", "max","S")
 

#slows 0.6 
sl6 <- cbind(ms6, sds6, mins6, maxs6, 0.6)
sl6 <- as.data.frame(sl6)
sl6["Type"] <- "Slow"
colnames(sl6)[1:5] <- c("mean", "sd", "min", "max", "S")

slow <- rbind(sl5,sl6)

#fasts 0.5

fa5 <- cbind(mf5, sdf5, minf5, maxf5, 0.5)
fa5 <- as.data.frame(fa5)
fa5["Type"] <- "Fast"
colnames(fa5)[1:5] <- c("mean", "sd", "min", "max","S")

#fasts 0.6
fa6 <- cbind(mf6, sdf6, minf6, maxf6, 0.6)
fa6 <- as.data.frame(fa6)
fa6["Type"] <- "Fast"
colnames(fa6)[1:5] <- c("mean", "sd", "min", "max","S")

fast <- rbind(fa5,fa6)

sum.stats <- rbind(slow, fast)
sum.stats

###########################################

### S  = 0.4 

#Figures to show the mean age at first reproduction 
#We can only have one more figure in the manuscript (6 is the limit)

#drop phi 
phis <- names(new.data) %in% c("X","pc")
new.data <- new.data[!phis]

#Subset for graphing  - values for Phi (pc) and gbar match other figures
## Trying to solve the problem 
## Two data frames, one with slow mean age and the other with fast 
d1s <- subset(new.data, gbar == 0.1, select = gbar:ages)
d2s <- subset(new.data, gbar == 0.3, select = gbar:ages)
d3s <- subset(new.data, gbar == 0.5, select = gbar:ages)
d4s <- subset(new.data, gbar == 0.7, select = gbar:ages)
d5s <- subset(new.data, gbar == 0.9, select = gbar:ages)

#slow
dat.alls <- rbind(d1s,d2s,d3s,d4s,d5s)

#fast
nam <- names(new.data) %in% c("ages")
new.data <- new.data[!nam]

d1f <- subset(new.data, gbar == 0.1)
d2f <- subset(new.data, gbar == 0.3)
d3f <- subset(new.data, gbar == 0.5)
d4f <- subset(new.data, gbar == 0.7)
d5f <- subset(new.data, gbar == 0.9)

#fast all
dat.allf <- rbind(d1f,d2f,d3f,d4f,d5f)
              
##################################
#Label function -  taken from stackoverflow
#http://stackoverflow.com/questions/14181234/facet-labels-involving-a-greek-symbol
my.label <- function (expr1 = gamma == .(x), expr2 = sigma == .(x)) 
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
###################################

# Phi doesn't matter so adjust figures to remove phi
#facet grid 


## Plots

ggplot() +
  geom_line(data=dat.alls, aes(x=sigma, y=ages), color = "green") +
  geom_line(data= dat.allf, aes(x=sigma, y=agef), color="blue")

#S = 0.4
pz <- ggplot(dat.alls, aes(x=sigma, y=ages)) + geom_line()
pz + facet_grid(.~gbar, labeller = my.label()) +
  geom_line(data=dat.allf,aes(x=sigma, y=agef)) +
  labs(x=expression(sigma), y="Mean Age at 1st Reproduction")


