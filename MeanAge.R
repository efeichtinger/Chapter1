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
library(reshape2)
library(scales)


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
for (h in 10:10){
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
names(eig.dat) <- c("gamma", "S", "sigma", "lam","eigen2", "phi", "r", "R0","DampR", "time", "ages", "agef")

##################################
#Mean of the means across all parameter values (phi doesn't matter here)
## May 12 2016 - S = 0.1 and 0.2 
#write.csv(eig.dat, file="May12.csv")

# S = 0.1 and 0.2 
data.12 <- read.csv("May12.csv", header=TRUE)
names(data.12) <- c("X", "gbar","jsur","sigma","lam","eigen2","phi",
                    "r","R0","damp","time","ages","agef")

may12 <- subset(data.12,  jsur == 0.1)
may12b <- subset(data.12, jsur == 0.2)

mean(may12$ages)
mean(may12$agef)
sd(may12$ages)
sd(may12$agef)
mean(may12b$ages)
mean(may12b$agef)
sd(may12b$ages)
sd(may12b$agef)
range(may12$ages)
range(may12$agef)
range(may12b$ages)
range(may12b$agef)

fml <- may12$ages
fmlm <- may12$agef
f <- data.frame(c(fml,fmlm))
colnames(f) <- "age"
mean(f$age)
sd(f$age)


sht <- may12b$ages
shty <- may12b$agef
s <- data.frame(c(sht,shty))
colnames(s) <- "age"
mean(s$age)
sd(s$age)

##############################################################

# May 3 2016, S = 0.4 
# Read in data from output 
#new.data <- write.csv(eig.dat, file= "May3b.csv")
new.data <- read.csv("May3b.csv", header=TRUE)
names(new.data) <- c("X", "gamma","S","sigma","lam","eigen2","phi",
            "r","R0","damp","time","ages","agef")

phis3 <- subset(new.data, new.data$phi == -0.3 & gamma == 0.5)
phis0 <- subset(new.data, new.data$phi == 0 & gamma == 0.5)
phis5 <- subset(new.data, new.data$phi == 0.5 & gamma == 0.5)
phis9 <- subset(new.data, new.data$phi == 0.9 & gamma == 0.5)

may17 <- rbind(phis3, phis0, phis5, phis9)


ggplot(data = may17, aes(x=sigma, y=lam)) + 
  geom_line() + 
  theme_bw() + 
  facet_grid(. ~phi) + 
  labs(x = expression(sigma), y=expression(lambda)) + 
  theme(axis.text.x = element_text(angle = 45, size = 8))

ggplot(data = may17, aes(x=sigma, y=R0)) + 
  geom_line() + 
  theme_bw() + 
  facet_grid(. ~phi) + 
  labs(x = expression(sigma), y="Net Reproductive Rate") + 
  theme(axis.text.x = element_text(angle = 45, size = 8))

ggplot(data = may17, aes(x=sigma, y=time)) + 
  geom_line() + 
  theme_bw() + 
  facet_grid(. ~phi) + 
  labs(x = expression(sigma), y="Generation Time") + 
  theme(axis.text.x = element_text(angle = 45, size = 8))


# May 10 2016, S = 0.5 & 0.6 
# Read in data from output
data.56 <- read.csv("May10.csv", header = TRUE)
names(data.56) <- c("X", "gamma","S","sigma","lam","eigen2","phi",
                    "r","R0","damp","time","ages","agef")

#### June 3 2016 - Happy birthday to self although this work fills me with rage
# Subset data where sigma = 0 to see relationship between g and lambda, R0, T
# Using data.56 and new.data
no.sig4 <- subset(new.data, sigma ==0)
no.sig5 <- subset(data.56, sigma == 0 & S == 0.5)
no.sig6 <- subset(data.56, sigma == 0 & S == 0.6)

new <- rbind(no.sig4,no.sig5,no.sig6)
new$S <- as.factor(new$S)

## plots 

# Lambda
ggplot(data = new, aes(x=gamma, y=lam, fill=S)) + geom_line()

## 5 15 2017 
ggplot(data = new, aes(x=gamma, y = lam, colour=S)) +
  geom_line(size=1) + 
  labs(x="Individual Growth Rate", y = expression(lambda)) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_color_manual(breaks = c("0.4","0.5","0.6"), 
                     values = c("seagreen4", "seagreen3", "seagreen2"))

ggplot(data = new, aes(x= gamma, y = R0, colour=S)) + 
  geom_line(size = 1) + 
  labs(x="Individual Growth Rate", y = "Net Reproductive Rate") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_color_manual(breaks = c("0.4","0.5","0.6"), 
                     values = c("seagreen4", "seagreen3", "seagreen2"))

ggplot(data = new, aes(x= gamma, y = time, colour=S)) + 
  geom_line(size = 1) + 
  labs(x="Individual Growth Rate", y = "Generation Time") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_color_manual(breaks = c("0.4","0.5","0.6"), 
                     values = c("seagreen4", "seagreen3", "seagreen2"))


#png("Rplot01.png", width = 7.48, height = 5.31 , units="in", res = 300)

ggplot() + 
  geom_line(data=no.sig5, aes(x=gamma, y=lam), color = "green") +
  geom_line(data=no.sig6, aes(x =gamma, y=lam), color= "blue") +
  geom_line(data=no.sig4, aes(x=gamma, y=lam), color= "red") +
  labs(x="Growth Rate", y =expression(lambda)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggplot() + 
  geom_line(data=no.sig5, aes(x=gamma, y=R0), color = "green") +
  geom_line(data=no.sig6, aes(x =gamma, y=R0), color= "blue") +
  geom_line(data=no.sig4, aes(x=gamma, y=R0), color= "red")

ggplot() +
  geom_line(data=no.sig5, aes(x=gamma, y=time), color = "green") +
  geom_line(data=no.sig6, aes(x =gamma, y=time), color= "blue") +
  geom_line(data=no.sig4, aes(x=gamma, y=time), color= "red")





#####################################
#pooled mean
#pms <- new.data$ages
#pmf <- new.data$agef
#pm <- data.frame(c(pms,pmf))
#colnames(pm) <- "age"
#mean(pm$age)
########################################

### Reformat data by adding indicator columns for each phenotype
### I used the brute force approach although there are probably functions
### that could do it 

# S= 0.5 and 0.6  
data.56["typef"] <- 1
data.56["types"] <- 0
#new.data$typef <- as.factor(new.data$typef)
#new.data$types <- as.factor(new.data$types)

# Make two new dataframes with a type indicator and the mean age for each type
# 0 for slow
slows56 <- cbind(data.56$types, data.56$ages)
colnames(slows56) <- c("type", "meanage")
# 1 for fast
fasts56 <- cbind(data.56$typef, data.56$agef)
colnames(fasts56) <- c("type", "meanage")

#bind them together so they are "stacked", all slows then all fasts 
both56 <- rbind(slows56, fasts56)

#subset the full data frame without the mean age columns, twice
fst56 <- subset(data.56, select = gamma:time)
scd56 <- subset(data.56, select = gamma:time)

#bind to duplicate 
doub56 <- rbind(scd56, fst56)

#Add columns of type and mean age 
stack56 <- cbind(doub56, both56)

#add factor levels 
stack56$type <- as.factor(stack56$type)
levels(stack56$type) <- c(levels(stack56$type), c("slow","fast"))

stack56$type[stack56$type ==0] <- "slow"
stack56$type[stack56$type ==1] <- "fast"



#######
# Just trying something
 # S  = 0.4 
# Make two new dataframes with a type indicator and the mean age for each type
new.data["types"] <- 0
new.data["typef"] <- 1

slows <- cbind(new.data$types, new.data$ages)
colnames(slows) <- c("type", "meanage")
fasts <- cbind(new.data$typef, new.data$agef)
colnames(fasts) <- c("type", "meanage")

#bind them together so they are "stacked", all slows then all fasts 
both <- rbind(slows, fasts)

#subset the full data frame without the mean age columns, twice
fst <- subset(new.data, select = gamma:time)
scd <- subset(new.data, select = gamma:time)

#bind to duplicate 
doub <- rbind(fst,scd)

#Add columns of type and mean age 
stack <- cbind(doub, both)
stack$type <- as.factor(stack$type)

#add factor levels 
levels(stack$type) <- c(levels(stack$type), c("slow","fast"))

stack$type[stack$type ==0] <- "slow"
stack$type[stack$type ==1] <- "fast"


dt.all <- rbind(stack,stack56)


#subset by g 
da <- subset(dt.all, gamma == 0.1)
db <- subset(dt.all, gamma == 0.3)
dc <- subset(dt.all, gamma == 0.5)
de <- subset(dt.all, gamma == 0.7)
df <- subset(dt.all, gamma == 0.9)

#bind back together 
togtr <- rbind(da,db,dc,dc,de,df)


#Plot so we have levels for g and phenotype  
#yikes <- ggplot(togtr, aes(x=sigma, y=meanage)) + geom_line()
#yikes + facet_grid(gbar~type, labeller = label_both) +
  #labs(x=expression(sigma), y="Mean Age 1st Reproduction")
  
# Now do this for the other values of S and add a level for S? 

####################################################################
#Subset for graphing  - values for gbar match other figures
##################################
#Label function -  taken from stackoverflow
#http://stackoverflow.com/questions/14181234/facet-labels-involving-a-greek-symbol
my.label <- function (expr1 = gamma == .(x), expr2 = S == .(x)) 
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


d41 <- subset(new.data, gamma == 0.1, select = gamma:agef)
d42 <- subset(new.data, gamma == 0.3, select = gamma:agef)
d43 <- subset(new.data, gamma == 0.5, select = gamma:agef)
d44 <- subset(new.data, gamma == 0.7, select = gamma:agef)
d45 <- subset(new.data, gamma == 0.9, select = gamma:agef)

d4all <- rbind(d41,d42,d43,d44,d45)

d51 <- subset(data.56, gamma == 0.1 & S == 0.5, select = gamma:agef)
d52 <- subset(data.56, gamma == 0.3 & S == 0.5, select = gamma:agef)
d53 <- subset(data.56, gamma == 0.5 & S == 0.5, select = gamma:agef)
d54 <- subset(data.56, gamma == 0.7 & S == 0.5, select = gamma:agef)
d55 <- subset(data.56, gamma == 0.9 & S == 0.5, select = gamma:agef)

d5all <- rbind(d51,d52,d53,d54,d55)

d61 <- subset(data.56, gamma == 0.1 & S == 0.6, select = gamma:agef)
d62 <- subset(data.56, gamma == 0.3 & S == 0.6, select = gamma:agef)
d63 <- subset(data.56, gamma == 0.5 & S == 0.6, select = gamma:agef)
d64 <- subset(data.56, gamma == 0.7 & S == 0.6, select = gamma:agef)
d65 <- subset(data.56, gamma == 0.9 & S == 0.6, select = gamma:agef)

d6all <- rbind(d61,d62,d63,d64,d65)

#bind all 
d456all <- rbind(d4all,d5all,d6all)

newnew <- subset(d456all, gamma==0.5)

px <- ggplot(newnew, aes(x=sigma, y=ages)) + geom_line()
px + facet_grid(S~gamma, labeller = my.label()) +
  geom_line(data=newnew, aes(x=sigma, y=agef), linetype=2) +
  labs(x=expression(sigma), y="Mean age at 1st reproduction") +
  theme(strip.text.x = element_text(size = 13)) +
  theme(strip.text.y = element_text(size = 13)) + 
  theme(axis.text.x = element_text(size = 11, angle = 45)) +
  theme(axis.text.y = element_text(size = 11)) +
  theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 13)) +
  scale_x_continuous(breaks=pretty_breaks(n=3)) +
  scale_y_continuous(breaks=pretty_breaks(n=4)) +
  theme_bw()


#####################################################
#S = 0.4, 0.5, 0.6
px <- ggplot(d456all, aes(x=sigma, y=ages)) + geom_line()
px + facet_grid(S~gamma, labeller = my.label()) +
  geom_line(data=d456all, aes(x=sigma, y=agef), linetype=2) +
  labs(x=expression(sigma), y="Mean age at 1st reproduction") +
  theme(strip.text.x = element_text(size = 13)) +
  theme(strip.text.y = element_text(size = 13)) + 
  theme(axis.text.x = element_text(size = 11, angle = 45)) +
  theme(axis.text.y = element_text(size = 11)) +
  theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 13)) +
  scale_x_continuous(breaks=pretty_breaks(n=3)) +
  scale_y_continuous(breaks=pretty_breaks(n=4)) +
  theme_bw()

############################################################

#togtr

pzz <- ggplot(togtr, 
        aes(x=sigma, y=meanage, linetype = type)) + geom_line() 
pzz + facet_grid(S~gamma, labeller = my.label()) +
labs(x=expression(sigma), y="Mean age at 1st reproduction") +
  theme(strip.text.x = element_text(size = 13)) +
  theme(strip.text.y = element_text(size = 13)) + 
  theme(axis.text.x = element_text(size = 11, angle = 45)) +
  theme(axis.text.y = element_text(size = 11)) +
  theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 13)) +
  scale_x_continuous(breaks=pretty_breaks(n=3)) +
  scale_y_continuous(breaks=pretty_breaks(n=4))

        
####################################################

#subset by S 
sur5 <- subset(data.56, S == 0.5)
sur6 <- subset(data.56, S == 0.6)

##################################
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
##############################################

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


