##UPDATED JANUARY 2016

##UPDATED October 2015 
#April 2015 - Most recently updated script for matrix model, pop dyn, & figures

#parameters s, gbar, sigma, phi, F, P

#sigma < gbar
#gbar is a probability - so it's positive and < 1
#0 < phi < 1
#S and P also on (0,1)
#F is fertility - positive

#IT WORKS!!!!!!!

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


P <- 0.7
F <- 1

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

gb <- 0.01 * (1:99)
S <- 0
for (h in 8:8){
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
#fund.mat <- matrix(fund, byrow=TRUE)
#fund.mat <- data.frame(fund.mat= fund1, ncol = 64, byrow=TRUE)
#var.time <- matrix(fund2)
#deatht <- matrix(fund3)
#vardeath <- matrix(fund4)
#fdm <- fund

##New data frame without T2/Bienvenu's T 
eig.new2 <- eig.dat[,1:10]
names(eig.new2) <- c("gbar", "jsur", "sigma", "lam","eigen2", "phi", "r", "R0","DampR", "time")
#eig.new2$r <- NULL
cor(eig.new2)
eigen2[1:10,c(1,2,4,3,5,6)]
eigen2 <- eigen2[,c(1,2,4,3,5,6)]
head(eig.new2)
#write.table(eig.new2, file = "eignew2.txt", sep="\t")
#write.csv(eig.new2, file="eignew3.txt", sep="")

cols <- rainbow(10, s= 1, v =1, start = 0, end = 1, alpha = 1)
cols2 <- heat.colors(100,alpha =1)
#ugly
cols6 <- brewer.pal(5, "Spectral")
#ok
cols7 <- brewer.pal(8, "PRGn")
#ok but adjust starting point
cols8 <- brewer.pal(4, "Blues")

#Figures of T, R0 and r
rplot <- scatter2D(eig.new2$R0, eig.new2$time, colvar=eig.new2$r,
            type = "points", xlab = "R0", ylab="T", clab = "r", xlim = c(0,2),
            ylim = c(3,9), pin = c(5,5), col = cols8)
#Gordon liked T on x axis 
rplotb <- scatter2D(eig.new2$time, eig.new2$R0, colvar=eig.new2$r,
          xlab = "T", ylab ="R0", clab="r",
          xlim = c(3.5,9), ylim = c(0,2), pin = c(5,5), col = cols)
rplotc <- scatter2D(eig.new2$time, eig.new2$r, colvar = eig.new2$R0)
rplotd <- scatter2D(eig.new2$R0, eig.new2$r, colvar = eig.new2$T,
                    xlab = "R0", ylab="r", clab="T")

ggplot(pres.cand, aes(pres.cand$Name, y="Polling results", color = "Media Outlet")) + 
  geom_point(aes(y = pres.cand$ABC_politicalpollresults , col = "ABC")) + 
  geom_point(aes(y = pres.cand$NBC_politicalpollresults, col = "NBC")) + 
  labs(x = "Candidate", y = "Polling Results")   +
  theme(
    axis.title=element_text(face="bold", size = "12", color="black"))
p <- ggplot(eig.new2, aes(eig.new2$time), y="r") +
  geom_point(aes(y = eig.new2$r))

#rplotc <- scatterplot3d(eig.new2$time, eig.new2$R0, eig.new2$r, color = 
    #par("col"), pch = 8)


#unnecessary 
lamplotb <- scatter2D(eig.new2$R0, eig.new2$time, colvar = eig.new2$lam,
              xlab = "R0", ylab = "T", clab= "Lambda")
lamplotc <- scatter2D(eig.new2$gbar, eig.new2$R0, 
              colvar = eig.new2$lam, ylab = "R0", xlab = "Growth", 
              clab="Lambda")
lamplotd <- scatter2D(eig.new2$time, eig.new2$R0, colvar=eig.new2$lam,
              xlab = "T", ylab = "R0", clab = "lambda")


#panel graphs for reports only 
par(mfrow=c(2,1))
par(mar=c(4,4,2,2))
rplotc <- scatter2D(eig.new2$r, eig.new2$R0, colvar = eig.new2$time, 
                    xlab = "r", ylab = "R0", clab = "T")
rplotd <- scatter2D(eig.new2$r, eig.new2$time, colvar = eig.new2$R0, 
                    xlab = "r", ylab = "T", clab = "R0")

lamplot2 <- scatter3D(eig.new2$gbar, eig.new2$R0, eig.new2$time, 
                      colvar=eig.new2$lam,
                      xlab = "g", ylab = "R0", zlab = "T", ticktype = "detailed",
                      clab=c("lambda"))


#Read in CSV file of simulation output
eig2 <- read.csv("eignew2.csv")
eig2$r <- c(log(eig2$lam))
eig2b <- subset(eig2, jsur == 0.5, select=gbar:r)

r.plot <- scatter2D(eig2b$time, eig2b$R0, colvar = eig2b$r, xlab = "T", 
          ylab = "R0", clab = "r", xlim = c(3.5,9))

#EDA 3D plots for pltting T and R0 
#lamplot <- scatter3D(eig.new2$R0, eig.new2$time, eig.new2$phi, 
            #colvar=eig.new2$lam,
            #xlab = "R0", ylab = "T", zlab = "Phi", ticktype = "detailed",
            #clab=c("lambda"))

##subsets for graphing
####

#phi=0
phi0.g1 <- subset(eig2, phi==0 & jsur==0.3 & gbar==0.1, select=gbar:time)
phi0.g2 <- subset(eig2, phi==0 & jsur==0.3 & gbar==0.2, select=gbar:time)
phi0.g4 <- subset(eig2, phi==0 & jsur==0.3 & gbar==0.4, select=gbar:time)
#newd <- subset(phi.0.dat, gbar == 0.2, select=gbar:time)
phi0.g3 <- subset(eig2, phi==0 & jsur==0.3 & gbar==0.3, select=gbar:time)
phi0.g5 <- subset(eig2, phi==0 & jsur==0.3 & gbar==0.5, select=gbar:time)
phi0.g6 <- subset(eig2, phi==0 & jsur==0.3 & gbar==0.6, select=gbar:time)
#### Need to make panel figures for R0, T and damping ratio
## Follow code format below 
#legend("bottomright", title ="g", c("0.1","0.3","0.5","0.7","0.9"),lty= c(1,2,3,4,5), cex = 0.75)
### And make a graph of the matplots for ssd 
par(mar=c(4,4,2,1))
par(mfrow=c(2,2))
plot(phi0.g1$sigma, phi0.g1$R0, type = "l", lty = 1, xlim=c(0, 0.49),
     ylim = c(0,1.0), xlab = expression(sigma = sigma), ylab = "R0",
     main = "phi=0")
lines(phi0.g2$sigma, phi0.g2$R0, lty = 2)
lines(phi0.g3$sigma, phi0.g3$R0, lty = 3)
lines(phi0.g4$sigma, phi0.g4$R0, lty = 4)
lines(phi0.g5$sigma, phi0.g5$R0, lty = 5)

#phi=0.3
phi3.g1 <- subset(eig2, phi==0.3 & jsur==0.3 & gbar==0.1, select=gbar:time)
phi3.g2 <- subset(eig2, phi==0.3 & jsur==0.3 & gbar==0.2, select=gbar:time)
phi3.g3 <- subset(eig2, phi==0.3 & jsur==0.3 & gbar==0.3, select=gbar:time)
phi3.g4 <- subset(eig2, phi==0.3 & jsur==0.3 & gbar==0.4, select=gbar:time)
phi3.g5 <- subset(eig2, phi==0.3 & jsur==0.3 & gbar==0.5, select=gbar:time)

plot(phi3.g1$sigma, phi3.g1$R0, type = "l", lty = 1, xlim=c(0, 0.49),
     ylim = c(0,1.0), xlab = expression(sigma = sigma), ylab = "R0",
     main = "phi=0.3")
lines(phi3.g2$sigma, phi3.g2$R0, lty = 2)
lines(phi3.g3$sigma, phi3.g3$R0, lty = 3)
lines(phi3.g4$sigma, phi3.g4$R0, lty = 4)
lines(phi3.g5$sigma, phi3.g5$R0, lty = 5)


#phi=0.5
phi5.g1 <- subset(eig2, phi==0.5 & jsur==0.3 & gbar==0.1, select=gbar:time)
phi5.g2 <- subset(eig2, phi==0.5 & jsur==0.3 & gbar==0.2, select=gbar:time)
phi5.g3 <- subset(eig2, phi==0.5 & jsur==0.3 & gbar==0.3, select=gbar:time)
phi5.g4 <- subset(eig2, phi==0.5 & jsur==0.3 & gbar==0.4, select=gbar:time)
phi5.g5 <- subset(eig2, phi==0.5 & jsur==0.3 & gbar==0.5, select=gbar:time)

plot(phi5.g1$sigma, phi5.g1$R0, type = "l", lty = 1, xlim=c(0, 0.49),
     ylim = c(0,1.0), xlab = expression(sigma = sigma), ylab = "R0",
     main = "phi = 0.5")
lines(phi5.g2$sigma, phi5.g2$R0, lty = 2)
lines(phi5.g3$sigma, phi5.g3$R0, lty = 3)
lines(phi5.g4$sigma, phi5.g4$R0, lty = 4)
lines(phi5.g5$sigma, phi5.g5$R0, lty = 5)


#phi=0.8
phi8.g1 <- subset(eig2, phi==0.8 & gbar==0.1 & jsur==0.3, select=gbar:time)
phi8.g2 <- subset(eig2, phi==0.8 & gbar==0.2 & jsur==0.3, select=gbar:time)
phi8.g3 <- subset(eig2, phi==0.8 & gbar==0.3 & jsur==0.3, select=gbar:time)
phi8.g4 <- subset(eig2, phi==0.8 & gbar==0.4 & jsur==0.3, select=gbar:time)
phi8.g5 <- subset(eig2, phi==0.8 & gbar==0.5 & jsur==0.3, select=gbar:time)

plot(phi8.g1$sigma, phi8.g1$R0, type = "l", lty = 1, xlim=c(0, 0.49),
     ylim = c(0,1.0), xlab = expression(sigma = sigma), ylab = "R0",
     main = "phi = 0.8")
lines(phi8.g2$sigma, phi8.g2$R0, lty = 2)
lines(phi8.g3$sigma, phi8.g3$R0, lty = 3)
lines(phi8.g4$sigma, phi8.g4$R0, lty = 4)
lines(phi8.g5$sigma, phi8.g5$R0, lty = 5)
legend("bottomright", title ="g", c("0.1","0.2","0.3","0.4","0.5"),lty= c(1,2,3,4,5), cex = 0.75)

#phi=1
phi1.g1 <- subset(eig2, phi==1 & gbar==0.1 & jsur==0.5, select=gbar:time)
phi1.g2 <- subset(eig2, phi==1 & gbar==0.2 & jsur==0.5, select=gbar:time)
phi1.g3 <- subset(eig2, phi==1 & gbar==0.3 & jsur==0.5, select=gbar:time)
phi1.g4 <- subset(eig2, phi==1 & gbar==0.4 & jsur==0.5, select=gbar:time)
phi1.g5 <- subset(eig2, phi==1 & gbar==0.5 & jsur==0.5, select=gbar:time)

plot(phi1.g1$sigma, phi1.g1$R0, type = "l", lty = 1, xlim=c(0, 0.49),
     ylim = c(0,1.6), xlab = expression(sigma = sigma), ylab = "R0",
     main = "phi=1")
lines(phi1.g2$sigma, phi1.g2$R0, lty = 2)
lines(phi1.g3$sigma, phi1.g3$R0, lty = 3)
lines(phi1.g4$sigma, phi1.g4$R0, lty = 4)
lines(phi1.g5$sigma, phi1.g5$R0, lty = 5)


#phi=-1
negphi.g1 <- subset(eig2, phi==-1 & gbar==0.4 & jsur==0.3, select=gbar:time)
negphi.g2 <- subset(eig2, phi==-1 & gbar==0.4 & jsur==0.3, select=gbar:time)
negphi.g3 <- subset(eig2, phi==-1 & gbar==0.4 & jsur==0.3, select=gbar:time)
negphi.g4 <- subset(eig2, phi==-1 & gbar==0.4 & jsur==0.3, select=gbar:time)
negphi.g5 <- subset(eig2, phi==-1 & gbar==0.4 & jsur==0.3, select=gbar:time)

## not advisable to use 
plot(negphi.g1$sigma, negphi.g1$R0, type = "l", lty = 1, xlim=c(0, 0.49),
     ylim = c(0,1.6), xlab = expression(sigma = sigma), ylab = "R0",
     main = expression(phi = phi))
lines(negphi.g2$sigma, negphi.g2$R0, lty = 2)
lines(negphi.g3$sigma, negphi.g3$R0, lty = 3)
lines(negphi.g4$sigma, negphi.g4$R0, lty = 4)
lines(negphi.g5$sigma, negphi.g5$R0, lty = 5)

### T
plot(phi1.g5$sigma, phi1.g5$time, type = "l", lty = 1, xlim=c(0, 0.49),
     ylim = c(0,10), xlab = expression(sigma = sigma), ylab = "T",
     main = "phi=1")
lines(phi1.g2$sigma, phi1.g2$time, lty = 2)
lines(phi1.g3$sigma, phi1.g3$time, lty = 3)
lines(phi1.g4$sigma, phi1.g4$time, lty = 4)
lines(phi1.g5$sigma, phi1.g5$time, lty = 5)

#### Make panel graphs of the damping ratio over sigma from different phi's
## Jan 2016
plot(phi0.g4$sigma, phi0.g4$DampR, ylim=c(1.22,1.25), xlab = expression(sigma=sigma), ylab = 
       "Damping Ratio", main = "Casey where phi=0 and g=0.4")
plot(phi5.g4$sigma, phi5.g4$DampR, xlab = expression(sigma=sigma), ylab = 
       "Damping Ratio", main="Case where phi=0.5 and g=0.4")
plot(phi1.g4$sigma, phi1.g4$DampR, ylim=c(1,2),xlab = expression(sigma=sigma), ylab = 
       "Damping Ratio", main="Case where phi=1, and g=0.4")
plot(phi8.g4$sigma, phi8.g4$DampR,ylim=c(1,2), xlab = expression(sigma=sigma), ylab = 
       "Damping Ratio", main="Case where phi=0.8, and g=0.4")
plot(negphi$sigma, negphi$DampR, xlab = expression(sigma=sigma), ylab = 
  "Damping Ratio", main="Phi =-1 and g = 0.4")
plot(negphi$sigma, negphi$lam, xlab = expression(sigma=sigma), ylab = 
       expression(lambda=lambda), main="Phi = -1 and g=0.4")
######




### 1/18/2016 Not sure if these figures will be helpful 
#first and second eigenvalues 
phi5 <- subset(eig2, phi==0.5 & gbar==0.5, select=gbar:time)
plot(phi5$sigma, phi5$eigen2, ylim =c(0.7, 1), xlab = expression(sigma=sigma), 
ylab = "Second eigenvalue", main = "Phi = 0.5")
lines(phi5$sigma, phi5$lam)
phi0 <- subset(eig2, phi==0 & gbar==0.5, select=gbar:time)
plot(phi0$sigma, phi0$eigen2, ylim=c(0.68,1),xlab = expression(sigma=sigma), 
ylab = "eigenvalue")
lines(phi0$sigma, phi0$lam)
phi1 <- subset(eig2, phi==1 & gbar==0.5, select=gbar:time)
plot(phi1$lam, phi1$eigen2, ylim=c(0.75,1.05),xlab = "dominant eigenvalue", 
ylab = "eigenvalue", main = "phi = 1")
lines(phi1$sigma, phi1$lam)
neg <- subset(eig2, phi==-1 & gbar==0.5, select=gbar:time)
plot(neg$sigma, neg$eigen2, ylim=c(0.45,1),xlab = expression(sigma=sigma), 
ylab = "eigenvalue")
lines(neg$sigma, neg$lam)

plot(phi5$lam, phi5$eigen2, ylim =c(0.7, 1), xlab = "Dominant Eigenvalue", 
     ylab = "Second eigenvalue", main = "Phi = 0.5")

##########################
#Fixed values to see what happens
## - 1/17/2016
## s = 0.5, g = 0.4, sigma = 0.05, phi = 0, P = 0.7
b1 <- c(0.325, 0.5, 0, 0.5, 0.175, 0.7, 0, 0, 0, 0.5, 0.275, 0.5, 0, 0, 0.225,0.7)
B1 <- matrix(b1, nrow=4, ncol=4, byrow=TRUE)

##s = 0.5, g = 0.4, sigma = 0.1, phi = 0, P = 0.7
b2 <- c(0.35,0.5,0,0.5,0.15,0.7,0,0,0,0.5,0.25,0.5,0,0,0.25,0.7)
B2 <- matrix(b2, nrow=4, ncol=4, byrow=TRUE)

## s= 0.5, g = 0.4, sigma = 0.05, phi = 0.8, P=0.7
b3 <- c(0.325, 0.9, 0, 0.1, 0.175, 0.7, 0, 0, 0, 0.1, 0.275, 0.9, 0,0,0.225,0.7)
B3 <- matrix(b3, nrow=4, ncol=4, byrow=TRUE)

## s= 0.5, g = 0.4, sigma = 0.1, phi = 0.8, P=0.7
b4 <- c(0.35,0.9,0,0.1,0.15,0.7,0,0,0,0.1,0.25,0.9,0,0,0.25,0.7)
B4 <- matrix(b4, nrow=4, ncol=4, byrow=TRUE)

## s = 0.5, g = 0.4, sigma = 0.05, phi = 1, P = 0.7
## Slow
sl <- c(0.325, 1, 0.175, 0.7)
slow <- matrix(sl, nrow=2, ncol =2, byrow =TRUE)
fa <- c(0.275, 1, 0.225, 0.7)
fast <- matrix(fa, nrow=2, ncol=2, byrow=TRUE)


## s = 0.5, g=0.4, sigma = 0.1, phi = 1, P = 0.7
sl1 <- c(0.3, 1, 0.15, 0.7)
slow1 <- matrix(sl1, nrow=2, ncol=2, byrow=TRUE)
fa1 <- c(0.25, 1, 0.25, 0.7)
fast1 <- matrix(fa1, nrow=2, ncol=2, byrow=TRUE)
        

##
colnames(B1) <- c("slow juvenile", "slow adult", "fast juvenile", "fast adult")
rownames(B1) <- c("slow juvenile", "slow adult", "fast juvenile", "fast adult")
colnames(B2) <- c("slow juvenile", "slow adult", "fast juvenile", "fast adult")
rownames(B2) <- c("slow juvenile", "slow adult", "fast juvenile", "fast adult")
colnames(B3) <- c("slow juvenile", "slow adult", "fast juvenile", "fast adult")
rownames(B3) <- c("slow juvenile", "slow adult", "fast juvenile", "fast adult")
colnames(B4) <- c("slow juvenile", "slow adult", "fast juvenile", "fast adult")
rownames(B4) <- c("slow juvenile", "slow adult", "fast juvenile", "fast adult")
colnames(slow) <- c("slow juvenile", "slow adult")
rownames(slow) <- c("slow juvenile", "slow adult")
colnames(slow1) <- c("slow juvenile", "slow adult")
rownames(slow1) <- c("slow juvenile", "slow adult")
colnames(fast) <- c("fast juvenile", "fast adult")
rownames(fast) <- c("fast juvenile", "fast adult")
colnames(fast1) <- c("fast juvenile", "fast adult")
rownames(fast1) <- c("fast juvenile", "fast adult")

eigen.analysis(B1)
eigen.analysis(B2)
eigen.analysis(B3)
eigen.analysis(B4)
eigen.analysis(slow)
eigen.analysis(fast)
eigen.analysis(slow1)
eigen.analysis(fast1)
net.reproductive.rate(fast)
net.reproductive.rate(slow)
generation.time(fast)
generation.time(slow)

#fundamental matrix
tf <- splitA(B1, r = c(1,3), c(2,4))
Tmat1 <- tf$T
f1 <- fundamental.matrix(Tmat1)

tf2 <- splitA(B2, r = c(1,3), c(2,4))
Tmat2 <- tf2$T
f2 <- fundamental.matrix(Tmat2)

tf3 <- splitA(B3, r = c(1,3), c(2,4))
Tmat3 <- tf3$T
f3 <- fundamental.matrix(Tmat3)

tf4 <- splitA(slow, r = 1, c=2)
Tmat4 <- tf4$T
f4 <- fundamental.matrix(Tmat4)

tf5 <- splitA(slow1, r = 1, c=2)
Tmat5 <- tf5$T
f5 <- fundamental.matrix(Tmat5)

tf6 <- splitA(slow1, r = 1, c=2)
Tmat6 <- tf6$T
f6 <- fundamental.matrix(Tmat6)

tf7 <- splitA(fast, r = 1, c=2)
Tmat7 <- tf7$T
f7 <- fundamental.matrix(Tmat7)

tf8 <- splitA(fast1, r = 1, c=2)
Tmat8 <- tf8$T
f8 <- fundamental.matrix(Tmat8)


par(mar=c(3.9,3.9,1,1))
par(mfrow=c(1,2))
matplot2(pop.projection(B1, c(1,1,1,1), 100)$stage.vectors, col= 10:14, 
         lwd = 3, proportions = TRUE, legend= "topright")
matplot2(pop.projection(B2, c(1,1,1,1), 100)$stage.vectors, col= 10:14, 
         lwd = 3, prop=FALSE, legend= "topright")
matplot2(pop.projection(B3, c(1,1,1,1), 100)$stage.vectors, col= 10:14, 
         lwd = 3, prop=FALSE, legend= "topright")

### not working 
matplot2(pop.projection(slow, c(1,1,1,1), 40)$stage.vectors, col= 10:14, 
         lwd = 3, prop=FALSE, legend= "topright")
matplot2(pop.projection(slow1, c(1,1,1,1), 40)$stage.vectors, col= 10:14, 
         lwd = 3, prop=FALSE, legend= "topright")
matplot2(pop.projection(fast, c(1,1,1,1), 40)$stage.vectors, col= 10:14, 
         lwd = 3, prop=FALSE, legend= "topright")
matplot2(pop.projection(fast1, c(1,1,1,1), 40)$stage.vectors, col= 10:14, 
         lwd = 3, prop=FALSE, legend= "topright")




#A1: s= 0.4, g = 0.4, sigma = 0.1, phi = 0, P = 0.7
x1 <- c(0.28, 0.5,0, 0.5, 0.12, 0.7, 0, 0, 0, 0.5, 0.2, 0.5, 0, 0, 0.2, 0.7)
A1 <- matrix(x1, nrow=4, ncol=4, byrow=TRUE)
A1
#A2: s = 0.4, g=0.4, sigma = 0.1, phi = 0.5, P = 0.7
x2 <- c(0.28, 0.75, 0, 0.25, 0.12, 0.7, 0, 0, 0, 0.25, 0.2, 0.75, 0, 0, 0.2, 0.7)
A2 <- matrix(x2, nrow=4, ncol=4, byrow=TRUE)
A2
#A3: s = 0.4, g = 0.4, sigma = 0.05, phi = 0.5, P = 0.7)
x3 <- c(0.26, 0.75, 0, 0.25, 0.14, 0.7, 0, 0, 0, 0.25, 0.22, 0.75, 0, 0, 0.18, 0.7)
A3 <- matrix(x3, nrow=4, ncol=4, byrow=TRUE)
A3
#A4: s = 0.4, g = 0.4, sigma = 0.05, phi = 0.5, P = 0.7
x4 <- c(0.26, 1, 0, 0, 0.14, 0.7, 0, 0, 0, 0, 0.22, 1,0, 0, 0.18, 0.7)
A4 <- matrix(x4, nrow=4, ncol=4, byrow=TRUE)
A4
#A5: s = 0.4, g = 0.4, sigma = 0.1, phi = -1, P = 0.7
x5 <- c(0.28,0,0,1,0.12,0.7,0,0,0,1,0.2,0,0,0,0.2,0.7)
A5 <- matrix(x5, nrow=4, ncol=4, byrow=TRUE)
colnames(A1) <- c("slow juvenile", "slow adult", "fast juvenile", "fast adult")
rownames(A1) <- c("slow juvenile", "slow adult", "fast juvenile", "fast adult")
colnames(A2) <- c("slow juvenile", "slow adult", "fast juvenile", "fast adult")
rownames(A2) <- c("slow juvenile", "slow adult", "fast juvenile", "fast adult")
colnames(A3) <- c("slow juvenile", "slow adult", "fast juvenile", "fast adult")
rownames(A3) <- c("slow juvenile", "slow adult", "fast juvenile", "fast adult")
colnames(A4) <- c("slow juvenile", "slow adult", "fast juvenile", "fast adult")
rownames(A4) <- c("slow juvenile", "slow adult", "fast juvenile", "fast adult")
colnames(A5) <- c("slow juvenile", "slow adult", "fast juvenile", "fast adult")
rownames(A5) <- c("slow juvenile", "slow adult", "fast juvenile", "fast adult")

tf <- splitA(A1, r = c(1,3), c(2,4))
Tmat <- tf$T
f1 <- fundamental.matrix(Tmat)

tf2 <- splitA(A2, r=c(1,3), c(2,4))
Tmat2 <- tf2$T
f2 <- fundamental.matrix(Tmat2)

tf3 <- splitA(A3, r = c(1,3), c(2,4))
Tmat3 <- tf3$T
f3 <- fundamental.matrix(Tmat3)

tfA4 <- splitA(A4,r = c(1,3), c = c(2,4))
Tmat4 <- tfA4$T
fundamental.matrix(Tmat4)

par(mar=c(3.9,3.9,1,1))
matplot2(pop.projection(A1, c(1,1,1,1), 25)$stage.vectors, col= 10:14, lwd = 3, prop=TRUE, legend= "topright")
matplot2(pop.projection(A2, c(1,1,1,1), 25)$stage.vectors, col = 10:14, lwd = 3, prop=TRUE, legend="topright")
matplot2(pop.projection(A3, c(1,1,1,1), 25)$stage.vectors, col = 10:14, lwd =3, prop=TRUE, legend="topright")
#matplot2(pop.projection(A4, c(1,1,1,1), 25)$stage.vectors, prop=TRUE, legend=NA)
matplot2(pop.projection(A5, c(1,1,1,1), 25)$stage.vectors,col= 10:14, lwd = 3, prop=TRUE, legend="topright")

#A6: s = 0.4, g = 0.4, sigma = 0.1, phi = 1, P = 0.7
x6 <- c(0.28,1,0,0,0.12,0.7,0,0,0,0,0.2,1,0,0,0.2,0.7)
A6 <- matrix(x6, nrow=4, ncol=4, byrow=TRUE)
colnames(A6) <- c("slow juvenile", "slow adult", "fast juvenile", "fast adult")
rownames(A6) <- c("slow juvenile", "slow adult", "fast juvenile", "fast adult")
eigen.analysis(A6)
net.reproductive.rate(A6)
generation.time(A6)

tfA6 <- splitA(A6,r = c(1,3), c = c(2,4))
Tmat6 <- tfA6$T
fundamental.matrix(Tmat6)
matplot2(pop.projection(A6, c(1,1,1,1), 25)$stage.vectors, prop=TRUE, col= 10:14, lwd = 3,legend="topright")

##12.14.2015
##Figures below probably not for publication
##3D plot with g- x axis, sigma- y axis, lambda- z axis and grouped by phi
#1 Nov 2015 - I think I've reached a dead end here
#grp.phi <- gl(21,52500, ordered =TRUE)
#cols <- 1:21
#add axes labels 
#col <- heat.colors(52500, alpha = 0.5)
#col2 <- rainbow(52500, s =1, v=1, start = 0, end = 1, alpha = 1)
windows()
#plot1 <- scatterplot3d(eigen2$gbar, eigen2$sigma, eigen2$lam, angle=30,
              #color = col2,xlab = "Growth", ylab= expression(sigma), 
                       #zlab=expression(lambda))
cols <- rainbow(10, s= 1, v =1, start = 0, end = 1, alpha = 1)
cols2 <- heat.colors(100,alpha =1)
cols3 <- terrain.colors(20, alpha=1)
greycol <- grey.colors(10,start = 0, end = 0.8, alpha = NULL)
cols6 <- brewer.pal(5, "Spectral")
cols7 <- brewer.pal(8, "PRGn")
cols8 <- brewer.pal(8, "Blues")
cols12 <- cm.colors(8)
##I like 14
cols14 <- matlab.like2(10)
cols15 <- brewer.pal(10, "RdBu")
cols24 <- diverge_hcl(10)
cols25 <- rainbow_hcl(10)

#######

par(mar=c(1.3,1.3,1.3,1.3))
par(mfrow = c(1,1))

##test for colors 
## theta rotates left and right, postive theta left, negative right
## phi rotates up and down
## theta  = 55 is a good starting place in that plane
plot.test <- scatter3D(eigen2$gbar, eigen2$sigma, eigen2$phi, 
                   colvar=eigen2$lam, xlab="growth", ylab ="sigma",
                   zlab="phi", ticktype="detailed",
                   clab=c("lambda"), theta = 0, phi = -60)
## Damping ratio plot
plot.dr <- scatter3D(eigen2$gbar, eigen2$sigma, eigen2$phi, 
                       colvar=eigen2$DampR, xlab="growth", ylab ="sigma",
                       zlab="phi", ticktype="detailed",
                       clab=c("lambda"), theta = 0, phi = -60)

#top view KEEP 
#Nov 20 2015 
plot.1 <- scatter3D(eig.new$gbar, eig.new$sigma, eig.new$phi, 
                    colvar=eig.new$lam, xlab="Growth (g)", ylab ="sigma",
                    zlab="phi", ticktype="detailed",
                    clab=c("lambda"), theta = 55, phi = 55)
#Bottom view 
#Nov 20 2015 - not that great but keep
plot.2 <- scatter3D(eigen2$gbar, eigen2$sigma, eigen2$phi, 
                    colvar=eigen2$lam, xlab="Growth (g)", ylab ="sigma",
                    zlab="phi", ticktype="detailed",
                    clab=c("lambda"), theta = 45, phi = -90)

#lambda grey scale
plot1 <- scatter3D(eigen2$gbar, eigen2$sigma, eigen2$phi, 
            colvar=eigen2$lam, xlab="Growth (g)", ylab = expression(sigma),
            zlab=expression(phi), theta = 55, phi = 55, ticktype="detailed",
            clab=c("lambda"), col= greycol)
#add text
text(-1:0.2,-1:0.2,expression(phi))

#lambda default color
plot2 <- scatter3D(eigen2$gbar, eigen2$sigma, eigen2$phi, 
                   colvar=eigen2$lam, xlab="Growth (g)", ylab = expression(sigma),
                   zlab=expression(phi), ticktype="detailed",
                   clab=c("lambda"), phi = 80, theta = -20)
## add text 
text(-1:3,-1:3, labels = expression(sigma))

##change rotation, top view, GF says save 
plot3 <- scatter3D(eigen2$gbar, eigen2$sigma, eigen2$phi, 
                   colvar=eigen2$lam, xlab="Growth (g)", ylab = expression(sigma),
                   zlab=expression(phi), ticktype="detailed",
                   clab=c("lambda"), phi = 80, theta = 0)

#change rotation
plot3 <- scatter3D(eigen2$gbar, eigen2$sigma, eigen2$phi, 
                   colvar=eigen2$lam, xlab="Growth (g)", ylab = expression(sigma),
                   zlab=expression(phi), ticktype="detailed",
                   clab=c("lambda"), phi = 0, theta = 180)


#other colors
plot3 <- scatter3D(eigen2$gbar, eigen2$sigma, eigen2$phi, 
                   colvar=eigen2$lam, xlab="Growth (g)", ylab = expression(sigma),
                   zlab=expression(phi), ticktype="detailed",
                   clab=c("lambda"), phi = 40, theta = 50, col=cols15)


##Plot for R0 grey scale
plot4 <-scatter3D(eigen2$gbar, eigen2$sigma, eigen2$phi, 
                  colvar=eigen2$R0, xlab="Growth (g)", ylab = expression(sigma),
                  zlab=expression(phi),phi = 40, theta = 50, ticktype="detailed",
                  clab=c("R0"), col=greycol)

#default colors
plot5 <-scatter3D(eigen2$gbar, eigen2$sigma, eigen2$phi, 
                  colvar=eigen2$R0, xlab="Growth (g)", ylab = expression(sigma),
                  zlab=expression(phi),phi = 45, theta = -80, ticktype="detailed",
                  clab=c("R0"))

plot6 <-scatter3D(eigen2$gbar, eigen2$sigma, eigen2$phi, 
                  colvar=eigen2$R0, xlab="Growth (g)", ylab = expression(sigma),
                  zlab=expression(phi),phi = 90, theta = 45, ticktype="detailed",
                  clab=c("R0"))

##R0 color rainbow color
plot6 <-scatter3D(eigen2$gbar, eigen2$sigma, eigen2$phi, 
                  colvar=eigen2$R0, xlab="Growth (g)", ylab = expression(sigma),
                  zlab=expression(phi),phi = 0, theta = 50, ticktype="detailed",
                  clab=c("R0"), col=cols14)

##Plot for T
plot7 <-scatter3D(eigen2$gbar, eigen2$sigma, eigen2$phi, 
                   colvar=eigen2$time, xlab="Growth (g)", ylab = expression(sigma),
                   zlab=expression(phi),phi = 40, theta = 50, ticktype="detailed",
                   clab=c("T"), col=greycol)
##T color with default color
plot8 <-scatter3D(eigen2$gbar, eigen2$sigma, eigen2$phi, 
                  colvar=eigen2$time, xlab="Growth (g)", ylab = expression(sigma),
                  zlab=expression(phi),phi = 45, theta = 90, ticktype="detailed",
                  clab=c("T"))
#T with rainbow
plot9 <-scatter3D(eigen2$gbar, eigen2$sigma, eigen2$phi, 
                  colvar=eigen2$time, xlab="Growth (g)", ylab = expression(sigma),
                  zlab=expression(phi),phi = 0, theta = 50, ticktype="detailed",
                  clab=c("T"), col=cols14)


test2 <- plot3d(eigen2$gbar, eigen2$sigma, eigen2$phi)












