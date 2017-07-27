#October figures manuscript

#data subset code - MUST USE AFTER RUNNING SIMULATION! Comes from "eig.dat"
#subsetting data to construct figures, probably a much better way to do it 
dat1a <- subset(eig.dat, gbar=="0.1" & pc=="-1")
dat1b <- subset(eig.dat, gbar=="0.1" & pc=="-0.5")
dat1c <- subset(eig.dat, gbar=="0.1" & pc=="-0.3")
dat1d <- subset(eig.dat, gbar=="0.1" & pc=="-0.2")
dat1e <- subset(eig.dat, gbar=="0.1" & pc=="-0.1")
dat1f <- subset(eig.dat, gbar=="0.1" & pc=="0")
dat1g <- subset(eig.dat, gbar=="0.1" & pc=="0.1")
dat1h <- subset(eig.dat, gbar=="0.1" & pc=="0.3")
dat1i <- subset(eig.dat, gbar=="0.1" & pc=="0.4")
dat1j <- subset(eig.dat, gbar=="0.1" & pc=="0.5")
dat1k <- subset(eig.dat, gbar=="0.1" & pc=="0.8")
dat1l <- subset(eig.dat, gbar=="0.1" & pc=="1")


dat2a <- subset(eig.dat, gbar=="0.2" & pc=="-1")
dat2b <- subset(eig.dat, gbar=="0.2" & pc=="-0.5")
dat2c <- subset(eig.dat, gbar=="0.2" & pc=="-0.3")
dat2d <- subset(eig.dat, gbar=="0.2" & pc=="-0.2")
dat2e <- subset(eig.dat, gbar=="0.2" & pc=="-0.1")
dat2f <- subset(eig.dat, gbar=="0.2" & pc=="0")
dat2g <- subset(eig.dat, gbar=="0.2" & pc=="0.1")
dat2h <- subset(eig.dat, gbar=="0.2" & pc=="0.2")
dat2i <- subset(eig.dat, gbar=="0.2" & pc=="0.4")
dat2j <- subset(eig.dat, gbar=="0.2" & pc=="0.5")
dat2k <- subset(eig.dat, gbar=="0.2" & pc=="0.8")
dat2l <- subset(eig.dat, gbar=="0.2" & pc=="1")


dat3a <- subset(eig.dat, gbar=="0.3" & pc=="-1")
dat3b <- subset(eig.dat, gbar=="0.3" & pc=="-0.5")
dat3c <- subset(eig.dat, gbar=="0.3" & pc=="-0.3")
dat3d <- subset(eig.dat, gbar=="0.3" & pc=="-0.2")
dat3e <- subset(eig.dat, gbar=="0.3" & pc=="-0.1")
dat3f <- subset(eig.dat, gbar=="0.3" & pc=="0")
dat3g <- subset(eig.dat, gbar=="0.3" & pc=="0.1")
dat3h <- subset(eig.dat, gbar=="0.3" & pc=="0.3")
dat3i <- subset(eig.dat, gbar=="0.3" & pc=="0.4")
dat3j <- subset(eig.dat, gbar=="0.3" & pc=="0.5")
dat3k <- subset(eig.dat, gbar=="0.3" & pc=="0.8")
dat3l <- subset(eig.dat, gbar=="0.3" & pc=="1")


dat4a <- subset(eig.dat, gbar=="0.4" & pc=="-1")
dat4b <- subset(eig.dat, gbar=="0.4" & pc=="-0.5")
dat4c <- subset(eig.dat, gbar=="0.4" & pc=="-0.3")
dat4d <- subset(eig.dat, gbar=="0.4" & pc=="-0.2")
dat4e <- subset(eig.dat, gbar=="0.4" & pc=="-0.1")
dat4f <- subset(eig.dat, gbar=="0.4" & pc=="0")
dat4g <- subset(eig.dat, gbar=="0.4" & pc=="0.1")
dat4h <- subset(eig.dat, gbar=="0.4" & pc=="0.3")
dat4i <- subset(eig.dat, gbar=="0.4" & pc=="0.4")
dat4j <- subset(eig.dat, gbar=="0.4" & pc=="0.5")
dat4k <- subset(eig.dat, gbar=="0.4" & pc=="0.8")
dat4l <- subset(eig.dat, gbar=="0.4" & pc=="1")


dat5a <- subset(eig.dat, gbar=="0.5" & pc=="-1")
dat5b <- subset(eig.dat, gbar=="0.5" & pc=="-0.5")
dat5c <- subset(eig.dat, gbar=="0.5" & pc=="-0.3")
dat5d <- subset(eig.dat, gbar=="0.5" & pc=="-0.2")
dat5e <- subset(eig.dat, gbar=="0.5" & pc=="-0.1")
dat5f <- subset(eig.dat, gbar=="0.5" & pc=="0")
dat5g <- subset(eig.dat, gbar=="0.5" & pc=="0.1")
dat5h <- subset(eig.dat, gbar=="0.5" & pc=="0.3")
dat5i <- subset(eig.dat, gbar=="0.5" & pc=="0.4")
dat5j <- subset(eig.dat, gbar=="0.5" & pc=="0.5")
dat5k <- subset(eig.dat, gbar=="0.5" & pc=="0.8")
dat5l <- subset(eig.dat, gbar=="0.5" & pc=="1")


dat6a <- subset(eig.dat, gbar=="0.6" & pc=="-1")
dat6b <- subset(eig.dat, gbar=="0.6" & pc=="-0.5")
dat6c <- subset(eig.dat, gbar=="0.6" & pc=="-0.3")
dat6d <- subset(eig.dat, gbar=="0.6" & pc=="-0.2")
dat6e <- subset(eig.dat, gbar=="0.6" & pc=="-0.1")
dat6f <- subset(eig.dat, gbar=="0.6" & pc=="0")
dat6g <- subset(eig.dat, gbar=="0.6" & pc=="0.1")
dat6h <- subset(eig.dat, gbar=="0.6" & pc=="0.3")
dat6i <- subset(eig.dat, gbar=="0.6" & pc=="0.4")
dat6j <- subset(eig.dat, gbar=="0.6" & pc=="0.5")
dat6k <- subset(eig.dat, gbar=="0.6" & pc=="0.8")
dat6l <- subset(eig.dat, gbar=="0.6" & pc=="1")


dat7a <- subset(eig.dat, gbar=="0.7" & pc=="-1")
dat7b <- subset(eig.dat, gbar=="0.7" & pc=="-0.5")
dat7c <- subset(eig.dat, gbar=="0.7" & pc=="-0.3")
dat7d <- subset(eig.dat, gbar=="0.7" & pc=="-0.2")
dat7e <- subset(eig.dat, gbar=="0.7" & pc=="-0.1")
dat7f <- subset(eig.dat, gbar=="0.7" & pc=="0")
dat7g <- subset(eig.dat, gbar=="0.7" & pc=="0.1")
dat7h <- subset(eig.dat, gbar=="0.7" & pc=="0.3")
dat7i <- subset(eig.dat, gbar=="0.7" & pc=="0.4")
dat7j <- subset(eig.dat, gbar=="0.7" & pc=="0.5")
dat7k <- subset(eig.dat, gbar=="0.7" & pc=="0.8")
dat7l <- subset(eig.dat, gbar=="0.7" & pc=="1")


dat8a <- subset(eig.dat, gbar=="0.8" & pc=="-1")
dat8b <- subset(eig.dat, gbar=="0.8" & pc=="-0.5")
dat8c <- subset(eig.dat, gbar=="0.8" & pc=="-0.3")
dat8d <- subset(eig.dat, gbar=="0.8" & pc=="-0.2")
dat8e <- subset(eig.dat, gbar=="0.8" & pc=="-0.1")
dat8f <- subset(eig.dat, gbar=="0.8" & pc=="0")
dat8g <- subset(eig.dat, gbar=="0.8" & pc=="0.1")
dat8h <- subset(eig.dat, gbar=="0.8" & pc=="0.3")
dat8i <- subset(eig.dat, gbar=="0.8" & pc=="0.4")
dat8j <- subset(eig.dat, gbar=="0.8" & pc=="0.5")
dat8k <- subset(eig.dat, gbar=="0.8" & pc=="0.8")
dat8l <- subset(eig.dat, gbar=="0.8" & pc=="1")


dat9a <- subset(eig.dat, gbar=="0.9" & pc=="-1")
dat9b <- subset(eig.dat, gbar=="0.9" & pc=="-0.5")
dat9c <- subset(eig.dat, gbar=="0.9" & pc=="-0.3")
dat9d <- subset(eig.dat, gbar=="0.9" & pc=="-0.2")
dat9e <- subset(eig.dat, gbar=="0.9" & pc=="-0.1")
dat9f <- subset(eig.dat, gbar=="0.9" & pc=="0")
dat9g <- subset(eig.dat, gbar=="0.9" & pc=="0.1")
dat9h <- subset(eig.dat, gbar=="0.9" & pc=="0.3")
dat9i <- subset(eig.dat, gbar=="0.9" & pc=="0.4")
dat9j <- subset(eig.dat, gbar=="0.9" & pc=="0.5")
dat9k <- subset(eig.dat, gbar=="0.9" & pc=="0.8")
dat9l <- subset(eig.dat, gbar=="0.9" & pc=="1")

######
#December 9 2015
#Figures for Bruce  - see email exchange from 12/3-12/9
par(mar=c(4,4,4,4))
par(mfrow=c(2,2))
plot(dat1a$sigma, dat1a$R0, type = "l", lty = 1, xlim=c(0, 0.49),
     ylim = c(0,1.0), xlab = expression(sigma = sigma), ylab = "R0",
     main = expression("Phi = -1"))
lines(dat3a$sigma, dat3a$R0, lty = 2)
lines(dat5a$sigma, dat5a$R0, lty = 3)
lines(dat7a$sigma, dat7a$R0, lty = 4)
lines(dat9a$sigma, dat9a$R0, lty = 5)

plot(dat1b$sigma, dat1b$R0, type = "l", lty = 1, xlim=c(0, 0.49), 
     ylim = c(0,1.0), xlab = expression(sigma = sigma), ylab = "R0",
     main = expression("Phi = -0.5"))
lines(dat3b$sigma, dat3b$R0, lty = 2)
lines(dat5b$sigma, dat5b$R0, lty = 3)
lines(dat7b$sigma, dat7b$R0, lty = 4)
lines(dat9b$sigma, dat9b$R0, lty = 5)

plot(dat1f$sigma, dat1f$R0, type = "l", lty = 1, xlim=c(0, 0.49),
     ylim = c(0, 1.0), xlab = expression(sigma = sigma), ylab = "R0",
     main = expression("Phi = 0"))
lines(dat3f$sigma, dat3f$R0, lty = 2)
lines(dat5f$sigma, dat5f$R0, lty = 3)
lines(dat7f$sigma, dat7f$R0, lty = 4)
lines(dat9f$sigma, dat9f$R0, lty = 5)

plot(dat1h$sigma, dat1h$R0, type = "l", lty = 1, xlim=c(0, 0.49),
     ylim = c(0, 1.0), xlab = expression(sigma = sigma), ylab = "R0",
     main = expression("Phi = 0.3"))
lines(dat3h$sigma, dat3h$R0, lty = 2)
lines(dat5h$sigma, dat5h$R0, lty = 3)
lines(dat7h$sigma, dat7h$R0, lty = 4)
lines(dat9h$sigma, dat9h$R0, lty = 5)

legend("bottomright", title ="g", c("0.1","0.3","0.5","0.7","0.9"),lty= c(1,2,3,4,5), cex = 0.75)
########

par(mar=c(4,4,4,4))
par(mfrow=c(2,2))
plot(dat1i$sigma, dat1i$R0, type = "l", lty = 1, xlim=c(0, 0.49),
     ylim = c(0, 1.0), xlab = expression(sigma = sigma), ylab = "R0",
     main = expression("Phi = 0.4"))
lines(dat3i$sigma, dat3i$R0, lty = 2)
lines(dat5i$sigma, dat5i$R0, lty = 3)
lines(dat7i$sigma, dat7i$R0, lty = 4)
lines(dat9i$sigma, dat9i$R0, lty = 5)

plot(dat1j$sigma, dat1j$R0, type = "l", lty = 1, xlim=c(0, 0.49), 
     ylim = c(0,1.0), xlab = expression(sigma = sigma), ylab = "R0",
     main = expression("Phi = 0.5"))
lines(dat3j$sigma, dat3j$R0, lty = 2)
lines(dat5j$sigma, dat5j$R0, lty = 3)
lines(dat7j$sigma, dat7j$R0, lty = 4)
lines(dat9j$sigma, dat9j$R0, lty = 5)

plot(dat1k$sigma, dat1k$R0, type = "l", lty = 1, xlim=c(0, 0.49),
     ylim = c(0, 1.0), xlab = expression(sigma = sigma), ylab = "R0",
     main = expression("Phi = 0.8"))
lines(dat3k$sigma, dat3k$R0, lty = 2)
lines(dat5k$sigma, dat5k$R0, lty = 3)
lines(dat7k$sigma, dat7k$R0, lty = 4)
lines(dat9k$sigma, dat9k$R0, lty = 5)

plot(dat1l$sigma, dat1l$R0, type = "l", lty = 1, xlim=c(0, 0.49),
     ylim = c(0, 1.0), xlab = expression(sigma = sigma), ylab = "R0",
     main = expression("Phi = 1"))
lines(dat3l$sigma, dat3l$R0, lty = 2)
lines(dat5l$sigma, dat5l$R0, lty = 3)
lines(dat7l$sigma, dat7l$R0, lty = 4)
lines(dat9l$sigma, dat9l$R0, lty = 5)

legend("bottomright", title ="g", c("0.1","0.3","0.5","0.7","0.9"),lty= c(1,2,3,4,5), cex = 0.75)

################################
##T
par(mar=c(4,4,4,4))
par(mfrow=c(2,2))
plot(dat1a$sigma, dat1a$time, type = "l", lty = 1, xlim=c(0, 0.49),
     ylim = c(3,9), xlab = expression(sigma = sigma), ylab = "T",
     main = expression("Phi = -1"))
lines(dat3a$sigma, dat3a$time, lty = 2)
lines(dat5a$sigma, dat5a$time, lty = 3)
lines(dat7a$sigma, dat7a$time, lty = 4)
lines(dat9a$sigma, dat9a$time, lty = 5)

plot(dat1b$sigma, dat1b$time, type = "l", lty = 1, xlim=c(0, 0.49), 
     ylim = c(3,9), xlab = expression(sigma = sigma), ylab = "T",
     main = expression("Phi = -0.5"))
lines(dat3b$sigma, dat3b$time, lty = 2)
lines(dat5b$sigma, dat5b$time, lty = 3)
lines(dat7b$sigma, dat7b$time, lty = 4)
lines(dat9b$sigma, dat9b$time, lty = 5)

plot(dat1f$sigma, dat1f$time, type = "l", lty = 1, xlim=c(0, 0.49),
     ylim = c(3,9), xlab = expression(sigma = sigma), ylab = "T",
     main = expression("Phi = 0"))
lines(dat3f$sigma, dat3f$time, lty = 2)
lines(dat5f$sigma, dat5f$time, lty = 3)
lines(dat7f$sigma, dat7f$time, lty = 4)
lines(dat9f$sigma, dat9f$time, lty = 5)

plot(dat1h$sigma, dat1h$time, type = "l", lty = 1, xlim=c(0, 0.49),
     ylim = c(3,9), xlab = expression(sigma = sigma), ylab = "T",
     main = expression("Phi = 0.3"))
lines(dat3h$sigma, dat3h$time, lty = 2)
lines(dat5h$sigma, dat5h$time, lty = 3)
lines(dat7h$sigma, dat7h$time, lty = 4)
lines(dat9h$sigma, dat9h$time, lty = 5)

legend("topright", title ="g", c("0.1","0.3","0.5","0.7","0.9"),lty= c(1,2,3,4,5), cex = 0.75)



#####
par(mar=c(4,4,4,4))
par(mfrow=c(2,2))
plot(dat1i$sigma, dat1i$time, type = "l", lty = 1, xlim=c(0, 0.49),
     ylim = c(3,8), xlab = expression(sigma = sigma), ylab = "T",
     main = expression("Phi = 0.4"))
lines(dat3i$sigma, dat3i$time, lty = 2)
lines(dat5i$sigma, dat5i$time, lty = 3)
lines(dat7i$sigma, dat7i$time, lty = 4)
lines(dat9i$sigma, dat9i$time, lty = 5)

plot(dat1j$sigma, dat1j$time, type = "l", lty = 1, xlim=c(0, 0.49), 
     ylim = c(3,8), xlab = expression(sigma = sigma), ylab = "T",
     main = expression("Phi = 0.5"))
lines(dat3j$sigma, dat3j$time, lty = 2)
lines(dat5j$sigma, dat5j$time, lty = 3)
lines(dat7j$sigma, dat7j$time, lty = 4)
lines(dat9j$sigma, dat9j$time, lty = 5)

plot(dat1k$sigma, dat1k$time, type = "l", lty = 1, xlim=c(0, 0.49),
     ylim = c(3,8), xlab = expression(sigma = sigma), ylab = "T",
     main = expression("Phi = 0.8"))
lines(dat3k$sigma, dat3k$time, lty = 2)
lines(dat5k$sigma, dat5k$time, lty = 3)
lines(dat7k$sigma, dat7k$time, lty = 4)
lines(dat9k$sigma, dat9k$time, lty = 5)

plot(dat1l$sigma, dat1l$time, type = "l", lty = 1, xlim=c(0, 0.49),
     ylim = c(3,8), xlab = expression(sigma = sigma), ylab = "T",
     main = expression("Phi = 1"))
lines(dat3l$sigma, dat3l$time, lty = 2)
lines(dat5l$sigma, dat5l$time, lty = 3)
lines(dat7l$sigma, dat7l$time, lty = 4)
lines(dat9l$sigma, dat9l$time, lty = 5)

legend("topright", title ="g", c("0.1","0.3","0.5","0.7","0.9"),lty= c(1,2,3,4,5), cex = 0.75)





















##Example code for figures 
#g = 0.5, 6 phis
par(mar=c(4,4,2,2))
plot(dat5k$sigma, dat5k$eigen, type = "l", lty=1, xlim=c(0, 0.49), ylim=c(1, 1.25), xlab = expression(sigma = sigma), ylab = expression(lambda= lambda), main = expression("g = 0.5"))
lines(dat5j$sigma, dat5j$eigen, lty =2)
lines(dat5i$sigma, dat5i$eigen, lty =3)
lines(dat5h$sigma, dat5h$eigen, lty = 4)
lines(dat5g$sigma, dat5g$eigen, lty =5)
lines(dat5e$sigma, dat5e$eigen, lty = 6)
legend("bottomleft", title ="phi", c("0.8","0.5","0.4","0.3","0.1","-0.1"),lty= c(1,2,3,4,5,6), cex = 0.75)

#R0
## Still example code 
#g= 0.5
par(mar=c(4,4,2,2))
plot(dat5k$sigma, dat5k$R0, type = "l", lty=1, xlim=c(0, 0.49), ylim=c(1.0,2.4), xlab = expression(sigma = sigma), ylab = "R0", main = expression("g = 0.5"))
lines(dat5j$sigma, dat5j$R0, lty=2)
lines(dat5i$sigma, dat5i$R0, lty=2)
lines(dat5h$sigma, dat5h$R0, lty=4)
lines(dat5g$sigma, dat5g$R0, lty=5)
lines(dat5e$sigma, dat5e$R0, lty=6)
legend("bottomleft", title ="phi", c("0.8","0.5","0.4","0.3","0.1","-0.1"),lty= c(1,2,3,4,5,6), cex = 0.75)

## Still example code 
#T, the classic T, called "time" in dataframes 
#g= 0.5
par(mar=c(4,4,2,2))
plot(dat5k$sigma, dat5k$time, type = "l", lty=1, xlim=c(0, 0.49), ylim=c(3.5,5), xlab = expression(sigma = sigma), ylab = "T", main = expression("g = 0.5"))
lines(dat5j$sigma, dat5j$time, lty=2)
lines(dat5i$sigma, dat5i$time, lty=2)
lines(dat5h$sigma, dat5h$time, lty=4)
lines(dat5g$sigma, dat5g$time, lty=5)
lines(dat5e$sigma, dat5e$time, lty=6)
legend("bottomleft", title ="phi", c("0.8","0.5","0.4","0.3","0.1","-0.1"),lty= c(1,2,3,4,5,6), cex = 0.75)

## Still example code 
#Bienvenu's T
#g= 0.5
par(mar=c(4,4,2,2))
plot(dat5k$sigma, dat5k$BVT, type = "l", lty=1, xlim=c(0, 0.49), ylim=c(3,5), xlab = expression(sigma = sigma), ylab = "Bienvenu's T", main = expression("g = 0.5"))
lines(dat5j$sigma, dat5j$BVT, lty=2)
lines(dat5i$sigma, dat5i$BVT, lty=2)
lines(dat5h$sigma, dat5h$BVT, lty=4)
lines(dat5g$sigma, dat5g$BVT, lty=5)
lines(dat5e$sigma, dat5e$BVT, lty=6)
legend("bottomleft", title ="phi", c("0.8","0.5","0.4","0.3","0.1","-0.1"),lty= c(1,2,3,4,5,6), cex = 0.75)


#####################
# Code below used on most recent manuscript draft
# Sent to Gordon on Wednesday, October 14, 2015


#figure for lambda and sigma across different g's and phi's
#combine figures 2-7 as in current manuscript
par(mfrow=c(2,1))
par(mar=c(4 ,4, 1.5, 1.5))
plot(dat1l$sigma, dat1l$eigen, type = "l",xlim=c(0, 0.49), ylim=c(0.9, 1.25),xlab="", ylab=expression(lambda = lambda))
lines(dat2l$sigma, dat2l$eigen)
lines(dat3l$sigma, dat3l$eigen)
lines(dat4l$sigma, dat4l$eigen)
lines(dat5l$sigma, dat5l$eigen)
lines(dat6l$sigma, dat6l$eigen)
lines(dat7l$sigma, dat7l$eigen)
lines(dat8l$sigma, dat8l$eigen)
lines(dat9l$sigma, dat9l$eigen)

plot(dat1f$sigma, dat1f$eigen, type = "l",xlim=c(0, 0.49), ylim=c(0.9, 1.25),xlab=expression(sigma = sigma), ylab=expression(lambda = lambda))
lines(dat2f$sigma, dat2f$eigen)
lines(dat3f$sigma, dat3f$eigen)
lines(dat4f$sigma, dat4f$eigen)
lines(dat5f$sigma, dat5f$eigen)
lines(dat6f$sigma, dat6f$eigen)
lines(dat7f$sigma, dat7f$eigen)
lines(dat8f$sigma, dat8f$eigen)
lines(dat9f$sigma, dat9f$eigen)


#lines below till next comment: Figure of sigma vs. lambda across g's over 4 negative phis
par(mfrow=c(2,2))
par(mar=c(4,4,2,2))
plot(dat1a$sigma, dat1a$eigen, type = "l", lty= 1, xlim=c(0, 0.49), ylim=c(0.75, 1.25),xlab='', ylab="",main = expression("phi = -1"))
lines(dat2a$sigma, dat2a$eigen, lty =1)
lines(dat3a$sigma, dat3a$eigen, lty=1)
lines(dat4a$sigma, dat4a$eigen, lty=1)
lines(dat5a$sigma, dat5a$eigen, lty=1)
lines(dat6a$sigma, dat6a$eigen, lty=1)
lines(dat7a$sigma, dat7a$eigen, lty=1)
lines(dat8a$eigen, dat8a$eigen, lty=1)
lines(dat9a$eigen, dat9a$eigen, lty =1)

par(mar=c(4,4,2,2))
plot(dat1b$sigma, dat1b$eigen, type="l", lty = 1, xlim=c(0, 0.49), ylim=c(0.75, 1.25), xlab="", ylab = '',main = expression("phi = -0.5"))
lines(dat2b$sigma, dat2b$eigen, lty=1)
lines(dat3b$sigma, dat3b$eigen, lty=1)
lines(dat4b$sigma, dat4b$eigen, lty=1)
lines(dat5b$sigma, dat5b$eigen, lty=1)
lines(dat6b$sigma, dat6b$eigen, lty=1)
lines(dat7b$sigma, dat7b$eigen, lty=1)
lines(dat8b$eigen, dat8b$eigen, lty=1)
lines(dat9b$eigen, dat9b$eigen, lty =1)

par(mar=c(4,4,2,2))
plot(dat1c$sigma, dat1c$eigen, type="l", lty = 1, xlim=c(0, 0.49), ylim=c(0.75, 1.25), xlab=expression(sigma = sigma), ylab = expression(lambda= lambda),main = expression("phi = -0.3"))
lines(dat2c$sigma, dat2c$eigen, lty=1)
lines(dat3c$sigma, dat3c$eigen, lty=1)
lines(dat4c$sigma, dat4c$eigen, lty=1)
lines(dat5c$sigma, dat5c$eigen, lty=1)
lines(dat6c$sigma, dat6c$eigen, lty=1)
lines(dat7c$sigma, dat7c$eigen, lty=1)
lines(dat8c$eigen, dat8c$eigen, lty=1)
lines(dat9c$eigen, dat9c$eigen, lty =1)

par(mar=c(4,4,2,2))
plot(dat1d$sigma, dat1d$eigen, type="l", lty = 1, xlim=c(0, 0.49), ylim=c(0.75, 1.25), xlab="", ylab = '',main = expression("phi = -0.2"))
lines(dat2d$sigma, dat2d$eigen, lty=1)
lines(dat3d$sigma, dat3d$eigen, lty=1)
lines(dat4d$sigma, dat4d$eigen, lty=1)
lines(dat5d$sigma, dat5d$eigen, lty=1)
lines(dat6d$sigma, dat6d$eigen, lty=1)
lines(dat7d$sigma, dat7d$eigen, lty=1)
lines(dat8d$eigen, dat8d$eigen, lty=1)
lines(dat9d$eigen, dat9d$eigen, lty =1)

#lines : Figure of sigma vs. lambda across g's over 4 positive phis
par(mfrow=c(2,2))
par(mar=c(4,4,2,2))
plot(dat1h$sigma, dat1h$eigen, type = "l", lty= 1, xlim=c(0, 0.49), ylim=c(0.9, 1.3),xlab="", ylab="",main = expression("phi = 0.3"))
lines(dat2h$sigma, dat2h$eigen, lty =1)
lines(dat3h$sigma, dat3h$eigen, lty=1)
lines(dat4h$sigma, dat4h$eigen, lty=1)
lines(dat5h$sigma, dat5h$eigen, lty=1)
lines(dat6h$sigma, dat6h$eigen, lty=1)
lines(dat7h$sigma, dat7h$eigen, lty=1)
lines(dat8h$eigen, dat8h$eigen, lty=1)
lines(dat9h$eigen, dat9h$eigen, lty =1)

par(mar=c(4,4,2,2))
plot(dat1i$sigma, dat1i$eigen, type = "l", lty= 1, xlim=c(0, 0.49), ylim=c(0.9, 1.3),xlab="", ylab="",main = expression("phi = 0.4"))
lines(dat2i$sigma, dat2i$eigen, lty =1)
lines(dat3i$sigma, dat3i$eigen, lty=1)
lines(dat4i$sigma, dat4i$eigen, lty=1)
lines(dat5i$sigma, dat5i$eigen, lty=1)
lines(dat6i$sigma, dat6i$eigen, lty=1)
lines(dat7i$sigma, dat7i$eigen, lty=1)
lines(dat8i$eigen, dat8i$eigen, lty=1)
lines(dat9i$eigen, dat9i$eigen, lty =1)

par(mar=c(4,4,2,2))
plot(dat1j$sigma, dat1j$eigen, type = "l", lty= 1, xlim=c(0, 0.49), ylim=c(0.75, 1.3),xlab=expression(sigma = sigma), ylab=expression(lambda = lambda),main = expression("phi = 0.5"))
lines(dat2j$sigma, dat2j$eigen, lty =1)
lines(dat3j$sigma, dat3j$eigen, lty=1)
lines(dat4j$sigma, dat4j$eigen, lty=1)
lines(dat5j$sigma, dat5j$eigen, lty=1)
lines(dat6j$sigma, dat6j$eigen, lty=1)
lines(dat7j$sigma, dat7j$eigen, lty=1)
lines(dat8j$eigen, dat8j$eigen, lty=1)
lines(dat9j$eigen, dat9j$eigen, lty =1)

par(mar=c(4,4,2,2))
plot(dat1k$sigma, dat1k$eigen, type = "l", lty= 1, xlim=c(0, 0.49), ylim=c(0.75, 1.3),xlab='', ylab='',main = expression("phi = 0.8"))
lines(dat2k$sigma, dat2k$eigen, lty =1)
lines(dat3k$sigma, dat3k$eigen, lty=1)
lines(dat4k$sigma, dat4k$eigen, lty=1)
lines(dat5k$sigma, dat5k$eigen, lty=1)
lines(dat6k$sigma, dat6k$eigen, lty=1)
lines(dat7k$sigma, dat7k$eigen, lty=1)
lines(dat8k$eigen, dat8k$eigen, lty=1)
lines(dat9k$eigen, dat9k$eigen, lty =1)

####Generation Time
plot(dat1l$sigma, dat1l$time, type = "l", xlim=c(0, 0.49), ylim = c(1,6.5), xlab=expression(sigma=sigma), ylab=expression(lambda=lambda))
lines(dat2l$sigma, dat2l$time, lty =1)
lines(dat3l$sigma, dat3l$time, lty=1)
lines(dat4l$sigma, dat4l$time, lty=1)
lines(dat5l$sigma, dat5l$time, lty=1)
lines(dat6l$sigma, dat6l$time, lty=1)
lines(dat7l$sigma, dat7l$time, lty=1)
lines(dat8l$sigma, dat8l$time, lty=1)
lines(dat9l$sigma, dat9l$time, lty=1)
lines(dat5l$sigma, dat5l$BVT, lty=2)
