## Lifetime events - in this case, age at first reproduction
## From Caswell Chapter 5 
## Orca whale example 
library(popbio)
data(whale)
whale
generation.time(whale)
fundamental.matrix(whale, r =1, c = c(2,3))

b <- c(0.8093,0.8279,1,1)
tran <- matrix(c(0,0,0,0,0.9775,0.911,0,0,0,0.0736,0.9534,0,0,0,0.0452,0.9804),
               nrow=4,ncol=4,byrow=TRUE)

tc1 <- solve(diag(b))
tc2 <- tran%*%(diag(b))
tc <- tc1 %*% tc2
tc[3,3] <-0
tc[4,4] <-0
tc[4,3] <-0

identity <- matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),nrow=4,ncol=4,byrow=TRUE)

e <- matrix(c(1,1,1,1),nrow=4,ncol=1)


expect <- (t(e)) %*% (solve(identity-tc))
####IT WORKS! Finally!!!

#Identity matrix for later use
I <- matrix(c(1,0,0,1), nrow=2, ncol=2, byrow=TRUE)

#2 type (fast and slow)
I2 <- matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),
      nrow = 4, ncol = 4, byrow=TRUE)

#transition matrix for single type model 
trans <- matrix(ncol=2, nrow=2, byrow=TRUE,
          data = c(S*(1-gbar), 0, S*gbar, P))

#transition matrix for 2 type model
trans2 <-matrix(ncol=4, nrow=4, byrow=TRUE, data=c(
       S*(1-(gbar-sigma)), 0, 0, 0,
        S * (gbar-sigma), P, 0, 0,
       0, 0, S*(1 - (gbar + sigma)), 0,
       0, 0, S * (gbar + sigma), P))

#Part of new transition matrix P' - Eq 5.47 Caswell, pg 124
Tprime <- matrix(ncol=2, nrow=2, byrow=TRUE, 
          data = c(S*(1-gbar), 0, S*gbar, 0))

#Two type model
Tprime2 <- matrix(ncol=4, nrow=4, byrow=TRUE, data = c(
  S*(1-(gbar-sigma)), 0, 0, 0,
  S*(gbar-sigma), 0, 0, 0, 0,
  0, 0, S*(1-(gabr + sigma)), 0,
  0, 0, S*(gbar + sigma), 0))

#Part of P'
Mprime <- matrix(ncol=2, nrow=2, byrow=TRUE,
          data = c(1-S, 0, 0, 1))

#Two type model 
Mprime2 <- matrix(ncol=4, nrow=2, byrow=TRUE, 
          data = c(1-S, 0, 1-S, 0,
                   0, 1, 0, 1))

#Find B' and extract the second row  - probabilities of reproducing before death
Bprime <- Mprime %*% (solve(I - Tprime))
b2 <- Bprime[2,1:2]

#b2 two type
b2 <- Bprime2[2,1:4]

#New Markov Chain conditional on absorption in "reproduced before dying"
Tc <- (solve(diag(b2))) %*% trans %*% diag(b2)

#define e - a column vector of 1's
e <- matrix(c(1,1),nrow=2,ncol=1)

#2 type mode
e <- matrix(c(1,1,1,1), nrow = 1, ncol = 4)

mean.age <- t(e) %*% (solve(I - Tc))
  
### Single type growth model
### Finding the age at first reproduction 

#s = 0.4
#g = 0.6
#F = 1
#P = 0.7
pop <- matrix(c(0.16,1,0.24,0.7),nrow=2, ncol=2, byrow=TRUE)
colnames(pop) <- c("juvenile","adult")
rownames(pop) <- c("juvenile", "adult")

#splitA
spA <- splitA(pop, r = 1, c = 2)
trans <- spA$T
iden <- matrix(c(1,0,0,1), nrow=2, ncol=2, byrow=TRUE)

Tprime <- matrix(c(0.16,0,0.24,0), nrow=2, ncol=2, byrow=TRUE)
Mprime <- matrix(c(0.6,0,0,1), nrow=2, ncol=2, byrow=TRUE)
Pprime <- matrix(c(0.16,0,0,0,0.24,0,0,0,0.6,0,1,0,0,1,0,1), nrow=4, ncol=4,
                 byrow=TRUE)

Bprime <- Mprime %*% (solve(iden-Tprime))
b2 <- Bprime[2,1:2]

e <- matrix(c(1,1),nrow=2,ncol=1)
et <- t(e)

tc <- (solve(diag(b2))) %*% trans %*% diag(b2)
expect <- et %*% (solve(iden-tc))
expect

#### Functions for each step 
### Need to set up T', M' and I 

#This function is Eq (5.51) in Caswell, all entries all matrices
#input for m is mprime
#input for i is identity matrix 
#input for t is tprime 
bprime <- function(m,i,t){
  b <- m %*% (solve(i - t))
  return(b) 
}

#This object is the second row of Bprime, it's a row vector or 1 X 4 matrix
b2prime <- bprime[2,1:2]


#Function computes T(c), inputs and output are matrices
#the input for b is b2prime
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
  ma <- (t(e)) %*% solve((i-t))
  return(ma)
}

########
# What if I did the 2 type? See what happens
#s=0.4, g = 0.6, P = 0.7, F = 1, sigma = 0.05, phi = 0

pop.mat <- matrix(c(0.18,1,0,0,0.22,0.7,0,0,0,0,0.14,1,0,0,0.26,0.7),
                  nrow = 4, ncol=4, byrow=TRUE)

colnames(pop.mat) <- c("slow juvenile","slow adult", "fast juvenile","fast adult")
rownames(pop.mat) <- c("slow juvenile","slow adult", "fast juvenile","fast adult")

#splitA
spA2 <- splitA(pop.mat, r = c(1,3), c = c(2,4))
trans2 <- spA2$T
iden2 <- matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1), nrow=4, ncol=4, byrow=TRUE)

Tprime2 <- matrix(c(0.18,0,0,0,0.22,0,0,0,0,0,0.14,0,0,0,0.26,0),
                  nrow=4, ncol=4, byrow=TRUE)
Mprime2 <- matrix(c(0.6,0,0.6,0,0,1,0,1), nrow=2, ncol=4, byrow=TRUE)

Bprime2 <- Mprime2 %*% (solve(iden2-Tprime2))
b2.2 <- Bprime2[2,1:4]

e2 <- matrix(c(1,1,1,1),nrow=4,ncol=1)

tc2 <- solve(diag(b2.2)) %*% trans2 %*% diag(b2.2)
expect2 <- (t(e2)) %*% (solve(iden2 - tc2))




