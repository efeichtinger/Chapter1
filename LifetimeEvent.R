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
et <- t(e)

expect <- et %*% (solve(identity-tc))
####IT WORKS! Finally!!!



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


