### Finding the age at first reproduction 

pop <- matrix(c(0.16,1,0.24,0.7),nrow=2, ncol=2, byrow=TRUE)

Tprime <- matrix(c(0.16,0,0.24,0), nrow=2, ncol=2, byrow=TRUE)
Mprime <- matrix(c(0.6,0,0,1), nrow=2, ncol=2, byrow=TRUE)
Pprime <- matrix(c(0.16,0,0,0,0.24,0,0,0,0.6,0,1,0,0,1,0,1), nrow=4, ncol=4,
                 byrow=TRUE)

trans <- matrix(c(0.16,0,0.24,0.7), nrow=2, ncol=2, byrow=TRUE)

iden <- matrix(c(1,0,0,1), nrow=2, ncol=2, byrow=TRUE)

#Calculate the second part of 5.51 (I-T')^-1
bb <- iden-tc
bbb <- solve(bb)
expec <- exp(trans)%*%bbb


Bprime <- Mprime%*%second

b2 <- c(0.2857,1)
tc <- solve(diag(b2))%*%(trans%*%(diag(b2)))

#The term with matrix subtraction then inverse eq 5.58
sub <- solve(iden-tc)
expect <- exp(trans)%*%sub


####Orca whale
library(popbio)
data(whale)
whale


b <- c(0.8093,0.8279,1,1)
tran <- matrix(c(0,0,0,0,0.9775,0.911,0,0,0,0.0736,0.9534,0,0,0,0.0452,0.9804),
            nrow=4,ncol=4,byrow=TRUE)

tc1 <- solve(diag(b))
tc2 <- tran%*%(diag(b))
tc <- tc1 %*% tc2

identity <- matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),nrow=4,ncol=4,byrow=TRUE)

paran.58 <- solve((identity-tc))

