# 9 14 2016
# code from MeanAge.R, EMM3.R, and May31Bruce.Rmd

#Figures for Dissertation Chapter

library(ggplot2)



#S = 0.4 
#Read in data from output 

new.data <- read.csv("May3b.csv", header=TRUE)
names(new.data) <- c("X", "gamma","S","sigma","lam","eigen2","phi",
                     "r","R0","damp","time","ages","agef")

#S = 0.5 & 0.6 
#Read in data from output
data.56 <- read.csv("May10.csv", header = TRUE)
names(data.56) <- c("X", "gamma","S","sigma","lam","eigen2","phi",
                    "r","R0","damp","time","ages","agef")

# Subset data where sigma = 0 to see relationship between g and lambda, R0, T
# Using data.56 and new.data
no.sig4 <- subset(new.data, sigma == 0)
no.sig5 <- subset(data.56, sigma == 0 & S == 0.5)
no.sig6 <- subset(data.56, sigma == 0 & S == 0.6)

#data frame for all so I can have a legend on graphs 
all <- rbind(no.sig4,no.sig5,no.sig6)

## Figure 2, monotypic population

# Lambda
ggplot() + 
  geom_line(data=no.sig5, aes(x=gamma, y=lam), color = "green") +
  geom_line(data=no.sig6, aes(x =gamma, y=lam), color= "blue") +
  geom_line(data=no.sig4, aes(x=gamma, y=lam), color= "red") +
  labs(x=expression(gamma), y =expression(lambda)) 

ggplot() + 
  geom_line(data=no.sig5, aes(x=gamma, y=R0), color = "green") +
  geom_line(data=no.sig6, aes(x =gamma, y=R0), color= "blue") +
  geom_line(data=no.sig4, aes(x=gamma, y=R0), color= "red") +
  labs(x=expression(gamma), y = "R0") 

ggplot() +
  geom_line(data=no.sig5, aes(x=gamma, y=time), color = "green") +
  geom_line(data=no.sig6, aes(x =gamma, y=time), color= "blue") +
  geom_line(data=no.sig4, aes(x=gamma, y=time), color= "red") +
  labs(x=expression(gamma), y = "Generation time")
