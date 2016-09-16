# Figures 

library(ggplot2)
new.data <- read.csv("May3b.csv", header=TRUE)
names(new.data) <- c("X", "gamma","S","sigma","lam","eigen2","phi",
                     "r","R0","damp","time","ages","agef")
data.56 <- read.csv("May10.csv", header = TRUE)
names(data.56) <- c("X", "gamma","S","sigma","lam","eigen2","phi",
                    "r","R0","damp","time","ages","agef")

no.sig4 <- subset(new.data, sigma ==0)
no.sig5 <- subset(data.56, sigma == 0 & S == 0.5)
no.sig6 <- subset(data.56, sigma == 0 & S == 0.6)

no.sig.df <- rbind(no.sig4,no.sig5,no.sig6)
no.sig.df$S <- as.factor(no.sig.df$S)


#Lambda - monotypic population or where phi = 1
#Want a legend for each line 
ggplot() + geom_line(data=no.sig.df, aes(x=gamma, y=lam, 
                colour=S)) +
  scale_colour_manual(values=c("blue","red","darkgreen"))+
labs(x=expression(gamma), y =expression(lambda)) +
  theme(axis.text.y = element_text(size = 11)) +
  theme(axis.text.x = element_text(size = 11)) +
  theme(axis.title.x = element_text(face = "bold.italic", size = 18)) +
  theme(axis.title.y = element_text(face ="bold.italic", size = 18))
  




##Won't work because need to map for legend 
ggplot() + 
  geom_line(data=no.sig5, aes(x=gamma, y=lam), color = "green") +
  geom_line(data=no.sig6, aes(x =gamma, y=lam), color= "blue") +
  geom_line(data=no.sig4, aes(x=gamma, y=lam), color= "red") +
  labs(x=expression(gamma), y =expression(lambda)) 
  
