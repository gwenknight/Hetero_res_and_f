#### OVERVIEW
## Runs diverse bacterial population 

##*** Libraries needed and baseline graphics
library(mvtnorm);library(plyr); library(ggplot2);library(reshape2);library(deSolve);library(grid);library(gtools); library(directlabels); library(mvtnorm)
theme_set(theme_gray(base_size = 24)); cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7");

##*** Locations
home<-"~/Documents/Hetero_res_and_f/"
setwd(home)
plots_home<-paste(home,"plots",sep="")

##*** Code needed
# Loads functions for generalised mean function and simulation model that uses the generalised function
source("ec_generalised_function_withr.R") 

##*** Setting up
# Number of discrete fitness levels? Resistance levels?
nfit = 5;   mres = 5;
# Array of distribution of fitness and resistance c[resistance, fitness]
M0 <- array(0,c(mres,nfit,10))

# Initial acquisition distribution - bivariate here normal distribution with mean 0.5 and deviation 0.05
x <- seq(1/mres,1,1/mres) # seq(from = 0, to = 1, length.out = mres)
y <- seq(1/nfit,1,1/nfit) #seq(from = 0, to = 1, length.out = nfit)
f <- function(x, y) dmvnorm(cbind(x, y), mean = c(0.6, 0.6),sigma = diag(2)/20)
z <- outer(x, y, FUN = f); z <- z/sum(z) # Generate and normalise
dev.off(); persp(x, y, z, theta = -30, phi = 30, ticktype = "detailed")
acqdistn<-z
plot(rowSums(z),type="l");plot(colSums(z),type="l") # Same as normal distribution
z<-as.data.frame(z)
rownames(z) <- seq(1/mres,1,1/mres);colnames(z) <- seq(1,nfit,1);
z2<-as.data.frame(melt(z)); z2$res<-seq(1/mres,1,1/mres); colnames(z2)<-c("fitness","value","res")
p<-ggplot(z2, aes(x=res, y=value, fill=factor(fitness))) + geom_bar(stat="identity",colour="black") + facet_grid(~fitness) 
p<-p + scale_x_continuous("Resistance level",breaks=c(0,0.2,0.4,0.6,0.8,1)) + scale_y_continuous("Proportion") + scale_fill_brewer("Fitness \nlevel",palette="Reds") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p
setwd(plots_home)
ggsave("acqdistn_06.pdf",width=14,height=10)

# Initial conditions
iniv<-c(98,1) # Mostly space
#iniv<-c(60,39,1)

####*** LOAD TO HERE ***######################################################################################################################################################################################################################################################################################################

####*** # (1) run with constant omega ***######################################################################################################################################################################################################################################################################################################






