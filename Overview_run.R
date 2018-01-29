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
source("ec_nosrdiff_generalised_function_withr.R") 

##*** Setting up
# Number of discrete fitness levels? Resistance levels?
nfit = 5;   mres = 5;
# Array of distribution of fitness and resistance c[resistance, fitness]
M0 <- array(0,c(mres,nfit,10))
# M0 now needs initial condition - mostly "susceptible"
wildtype <- c(1,5) # most fit least resistant
susr <- wildtype[1]; susf <- wildtype[2]
M0[susr,susf,1] <- 1

# Initial acquisition distribution - bivariate here normal distribution with mean 0.5 and deviation 0.05
x <- seq(1/mres,1,1/mres) # seq(from = 0, to = 1, length.out = mres)
y <- seq(1/nfit,1,1/nfit) #seq(from = 0, to = 1, length.out = nfit)
f <- function(x, y) dmvnorm(cbind(x, y), mean = c(0.6, 0.6),sigma = diag(2)/20)
z <- outer(x, y, FUN = f); z <- z/sum(z) # Generate and normalise
dev.off(); persp(x, y, z, theta = -30, phi = 30, ticktype = "detailed")
acqdistn<-z

# Initial conditions
iniv<-c(98,1) # Mostly space
#iniv<-c(60,39,1)
dt=0.1
tsteps<-500*(1/dt)


####*** LOAD TO HERE ***######################################################################################################################################################################################################################################################################################################

####*** # (1) run with constant omega ***######################################################################################################################################################################################################################################################################################################
omega_M <- c(24,16,2,0.4)
submic_M <- c(1,1,1,1) # Linear decline
pref <- "Constant_omega_"
plot_diff_acd_output(acqdistn,"~/Documents/Hetero_res_and_f/plots",01, omega_M, submic_M, wildtype,pref)



