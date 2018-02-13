#### OVERVIEW
## Runs diverse bacterial population 

##*** Libraries needed and baseline graphics
library(mvtnorm);library(plyr); library(ggplot2);library(reshape2);library(deSolve);library(grid);library(gtools); library(directlabels); library(mvtnorm)
theme_set(theme_gray(base_size = 12)); cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7");

##*** Locations
home<-"~/Documents/Hetero_res_and_f/"
setwd(home)
plots_home<-paste(home,"plots",sep="")

##*** Code needed
# Loads functions for generalised mean function and simulation model that uses the generalised function
source("ec_nosrdiff_generalised_function_withr.R") 

##*** Setting up
# Number of discrete fitness levels? Resistance levels?
nfit = 10;   mres = 10;
# Array of distribution of fitness and resistance c[resistance, fitness]
M0 <- array(0,c(mres,nfit,10))
# M0 now needs initial condition - mostly "susceptible"
wildtype <- c(1,5) # most fit least resistant
susr <- wildtype[1]; susf <- wildtype[2]
M0[susr,susf,1] <- 1

# Initial acquisition distribution - bivariate here normal distribution with mean 0.6 and deviation 0.05
x <- seq(1/mres,1,1/mres) # seq(from = 0, to = 1, length.out = mres)
y <- seq(1/nfit,1,1/nfit) #seq(from = 0, to = 1, length.out = nfit)
f <- function(x, y) dmvnorm(cbind(x, y), mean = c(0.6, 0.6),sigma = diag(2)/20)
z <- outer(x, y, FUN = f); z <- z/sum(z) # Generate and normalise
dev.off(); persp(x, y, z, theta = -30, phi = 30, ticktype = "detailed")
acqdistn<-z

# Initial conditions
iniv<-c(80,20) # Mostly space
#iniv<-c(60,39,1)
dt=0.1
tsteps<-500*(1/dt)
kk <- 500


####*** LOAD TO HERE ***######################################################################################################################################################################################################################################################################################################

####*** # (1) run with constant omega ***######################################################################################################################################################################################################################################################################################################
omega_M <- c(24,16,2,0.4)
submic_M <- c(1,1,1,1) # Linear decline
pref <- "Constant_omega_"
plot_diff_acd_output(acqdistn,"~/Documents/Hetero_res_and_f/plots",01, omega_M, submic_M, wildtype,pref)

####*** # (2) run with timevarying omega ***######################################################################################################################################################################################################################################################################################################
omega_M <- rbind(c(matrix(16,1,20),matrix(0,1,tsteps - 20)),c(16,16*exp(-0.1*(2:tsteps))),
                 c(matrix(16,1,20),matrix(0.4,1,20),matrix(16,1,20),matrix(0,1,tsteps - 60)), c(matrix(0,1,20),matrix(16,1,tsteps - 20)))
submic_M <- c(1,1,1,1) # Linear decline
pref <- "TV_omega_"
plot_diff_acd_output(acqdistn,"~/Documents/Hetero_res_and_f/plots",01, omega_M, submic_M, wildtype,pref)

####*** # (3) run with stepped omega ***######################################################################################################################################################################################################################################################################################################
omega_M <- rbind(c(matrix(16,1,10*(1/dt)),matrix(0,1,tsteps - 10*(1/dt))),
                 c(matrix(16,1,20*(1/dt)),matrix(0,1,tsteps - 20*(1/dt))),
                 c(matrix(16,1,30*(1/dt)),matrix(0,1,tsteps - 30*(1/dt))),
                 c(matrix(16,1,40*(1/dt)),matrix(0,1,tsteps - 40*(1/dt))))
submic_M <- c(1,1,1,1) # Linear decline
pref <- "Stepped_omega_"
plot_diff_acd_output(acqdistn,"~/Documents/Hetero_res_and_f/plots",01, omega_M, submic_M, wildtype,pref)

####*** # (4) run with high short or long low omega ***######################################################################################################################################################################################################################################################################################################
omega_M <- rbind(c(matrix(30,1,10*(1/dt)),matrix(0,1,tsteps - 10*(1/dt))),
                 c(matrix(30,1,20*(1/dt)),matrix(0,1,tsteps - 20*(1/dt))),
                 c(matrix(0.4,1,50*(1/dt)),matrix(0,1,tsteps - 50*(1/dt))),
                 c(matrix(0.4,1,100*(1/dt)),matrix(0,1,tsteps - 100*(1/dt))))
submic_M <- c(1,1,1,1) # Linear decline
pref <- "HSoLL_"
plot_diff_acd_output(acqdistn,"~/Documents/Hetero_res_and_f/plots",01, omega_M, submic_M, wildtype,pref)

####*** # (5) run with PULSING high omega ***######################################################################################################################################################################################################################################################################################################
exp_dec <- -0.01
gap_doses <- 120
dose_h <- 16 
fn_pulse <- dose_h*exp(exp_dec*(0:gap_doses))
omega_M <- rbind(c(dose_h*exp(exp_dec*(0:(tsteps-1))) ),
                 c(fn_pulse,fn_pulse,                                     dose_h*exp(exp_dec*(0:(tsteps-2*(gap_doses+1)-1)))  ),
                 c(fn_pulse,fn_pulse,fn_pulse,fn_pulse,fn_pulse,          dose_h*exp(exp_dec*(0:(tsteps-5*(gap_doses+1)-1)))  ),
                 c(fn_pulse,fn_pulse,fn_pulse,fn_pulse,fn_pulse,
                   fn_pulse,fn_pulse,fn_pulse,fn_pulse,fn_pulse,fn_pulse, dose_h*exp(exp_dec*(0:(tsteps-11*(gap_doses+1)-1))) ) )
                 
pref <- "Pulse_high_"
plot_diff_acd_output(acqdistn,"~/Documents/Hetero_res_and_f/plots",01, omega_M, submic_M, wildtype,pref)


####*** # (6) run with PULSING low omega ***######################################################################################################################################################################################################################################################################################################
exp_dec <- -0.01
gap_doses <- 120
dose_h <- 2
fn_pulse <- dose_h*exp(exp_dec*(0:gap_doses))
omega_M <- rbind(c(dose_h*exp(exp_dec*(0:(tsteps-1))) ),
                 c(fn_pulse,fn_pulse,                                     dose_h*exp(exp_dec*(0:(tsteps-2*(gap_doses+1)-1)))  ),
                 c(fn_pulse,fn_pulse,fn_pulse,fn_pulse,fn_pulse,          dose_h*exp(exp_dec*(0:(tsteps-5*(gap_doses+1)-1)))  ),
                 c(fn_pulse,fn_pulse,fn_pulse,fn_pulse,fn_pulse,
                   fn_pulse,fn_pulse,fn_pulse,fn_pulse,fn_pulse,fn_pulse, dose_h*exp(exp_dec*(0:(tsteps-11*(gap_doses+1)-1))) ) )

pref <- "Pulse_low_"
plot_diff_acd_output(acqdistn,"~/Documents/Hetero_res_and_f/plots",01, omega_M, submic_M, wildtype,pref)

####*** # (7) run with high-low or low-high ***######################################################################################################################################################################################################################################################################################################
omega_M <- rbind(c(matrix(16,1,20*(1/dt)), matrix(0.4,1,20*(1/dt)),matrix(0,1,tsteps - 40*(1/dt))),
                 c(matrix(0.4,1,20*(1/dt)),matrix(16,1,20*(1/dt)), matrix(0,1,tsteps - 40*(1/dt))),
                 c(matrix(16,1,50*(1/dt)), matrix(0.4,1,50*(1/dt)),matrix(0,1,tsteps - 100*(1/dt))),
                 c(matrix(0.4,1,50*(1/dt)),matrix(16,1,50*(1/dt)), matrix(0,1,tsteps - 100*(1/dt))))
submic_M <- c(1,1,1,1) # Linear decline
pref <- "Highlo_lohigh_"
plot_diff_acd_output(acqdistn,"~/Documents/Hetero_res_and_f/plots",01, omega_M, submic_M, wildtype,pref)

####******######################################################################################################################################################################################################################################################################################################
####******######################################################################################################################################################################################################################################################################################################
####*** BIO ACQUISITION DISTBN ***######################################################################################################################################################################################################################################################################################################
####******######################################################################################################################################################################################################################################################################################################
##*** Setting up
# Number of discrete fitness levels? Resistance levels?
nfit = 10;   mres = 10;
# Array of distribution of fitness and resistance c[resistance, fitness]
M0 <- array(0,c(mres,nfit,10))
# M0 now needs initial condition - mostly "susceptible"
wildtype <- c(1,8) # most fit least resistant
susr <- wildtype[1]; susf <- wildtype[2]
M0[susr,susf,1] <- 1

# most low level acquisition distribution 
x <- seq(1/mres,1,1/mres)
y <- seq(1/nfit,1,1/nfit)
f <- function(x, y) dmvnorm(cbind(x, y), mean = c(0.05, 0.8),sigma = diag(2)/10)
z <- outer(x, y, FUN = f); z <- z/sum(z) # Generate and normalise
dev.off(); persp(x, y, z, theta = 60, phi = 30, ticktype = "detailed")
acqdistn<-z


###### MIC SHAPE
submic_M <- c(1,1,1,1) # Linear decline


####*** # (0) run with zero omega ***######################################################################################################################################################################################################################################################################################################
Sv20<-ec_srs_funcf_mean_varsr(tsteps,home, c(kk,0),iniv,M0,acqdistn,dt,1)
# What happens? 
mm20<-c(); ll<-dim(Sv20$M)[3];
ss<-seq(0,ll,1/dt) # Don't want to grab all 
for(i in 2:length(ss)){
  mm220<-as.data.frame(melt(Sv20$M[,,ss[i]])); mm220$tstep=ss[i]*dt; mm20<-rbind(mm20,mm220);
} 
colnames(mm20)<-c("x","y","z","t"); 
setwd(plots)
# Grab a subset
sub<-c(1/dt,100,250,500,750,1000,1500,2000,2500,seq(3000,4001,500))*dt
w<-which(mm20[,"t"] %in% sub)
# plots
p1<-ggplot(mm20[w,],aes(x,y,fill=z))  + facet_wrap( ~ t, ncol=3) + ggtitle(paste("w = ", 0,sep="")) + 
  scale_fill_gradient("Proportion", limits=c(0,1),low="white", high="red") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_tile() + 
  scale_y_continuous(breaks=c(1,nfit),"Relative fitness levels",labels=c("Least","Most")) + 
  scale_x_continuous(breaks=c(mres,1),"Resistance levels",labels=c("Most","Least"))
p1
ggsave(paste("norsd_Array_w=",0,"_06.pdf",sep=""))

# plots only proportion with "resistance" - not background initial resistance
M20n <- Sv20$M;
M20n[susr,susf,] <- 0; 
for (i in 1:dim(Sv20$M)[3]){ M20n[,,i] <- M20n[,,i]/sum(M20n[,,i])} # re-normalise

# What happens? 
mm20<-c() 
ll<-dim(Sv20$M)[3];
ss<-seq(0,ll,1/dt) # Don't want to grab all 
for(i in 2:length(ss)){
  mm220<-as.data.frame(melt(M20n[,,ss[i]])); 
  mm220$tstep=ss[i]*dt; 
  mm20<-rbind(mm20,mm220);
} 
colnames(mm20)<-c("x","y","z","t"); colnames(mm15)<-c("x","y","z","t"); colnames(mm10)<-c("x","y","z","t") ;colnames(mm05)<-c("x","y","z","t")

# plots
p1<-ggplot(mm20[w,],aes(x,y,fill=z))  + facet_wrap( ~ t, ncol=3) + ggtitle(paste("w = ", 0,sep=""))
p1<-p1 + scale_fill_gradient("Proportion", limits=c(0,1),low="white", high="red",guide = FALSE)
p1<-p1 + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p1<-p1 + geom_tile() + scale_y_continuous(breaks=c(1,nfit),"Relative fitness levels",labels=c("Least","Most")) + scale_x_continuous(breaks=c(mres,1),"Resistance levels",labels=c("Most","Least"))
p1
ggsave(paste("norsdR_Array_w=",0,"_06.pdf",sep=""))

pdf("B_no_omega.pdf")
plot(Sv20$B,type="l", ylim = c(0,100), ylab = c("Percentage"), xlab = "Time")
dev.off()


####*** # (1) run with constant omega ***######################################################################################################################################################################################################################################################################################################
omega_M <- c(24,16,2,0.4)

pref <- "Constant_omega_"
plot_diff_acd_output(acqdistn,"~/Documents/Hetero_res_and_f/plots",02, omega_M, submic_M, wildtype,pref)

####*** # (2) run with timevarying omega ***######################################################################################################################################################################################################################################################################################################
omega_M <- rbind(c(matrix(16,1,20),matrix(0,1,tsteps - 20)),c(16,16*exp(-0.1*(2:tsteps))),
                 c(matrix(16,1,20),matrix(0.4,1,20),matrix(16,1,20),matrix(0,1,tsteps - 60)), c(matrix(0,1,20),matrix(16,1,tsteps - 20)))

pref <- "TV_omega_"
plot_diff_acd_output(acqdistn,"~/Documents/Hetero_res_and_f/plots",02, omega_M, submic_M, wildtype,pref)

####*** # (3) run with stepped omega ***######################################################################################################################################################################################################################################################################################################
omega_M <- rbind(c(matrix(16,1,10*(1/dt)),matrix(0,1,tsteps - 10*(1/dt))),
                 c(matrix(16,1,20*(1/dt)),matrix(0,1,tsteps - 20*(1/dt))),
                 c(matrix(16,1,30*(1/dt)),matrix(0,1,tsteps - 30*(1/dt))),
                 c(matrix(16,1,40*(1/dt)),matrix(0,1,tsteps - 40*(1/dt))))

pref <- "Stepped_omega_"
plot_diff_acd_output(acqdistn,"~/Documents/Hetero_res_and_f/plots",02, omega_M, submic_M, wildtype,pref)

####*** # (4) run with high short or long low omega ***######################################################################################################################################################################################################################################################################################################
omega_M <- rbind(c(matrix(30,1,10*(1/dt)),matrix(0,1,tsteps - 10*(1/dt))),
                 c(matrix(30,1,20*(1/dt)),matrix(0,1,tsteps - 20*(1/dt))),
                 c(matrix(0.4,1,50*(1/dt)),matrix(0,1,tsteps - 50*(1/dt))),
                 c(matrix(0.4,1,100*(1/dt)),matrix(0,1,tsteps - 100*(1/dt))))

pref <- "HSoLL_"
plot_diff_acd_output(acqdistn,"~/Documents/Hetero_res_and_f/plots",02, omega_M, submic_M, wildtype,pref)

####*** # (5) run with PULSING high omega ***######################################################################################################################################################################################################################################################################################################
exp_dec <- -0.01
gap_doses <- 120
dose_h <- 16 
fn_pulse <- dose_h*exp(exp_dec*(0:gap_doses))
omega_M <- rbind(c(dose_h*exp(exp_dec*(0:(tsteps-1))) ),
                 c(fn_pulse,fn_pulse,                                     dose_h*exp(exp_dec*(0:(tsteps-2*(gap_doses+1)-1)))  ),
                 c(fn_pulse,fn_pulse,fn_pulse,fn_pulse,fn_pulse,          dose_h*exp(exp_dec*(0:(tsteps-5*(gap_doses+1)-1)))  ),
                 c(fn_pulse,fn_pulse,fn_pulse,fn_pulse,fn_pulse,
                   fn_pulse,fn_pulse,fn_pulse,fn_pulse,fn_pulse,fn_pulse, dose_h*exp(exp_dec*(0:(tsteps-11*(gap_doses+1)-1))) ) )

pref <- "Pulse_high_"
plot_diff_acd_output(acqdistn,"~/Documents/Hetero_res_and_f/plots",02, omega_M, submic_M, wildtype,pref)


####*** # (6) run with PULSING low omega ***######################################################################################################################################################################################################################################################################################################
exp_dec <- -0.01
gap_doses <- 120
dose_h <- 2
fn_pulse <- dose_h*exp(exp_dec*(0:gap_doses))
omega_M <- rbind(c(dose_h*exp(exp_dec*(0:(tsteps-1))) ),
                 c(fn_pulse,fn_pulse,                                     dose_h*exp(exp_dec*(0:(tsteps-2*(gap_doses+1)-1)))  ),
                 c(fn_pulse,fn_pulse,fn_pulse,fn_pulse,fn_pulse,          dose_h*exp(exp_dec*(0:(tsteps-5*(gap_doses+1)-1)))  ),
                 c(fn_pulse,fn_pulse,fn_pulse,fn_pulse,fn_pulse,
                   fn_pulse,fn_pulse,fn_pulse,fn_pulse,fn_pulse,fn_pulse, dose_h*exp(exp_dec*(0:(tsteps-11*(gap_doses+1)-1))) ) )

pref <- "Pulse_low_"
plot_diff_acd_output(acqdistn,"~/Documents/Hetero_res_and_f/plots",02, omega_M, submic_M, wildtype,pref)

####*** # (7) run with high-low or low-high ***######################################################################################################################################################################################################################################################################################################
omega_M <- rbind(c(matrix(16,1,20*(1/dt)), matrix(0.4,1,20*(1/dt)),matrix(0,1,tsteps - 40*(1/dt))),
                 c(matrix(0.4,1,20*(1/dt)),matrix(16,1,20*(1/dt)), matrix(0,1,tsteps - 40*(1/dt))),
                 c(matrix(16,1,50*(1/dt)), matrix(0.4,1,50*(1/dt)),matrix(0,1,tsteps - 100*(1/dt))),
                 c(matrix(0.4,1,50*(1/dt)),matrix(16,1,50*(1/dt)), matrix(0,1,tsteps - 100*(1/dt))))

pref <- "Highlo_lohigh_"
plot_diff_acd_output(acqdistn,"~/Documents/Hetero_res_and_f/plots",02, omega_M, submic_M, wildtype,pref)


#####################################################################################################################################################################################################################################################################################################
#####################################################################################################################################################################################################################################################################################################
####**** DIFF MIC SHAPE ***###############################################################################################################################################################################################################################################################################################################
#####################################################################################################################################################################################################################################################################################################
submic_M <- c(2,2,2,2) # flat then drop at 0.9*MIC

####*** # (0) run with zero omega ***######################################################################################################################################################################################################################################################################################################
Sv20<-ec_srs_funcf_mean_varsr(tsteps,home, c(kk,0),iniv,M0,acqdistn,dt,2)
# What happens? 
mm20<-c(); ll<-dim(Sv20$M)[3];
ss<-seq(0,ll,1/dt) # Don't want to grab all 
for(i in 2:length(ss)){
  mm220<-as.data.frame(melt(Sv20$M[,,ss[i]])); mm220$tstep=ss[i]*dt; mm20<-rbind(mm20,mm220);
} 
colnames(mm20)<-c("x","y","z","t"); 
setwd(plots)
# Grab a subset
sub<-c(1/dt,100,250,500,750,1000,1500,2000,2500,seq(3000,4001,500))*dt
w<-which(mm20[,"t"] %in% sub)
# plots
p1<-ggplot(mm20[w,],aes(x,y,fill=z))  + facet_wrap( ~ t, ncol=3) + ggtitle(paste("w = ", 0,sep="")) + 
  scale_fill_gradient("Proportion", limits=c(0,1),low="white", high="red") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_tile() + 
  scale_y_continuous(breaks=c(1,nfit),"Relative fitness levels",labels=c("Least","Most")) + 
  scale_x_continuous(breaks=c(mres,1),"Resistance levels",labels=c("Most","Least"))
p1
ggsave(paste("2_norsd_Array_w=",0,"_06.pdf",sep=""))

# plots only proportion with "resistance" - not background initial resistance
M20n <- Sv20$M;
M20n[susr,susf,] <- 0; 
for (i in 1:dim(Sv20$M)[3]){ M20n[,,i] <- M20n[,,i]/sum(M20n[,,i])} # re-normalise

# What happens? 
mm20<-c() 
ll<-dim(Sv20$M)[3];
ss<-seq(0,ll,1/dt) # Don't want to grab all 
for(i in 2:length(ss)){
  mm220<-as.data.frame(melt(M20n[,,ss[i]])); 
  mm220$tstep=ss[i]*dt; 
  mm20<-rbind(mm20,mm220);
} 
colnames(mm20)<-c("x","y","z","t"); 

# plots
p1<-ggplot(mm20[w,],aes(x,y,fill=z))  + facet_wrap( ~ t, ncol=3) + ggtitle(paste("w = ", 0,sep=""))
p1<-p1 + scale_fill_gradient("Proportion", limits=c(0,1),low="white", high="red",guide = FALSE)
p1<-p1 + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p1<-p1 + geom_tile() + scale_y_continuous(breaks=c(1,nfit),"Relative fitness levels",labels=c("Least","Most")) + scale_x_continuous(breaks=c(mres,1),"Resistance levels",labels=c("Most","Least"))
p1
ggsave(paste("2_norsdR_Array_w=",0,"_06.pdf",sep=""))

pdf("2_B_no_omega.pdf")
plot(Sv20$B,type="l", ylim = c(0,100), ylab = c("Percentage"), xlab = "Time")
dev.off()


####*** # (1) run with constant omega ***######################################################################################################################################################################################################################################################################################################
omega_M <- c(24,16,2,0.4)
pref <- "Constant_omega_"
plot_diff_acd_output(acqdistn,"~/Documents/Hetero_res_and_f/plots",02, omega_M, submic_M, wildtype,pref)

####*** # (2) run with timevarying omega ***######################################################################################################################################################################################################################################################################################################
omega_M <- rbind(c(matrix(16,1,20),matrix(0,1,tsteps - 20)),c(16,16*exp(-0.1*(2:tsteps))),
                 c(matrix(16,1,20),matrix(0.4,1,20),matrix(16,1,20),matrix(0,1,tsteps - 60)), c(matrix(0,1,20),matrix(16,1,tsteps - 20)))
pref <- "TV_omega_"
plot_diff_acd_output(acqdistn,"~/Documents/Hetero_res_and_f/plots",02, omega_M, submic_M, wildtype,pref)

####*** # (3) run with stepped omega ***######################################################################################################################################################################################################################################################################################################
omega_M <- rbind(c(matrix(16,1,10*(1/dt)),matrix(0,1,tsteps - 10*(1/dt))),
                 c(matrix(16,1,20*(1/dt)),matrix(0,1,tsteps - 20*(1/dt))),
                 c(matrix(16,1,30*(1/dt)),matrix(0,1,tsteps - 30*(1/dt))),
                 c(matrix(16,1,40*(1/dt)),matrix(0,1,tsteps - 40*(1/dt))))
pref <- "Stepped_omega_"
plot_diff_acd_output(acqdistn,"~/Documents/Hetero_res_and_f/plots",02, omega_M, submic_M, wildtype,pref)

####*** # (4) run with high short or long low omega ***######################################################################################################################################################################################################################################################################################################
omega_M <- rbind(c(matrix(30,1,10*(1/dt)),matrix(0,1,tsteps - 10*(1/dt))),
                 c(matrix(30,1,20*(1/dt)),matrix(0,1,tsteps - 20*(1/dt))),
                 c(matrix(0.4,1,50*(1/dt)),matrix(0,1,tsteps - 50*(1/dt))),
                 c(matrix(0.4,1,100*(1/dt)),matrix(0,1,tsteps - 100*(1/dt))))
pref <- "HSoLL_"
plot_diff_acd_output(acqdistn,"~/Documents/Hetero_res_and_f/plots",02, omega_M, submic_M, wildtype,pref)

####*** # (5) run with PULSING high omega ***######################################################################################################################################################################################################################################################################################################
exp_dec <- -0.01
gap_doses <- 120
dose_h <- 16 
fn_pulse <- dose_h*exp(exp_dec*(0:gap_doses))
omega_M <- rbind(c(dose_h*exp(exp_dec*(0:(tsteps-1))) ),
                 c(fn_pulse,fn_pulse,                                     dose_h*exp(exp_dec*(0:(tsteps-2*(gap_doses+1)-1)))  ),
                 c(fn_pulse,fn_pulse,fn_pulse,fn_pulse,fn_pulse,          dose_h*exp(exp_dec*(0:(tsteps-5*(gap_doses+1)-1)))  ),
                 c(fn_pulse,fn_pulse,fn_pulse,fn_pulse,fn_pulse,
                   fn_pulse,fn_pulse,fn_pulse,fn_pulse,fn_pulse,fn_pulse, dose_h*exp(exp_dec*(0:(tsteps-11*(gap_doses+1)-1))) ) )

pref <- "Pulse_high_"
plot_diff_acd_output(acqdistn,"~/Documents/Hetero_res_and_f/plots",02, omega_M, submic_M, wildtype,pref)


####*** # (6) run with PULSING low omega ***######################################################################################################################################################################################################################################################################################################
exp_dec <- -0.01
gap_doses <- 120
dose_h <- 2
fn_pulse <- dose_h*exp(exp_dec*(0:gap_doses))
omega_M <- rbind(c(dose_h*exp(exp_dec*(0:(tsteps-1))) ),
                 c(fn_pulse,fn_pulse,                                     dose_h*exp(exp_dec*(0:(tsteps-2*(gap_doses+1)-1)))  ),
                 c(fn_pulse,fn_pulse,fn_pulse,fn_pulse,fn_pulse,          dose_h*exp(exp_dec*(0:(tsteps-5*(gap_doses+1)-1)))  ),
                 c(fn_pulse,fn_pulse,fn_pulse,fn_pulse,fn_pulse,
                   fn_pulse,fn_pulse,fn_pulse,fn_pulse,fn_pulse,fn_pulse, dose_h*exp(exp_dec*(0:(tsteps-11*(gap_doses+1)-1))) ) )

pref <- "Pulse_low_"
plot_diff_acd_output(acqdistn,"~/Documents/Hetero_res_and_f/plots",02, omega_M, submic_M, wildtype,pref)

####*** # (7) run with high-low or low-high ***######################################################################################################################################################################################################################################################################################################
omega_M <- rbind(c(matrix(16,1,20*(1/dt)), matrix(0.4,1,20*(1/dt)),matrix(0,1,tsteps - 40*(1/dt))),
                 c(matrix(0.4,1,20*(1/dt)),matrix(16,1,20*(1/dt)), matrix(0,1,tsteps - 40*(1/dt))),
                 c(matrix(16,1,50*(1/dt)), matrix(0.4,1,50*(1/dt)),matrix(0,1,tsteps - 100*(1/dt))),
                 c(matrix(0.4,1,50*(1/dt)),matrix(16,1,50*(1/dt)), matrix(0,1,tsteps - 100*(1/dt))))

pref <- "Highlo_lohigh_"
plot_diff_acd_output(acqdistn,"~/Documents/Hetero_res_and_f/plots",02, omega_M, submic_M, wildtype,pref)

