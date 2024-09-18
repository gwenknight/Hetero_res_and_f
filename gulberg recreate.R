### Gulberg plot recreate

##*** Libraries needed
library(mvtnorm);library(plyr); library(ggplot2);library(reshape2); library(tidyverse)
library(deSolve);library(grid);library(gtools); library(directlabels);
theme_set(theme_gray(base_size = 24)); 
##*** Locations
home<-"~/Documents/Hetero_res_and_f/"
plots<-paste(home,"plots",sep="")
setwd(home)

theme_set(theme_bw(base_size = 34))
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

##*** Code needed
# Loads functions for generalised mean function and simulation model that uses the generalised function
# Also 2 fitness level ode function, original Sourya model as a difference model and multiplots function
source("ec_generalised_function_withr.R") 

##*** Setting up
# Number of discrete fitness levels? Resistance levels?
nfit = 5; 
mres = 5;
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
setwd(plots)
ggsave("acqdistn_06.pdf",width=14,height=10)

# ## Try another
# aa <- matrix(0,25,25)
# aa[1:5,21:25] <- acqdistn
# acqdistn <- aa
# nfit = 25; 
# mres = 25;
# # Array of distribution of fitness and resistance c[resistance, fitness]
# M0 <- array(0,c(mres,nfit,10))

# most low level acquisition distribution 
x <- seq(1/mres,1,1/mres)
y <- seq(1/nfit,1,1/nfit)
f <- function(x, y) dmvnorm(cbind(x, y), mean = c(0.0000001, 0.5),sigma = diag(2)/10)
z <- outer(x, y, FUN = f); z <- z/sum(z) # Generate and normalise
dev.off(); persp(x, y, z, theta = 60, phi = 30, ticktype = "detailed")
acqdistn<-z

# Initial conditions
iniv<-c(99,1,0)
#iniv<-c(60,39,1)

#############********************************************** LOAD UP TO HERE *********************************************########
dt=0.1
tsteps<-500*(1/dt)
omega1 <- 1
omega2 <- 0.4
omega3 <- 0.2
omega4 <- 0.1
Sv20<-ec_funcf_mean_varsr(tsteps,home, c(omega1),iniv,M0,acqdistn,dt,500)
Sv15<-ec_funcf_mean_varsr(tsteps,home, c(omega2),iniv,M0,acqdistn,dt,500)
Sv10<-ec_funcf_mean_varsr(tsteps,home, c(omega3),iniv,M0,acqdistn,dt,500)
Sv05<-ec_funcf_mean_varsr(tsteps,home, c(omega4),iniv,M0,acqdistn,dt,500)
## NEED TO SPEED IT UP?? Fast for 5 x 5... ~6 sec on laptop 

# Plot proportions in each fitness / resistance level over time
pp<-c();
ll<-dim(Sv20$M)[3];
ss<-seq(0,ll,1/dt) # Don't want to grab all 
for(i in 2:length(ss)){
  pp220<-c(ss[i]*dt,colSums(Sv20$M[,,ss[i]]), rowSums(Sv20$M[,,ss[i]]),omega1)
  pp215<-c(ss[i]*dt,colSums(Sv15$M[,,ss[i]]), rowSums(Sv15$M[,,ss[i]]),omega2)
  pp210<-c(ss[i]*dt,colSums(Sv10$M[,,ss[i]]), rowSums(Sv10$M[,,ss[i]]),omega3)
  pp205<-c(ss[i]*dt,colSums(Sv05$M[,,ss[i]]), rowSums(Sv05$M[,,ss[i]]),omega4)
  pp<-rbind(pp,pp220,pp215,pp210,pp205);
} 
pp<-as.data.frame(pp);
colnames(pp) <- c("t","F1","F2","F3","F4","F5","R1","R2","R3","R4","R5","w")
#colnames(pp)<-c("t","Fitness level 1\n(low)","Fitness level 2","Fitness level 3","Fitness level 4","Fitness level 5\n(high)","Res. level 1\n(low)","Res. level 2","Res. level 3","Res. level 4","Res. level 5\n(high)","w"); 
pp2<-pp %>% pivot_longer(cols = F1:R5)
pp2$type  <- substring(pp2$name,1,1)
pp2$level  <- substring(pp2$name,2,2)
theme_set(theme_bw(base_size = 34)); 

pp2 <- pp2 %>% group_by(type,w,t) %>% arrange(desc(level),.by_group = TRUE) %>% mutate(cs = cumsum(value))

#g<-ggplot(pp2,aes(x=t,y=value,colour=factor(w))) + facet_wrap(~variable,ncol=5) + geom_line(size=2) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
#g<-g + scale_x_continuous("Generation") + scale_y_continuous("Proportion") + scale_colour_manual(values=cbPalette,"Abx\nlevel",labels=c(0.05,0.1,0.2))
#g # Suggests that although v similar proportions in the most fit fewer are in the higher resistance levels with low level antibiotics use. In fact with this model
# the same rate of selection for no cost mutations is seen whether there is high or low anitbiotic use 


## Plot like in Gulberg
#pp2n<-pp2[7501:15000,]
#w<-intersect(which(pp2n$w==5),c(which(pp2n$t==100),which(pp2n$t==200),which(pp2n$t==300),which(pp2n$t==400),which(pp2n$t==500)))
#ggplot(pp2n[w,], aes(x=t, y= value,colour=factor(variable))) + geom_point(aes(shape = factor(variable)),size=5)

ggplot(pp2 %>% filter(t %in% c(100,200,300,400,500), grepl('Res.', variable)), aes(x=t, y= value, colour=factor(variable))) + 
  scale_color_brewer(palette = "RdYlBu") + 
  geom_point(size = 5) + #aes(shape = factor(variable)),size=5) + 
  facet_wrap(~w) + 
  scale_y_continuous(trans='log10', limits = c(0.0001,1), breaks = c(1,0.1,0.01,0.001)) + 
  geom_line()


ggplot(pp2 %>% filter(type == "R"), aes(x=t, y=value, colour=factor(level))) + 
  scale_color_brewer(palette = "RdYlBu") + 
  geom_point(size = 3) + #aes(shape = factor(variable)),size=5) + 
  facet_wrap(~w) + 
  scale_y_continuous(limits = c(0.0001,1), breaks = c(1,0.1,0.01,0.001)) + 
  geom_line()

