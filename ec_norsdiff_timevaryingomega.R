# Time varying omega play

##### No difference betwen susceptible and resistant populations

### Run as for simple e coli model

##*** Libraries needed
library(mvtnorm);library(plyr); library(ggplot2);library(reshape2);library(deSolve);library(grid);library(gtools); library(directlabels); library(cowplot)
theme_set(theme_gray(base_size = 24)); 
##*** Locations
home<-"~/Documents/Hetero_res_and_f/"
plots<-paste(home,"plots/timevarying",sep="")
setwd(home)

theme_set(theme_bw(base_size = 34))
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

##*** Code needed
# Loads functions for generalised mean function and simulation model that uses the generalised function
# Also 2 fitness level ode function, original Sourya model as a difference model and multiplots function
source("ec_nosrdiff_generalised_function_withr.R") 

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
ggsave("norsd_acqdistn_06.pdf",width=14,height=10)

# ## Try another
# aa <- matrix(0,25,25)
# aa[1:5,21:25] <- acqdistn
# acqdistn <- aa
# nfit = 25; 
# mres = 25;
# # Array of distribution of fitness and resistance c[resistance, fitness]
# M0 <- array(0,c(mres,nfit,10))


# Initial conditions
iniv<-c(20,80) # resource vs. bug #### TREATMENT => AT EQUILIBRIUM?
#iniv<-c(60,39,1)
# M0 now needs initial condition - mostly "susceptible"
# vs_mic <- seq(0.1,32,length.out = m) MIC of susceptible depends on available resistance levels
susr <- 1; susf <- 5
M0[susr,susf,1] <- 1

#############********************************************** LOAD UP TO HERE *********************************************########
dt=0.1
tsteps<-500*(1/dt)
omega1 <- c(matrix(16,1,20),matrix(0,1,tsteps - 20))
omega2 <- c(16,16*exp(-1*(2:20)),matrix(0,1,tsteps-20))
omega3 <- c(matrix(16,1,20),matrix(0.4,1,20),matrix(16,1,20),matrix(0,1,tsteps - 60))
omega4 <- c(matrix(0,1,20),matrix(16,1,tsteps - 20))
omegam <- as.data.frame(cbind(seq(1,tsteps),omega1,omega2,omega3, omega4)); colnames(omegam) <- c("time","1","2","3","4"); momegam <- melt(omegam,id.vars = "time")
go<-ggplot(momegam,aes(x=time,y=value,colour=variable)) + geom_line(size=2) + scale_x_continuous(lim=c(0,200), breaks = seq(0,200,50), labels = dt*seq(0,200,50)) + scale_y_continuous("Omega") + scale_color_discrete("")
ggsave(paste("norsd_omegaconcs_06.pdf",sep=""))

kk <- 500
Sv20<-ec_srs_funcf_mean_varsr(tsteps,home, c(kk,omega1),iniv,M0,acqdistn,dt)
Sv15<-ec_srs_funcf_mean_varsr(tsteps,home, c(kk,omega2),iniv,M0,acqdistn,dt)
Sv10<-ec_srs_funcf_mean_varsr(tsteps,home, c(kk,omega3),iniv,M0,acqdistn,dt)
Sv05<-ec_srs_funcf_mean_varsr(tsteps,home, c(kk,omega4),iniv,M0,acqdistn,dt)
## NEED TO SPEED IT UP?? Fast for 5 x 5... ~6 sec on laptop 

# What happens? 
mm20<-c() ; mm10<-c() ; mm05<-c() ; mm15<-c() 
ll<-dim(Sv20$M)[3];
ss<-seq(0,ll,1/dt) # Don't want to grab all 
for(i in 2:length(ss)){
  mm220<-as.data.frame(melt(Sv20$M[,,ss[i]])); 
  mm215<-as.data.frame(melt(Sv15$M[,,ss[i]])); 
  mm210<-as.data.frame(melt(Sv10$M[,,ss[i]])); 
  mm205<-as.data.frame(melt(Sv05$M[,,ss[i]])); 
  mm220$tstep=ss[i]*dt; mm215$tstep=ss[i]*dt; mm210$tstep=ss[i]*dt; mm205$tstep=ss[i]*dt # To go to generations
  mm20<-rbind(mm20,mm220);mm15<-rbind(mm15,mm215); mm10<-rbind(mm10,mm210) ; mm05<-rbind(mm05,mm205) 
} 
colnames(mm20)<-c("x","y","z","t"); colnames(mm15)<-c("x","y","z","t"); colnames(mm10)<-c("x","y","z","t") ;colnames(mm05)<-c("x","y","z","t")
#mm20$x<-seq(mres,1,-1); mm10$x<-seq(mres,1,-1); mm05$x<-seq(mres,1,-1)
setwd(plots)
# Grab a subset
#sub<-c(1/dt,250,500,750,1000,1500,2000,2500,seq(3000,4001,1000),4500,tsteps)*dt
sub<-c(1/dt,100,250,500,750,1000,1500,2000,2500,seq(3000,4001,500))*dt
w<-which(mm20[,"t"] %in% sub)
# plots
p1<-ggplot(mm20[w,],aes(x,y,fill=z))  + facet_wrap( ~ t, ncol=3) + ggtitle(paste("w = ", 1,sep=""))
p1<-p1 + scale_fill_gradient("Proportion", limits=c(0,1),low="white", high="red",guide = FALSE)
p1<-p1 + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p1<-p1 + geom_tile() + scale_y_continuous(breaks=c(1,nfit),"Relative fitness levels",labels=c("Least","Most")) + scale_x_continuous(breaks=c(mres,1),"Resistance levels",labels=c("Most","Least"))
p1
ggsave(paste("norsd_Array_w=",omega1,"_06.pdf",sep=""))
p2<-ggplot(mm15[w,],aes(x,y,fill=z))  + facet_wrap( ~ t, ncol=3) + ggtitle(paste("w = ", 2,sep=""))
p2<-p2 + scale_fill_gradient("Proportion", limits=c(0,1),low="white", high="red",guide=FALSE)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p2<-p2 + geom_tile() + scale_y_continuous(breaks=c(1,nfit),"Relative fitness levels",labels=c("Least","Most")) + scale_x_continuous(breaks=c(mres,1),"Resistance levels",labels=c("Most","Least"))
p2
ggsave(paste("norsd_Array_w=",omega2,"_06.pdf",sep=""))
p3<-ggplot(mm10[w,],aes(x,y,fill=z))  + facet_wrap( ~ t, ncol=3) + ggtitle(paste("w = ", 3,sep=""))
p3<-p3 + scale_fill_gradient("Proportion", limits=c(0,1),low="white", high="red",guide = FALSE)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p3<-p3 + geom_tile() + scale_y_continuous(breaks=c(1,nfit),"Relative fitness levels",labels=c("Least","Most")) + scale_x_continuous(breaks=c(mres,1),"Resistance levels",labels=c("Most","Least"))
p3
ggsave(paste("norsd_Array_w=",omega3,"_06.pdf",sep=""))
p4<-ggplot(mm05[w,],aes(x,y,fill=z))  + facet_wrap( ~ t, ncol=3) + ggtitle(paste("w = ", 4,sep=""))
p4<-p4 + scale_fill_gradient("Proportion", limits=c(0,1),low="white", high="red",guide=FALSE)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p4<-p4 + geom_tile() + scale_y_continuous(breaks=c(1,nfit),"Relative fitness levels",labels=c("Least","Most")) + scale_x_continuous(breaks=c(mres,1),"Resistance levels",labels=c("Most","Least"))
p4
ggsave(paste("norsd_Array_w=",omega4,"_06.pdf",sep=""))


setwd(plots)
pdf("norsd_Array_w_all.pdf",width=18,height=18)
multiplot(p1,p2,p4,p3,cols=2)
dev.off()

## plot U, S & R over time
Mu <- as.data.frame(cbind(seq(0, tsteps, 1),matrix(0,tsteps+1,1), Sv05$U, Sv10$U,Sv15$U, Sv20$U))
Mb <- as.data.frame(cbind(seq(0, tsteps, 1),matrix(1,tsteps+1,1), Sv05$B, Sv10$B,Sv15$B, Sv20$B))
Mmf <- as.data.frame(cbind(seq(0, tsteps, 1),matrix(3,tsteps+1,1), Sv05$meanf[,1], Sv10$meanf[,1],Sv15$meanf[,1], Sv20$meanf[,1]))
Mmr <- as.data.frame(cbind(seq(0, tsteps, 1),matrix(4,tsteps+1,1), Sv05$meanf[,2], Sv10$meanf[,2],Sv15$meanf[,2], Sv20$meanf[,2]))
Mhigh <- as.data.frame(cbind(seq(0, tsteps, 1),matrix(5,tsteps+1,1), Sv05$M[5,5,], Sv10$M[5,5,], Sv15$M[5,5,], Sv20$M[5,5,]))
Musr <- rbind(Mu, Mb,Mmf,Mmr,Mhigh)
colnames(Musr) <- c("t", "pop",4,3,2,1)
Msrm <- melt(Musr, id.vars = c("t","pop"))
facet_names <- c(`0` = "U", `1` = "B", `3` = "mean fit", `4` = "mean res", `5` = "Highest fit/res")
ggplot(Msrm, aes(x=t, y = value, colour = variable)) + geom_line() + facet_wrap(~pop,labeller = as_labeller(facet_names), scales = "free")
ggsave("norsd_TimeSeries_output_06.pdf")
# number in highest fitness changes but mean r and f don't? 

# How is the bacterial load affected? want to control.
colnames(Mb) <- c("t", "pop",4,3,2,1)
Mbm <- melt(Mb, id.vars = c("t","pop"))
pmb<-ggplot(Mbm, aes(x=t, y = value, colour = variable)) + geom_line() + scale_x_continuous(lim = c(0,1200)) + scale_y_continuous(lim = c(0,100))
ggsave(paste("norsd_TimeSeries_Boutput_06.pdf",sep=""))

# plots only proportion with "resistance" - not background initial resistance
M20n <- Sv20$M;M10n <- Sv10$M;M15n <- Sv15$M;M05n <- Sv05$M
M20n[susr,susf,] <- 0; M10n[susr,susf,] <- 0; M05n[susr,susf,] <- 0; M15n[susr,susf,] <- 0; 
for (i in 1:dim(Sv20$M)[3]){ M20n[,,i] <- M20n[,,i]/sum(M20n[,,i]) # re-normalise
M15n[,,i] <- M15n[,,i]/sum(M15n[,,i]) # re-normalise
M10n[,,i] <- M10n[,,i]/sum(M10n[,,i]) # re-normalise
M05n[,,i] <- M05n[,,i]/sum(M05n[,,i])} # re-normalise

# What happens? 
mm20<-c() ; mm10<-c() ; mm05<-c() ; mm15<-c() 
ll<-dim(Sv20$M)[3];
ss<-seq(0,ll,1/dt) # Don't want to grab all 
for(i in 2:length(ss)){
  mm220<-as.data.frame(melt(M20n[,,ss[i]])); 
  mm215<-as.data.frame(melt(M15n[,,ss[i]])); 
  mm210<-as.data.frame(melt(M10n[,,ss[i]])); 
  mm205<-as.data.frame(melt(M05n[,,ss[i]])); 
  mm220$tstep=ss[i]*dt; mm215$tstep=ss[i]*dt; mm210$tstep=ss[i]*dt; mm205$tstep=ss[i]*dt # To go to generations
  mm20<-rbind(mm20,mm220);mm15<-rbind(mm15,mm215); mm10<-rbind(mm10,mm210) ; mm05<-rbind(mm05,mm205) 
} 
colnames(mm20)<-c("x","y","z","t"); colnames(mm15)<-c("x","y","z","t"); colnames(mm10)<-c("x","y","z","t") ;colnames(mm05)<-c("x","y","z","t")

# plots
p1<-ggplot(mm20[w,],aes(x,y,fill=z))  + facet_wrap( ~ t, ncol=3) + ggtitle(paste("w = ", 1,sep=""))
p1<-p1 + scale_fill_gradient("Proportion", limits=c(0,1),low="white", high="red",guide = FALSE)
p1<-p1 + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p1<-p1 + geom_tile() + scale_y_continuous(breaks=c(1,nfit),"Relative fitness levels",labels=c("Least","Most")) + scale_x_continuous(breaks=c(mres,1),"Resistance levels",labels=c("Most","Least"))
p1
ggsave(paste("norsdR_Array_w=",omega1,"_06.pdf",sep=""))
p2<-ggplot(mm15[w,],aes(x,y,fill=z))  + facet_wrap( ~ t, ncol=3) + ggtitle(paste("w = ", 2,sep=""))
p2<-p2 + scale_fill_gradient("Proportion", limits=c(0,1),low="white", high="red",guide=FALSE)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p2<-p2 + geom_tile() + scale_y_continuous(breaks=c(1,nfit),"Relative fitness levels",labels=c("Least","Most")) + scale_x_continuous(breaks=c(mres,1),"Resistance levels",labels=c("Most","Least"))
p2
ggsave(paste("norsdR_Array_w=",omega2,"_06.pdf",sep=""))
p3<-ggplot(mm10[w,],aes(x,y,fill=z))  + facet_wrap( ~ t, ncol=3) + ggtitle(paste("w = ", 3,sep=""))
p3<-p3 + scale_fill_gradient("Proportion", limits=c(0,1),low="white", high="red",guide = FALSE)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p3<-p3 + geom_tile() + scale_y_continuous(breaks=c(1,nfit),"Relative fitness levels",labels=c("Least","Most")) + scale_x_continuous(breaks=c(mres,1),"Resistance levels",labels=c("Most","Least"))
p3
ggsave(paste("norsdR_Array_w=",omega3,"_06.pdf",sep=""))
p4<-ggplot(mm05[w,],aes(x,y,fill=z))  + facet_wrap( ~ t, ncol=3) + ggtitle(paste("w = ", 4,sep=""))
p4<-p4 + scale_fill_gradient("Proportion", limits=c(0,1),low="white", high="red",guide=FALSE)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p4<-p4 + geom_tile() + scale_y_continuous(breaks=c(1,nfit),"Relative fitness levels",labels=c("Least","Most")) + scale_x_continuous(breaks=c(mres,1),"Resistance levels",labels=c("Most","Least"))
p4
ggsave(paste("norsdR_Array_w=",omega4,"_06.pdf",sep=""))

setwd(plots)
pdf("norsdR_Array_w_all.pdf",width=18,height=18)
multiplot(p1,p2,p4,p3,cols=2)
dev.off()

pdf("all_norsd_all_Array_w_all.pdf",width=18,height=18)
multiplot(go,p1,p2,pmb,p4,p3,cols=2)
dev.off()

#############********************************************** STEP DOWN FUNCTIONS *********************************************########
## prefix = StD
pref <- "Fn_STD_"
dt=0.1
tsteps<-500*(1/dt)
omega1 <- c(matrix(16,1,10*(1/dt)),matrix(0,1,tsteps - 10*(1/dt)))
omega2 <- c(matrix(16,1,20*(1/dt)),matrix(0,1,tsteps - 20*(1/dt)))
omega3 <- c(matrix(16,1,30*(1/dt)),matrix(0,1,tsteps - 30*(1/dt)))
omega4 <- c(matrix(16,1,40*(1/dt)),matrix(0,1,tsteps - 40*(1/dt)))
omegam <- as.data.frame(cbind(seq(1,tsteps),omega1,omega2,omega3, omega4)); colnames(omegam) <- c("time","1","2","3","4"); momegam <- melt(omegam,id.vars = "time")
go<-ggplot(momegam,aes(x=time,y=value,colour=variable)) + geom_line(size=2) + scale_x_continuous(lim=c(0,200), breaks = seq(0,200,50), labels = dt*seq(0,200,50)) + scale_y_continuous("Omega") + scale_color_discrete("")
ggsave(paste(pref,"norsd_omegaconcs_06.pdf",sep=""))

kk <- 500
Sv20<-ec_srs_funcf_mean_varsr(tsteps,home, c(kk,omega1),iniv,M0,acqdistn,dt)
Sv15<-ec_srs_funcf_mean_varsr(tsteps,home, c(kk,omega2),iniv,M0,acqdistn,dt)
Sv10<-ec_srs_funcf_mean_varsr(tsteps,home, c(kk,omega3),iniv,M0,acqdistn,dt)
Sv05<-ec_srs_funcf_mean_varsr(tsteps,home, c(kk,omega4),iniv,M0,acqdistn,dt)
## NEED TO SPEED IT UP?? Fast for 5 x 5... ~6 sec on laptop 

# What happens? 
mm20<-c() ; mm10<-c() ; mm05<-c() ; mm15<-c() 
ll<-dim(Sv20$M)[3];
ss<-seq(0,ll,1/dt) # Don't want to grab all 
for(i in 2:length(ss)){
  mm220<-as.data.frame(melt(Sv20$M[,,ss[i]])); 
  mm215<-as.data.frame(melt(Sv15$M[,,ss[i]])); 
  mm210<-as.data.frame(melt(Sv10$M[,,ss[i]])); 
  mm205<-as.data.frame(melt(Sv05$M[,,ss[i]])); 
  mm220$tstep=ss[i]*dt; mm215$tstep=ss[i]*dt; mm210$tstep=ss[i]*dt; mm205$tstep=ss[i]*dt # To go to generations
  mm20<-rbind(mm20,mm220);mm15<-rbind(mm15,mm215); mm10<-rbind(mm10,mm210) ; mm05<-rbind(mm05,mm205) 
} 
colnames(mm20)<-c("x","y","z","t"); colnames(mm15)<-c("x","y","z","t"); colnames(mm10)<-c("x","y","z","t") ;colnames(mm05)<-c("x","y","z","t")
#mm20$x<-seq(mres,1,-1); mm10$x<-seq(mres,1,-1); mm05$x<-seq(mres,1,-1)
setwd(plots)
# Grab a subset
#sub<-c(1/dt,250,500,750,1000,1500,2000,2500,seq(3000,4001,1000),4500,tsteps)*dt
sub<-c(1/dt,100,250,500,750,1000,1500,2000,2500,seq(3000,4001,500))*dt
w<-which(mm20[,"t"] %in% sub)
# plots
p1<-ggplot(mm20[w,],aes(x,y,fill=z))  + facet_wrap( ~ t, ncol=3) + ggtitle(paste("w = ", 1,sep=""))
p1<-p1 + scale_fill_gradient("Proportion", limits=c(0,1),low="white", high="red",guide = FALSE)
p1<-p1 + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p1<-p1 + geom_tile() + scale_y_continuous(breaks=c(1,nfit),"Relative fitness levels",labels=c("Least","Most")) + scale_x_continuous(breaks=c(mres,1),"Resistance levels",labels=c("Most","Least"))
p1
ggsave(paste(pref,"norsd_Array_w=",omega1,"_06.pdf",sep=""))
p2<-ggplot(mm15[w,],aes(x,y,fill=z))  + facet_wrap( ~ t, ncol=3) + ggtitle(paste("w = ", 2,sep=""))
p2<-p2 + scale_fill_gradient("Proportion", limits=c(0,1),low="white", high="red",guide=FALSE)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p2<-p2 + geom_tile() + scale_y_continuous(breaks=c(1,nfit),"Relative fitness levels",labels=c("Least","Most")) + scale_x_continuous(breaks=c(mres,1),"Resistance levels",labels=c("Most","Least"))
p2
ggsave(paste(pref,"norsd_Array_w=",omega2,"_06.pdf",sep=""))
p3<-ggplot(mm10[w,],aes(x,y,fill=z))  + facet_wrap( ~ t, ncol=3) + ggtitle(paste("w = ", 3,sep=""))
p3<-p3 + scale_fill_gradient("Proportion", limits=c(0,1),low="white", high="red",guide = FALSE)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p3<-p3 + geom_tile() + scale_y_continuous(breaks=c(1,nfit),"Relative fitness levels",labels=c("Least","Most")) + scale_x_continuous(breaks=c(mres,1),"Resistance levels",labels=c("Most","Least"))
p3
ggsave(paste(pref,"norsd_Array_w=",omega3,"_06.pdf",sep=""))
p4<-ggplot(mm05[w,],aes(x,y,fill=z))  + facet_wrap( ~ t, ncol=3) + ggtitle(paste("w = ", 4,sep=""))
p4<-p4 + scale_fill_gradient("Proportion", limits=c(0,1),low="white", high="red",guide=FALSE)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p4<-p4 + geom_tile() + scale_y_continuous(breaks=c(1,nfit),"Relative fitness levels",labels=c("Least","Most")) + scale_x_continuous(breaks=c(mres,1),"Resistance levels",labels=c("Most","Least"))
p4
ggsave(paste(pref,"norsd_Array_w=",omega4,"_06.pdf",sep=""))


setwd(plots)
pdf(paste(pref,"norsd_Array_w_all.pdf",sep=""),width=18,height=18)
multiplot(p1,p2,p4,p3,cols=2)
dev.off()

## plot U, S & R over time
Mu <- as.data.frame(cbind(seq(0, tsteps, 1),matrix(0,tsteps+1,1), Sv05$U, Sv10$U,Sv15$U, Sv20$U))
Mb <- as.data.frame(cbind(seq(0, tsteps, 1),matrix(1,tsteps+1,1), Sv05$B, Sv10$B,Sv15$B, Sv20$B))
Mmf <- as.data.frame(cbind(seq(0, tsteps, 1),matrix(3,tsteps+1,1), Sv05$meanf[,1], Sv10$meanf[,1],Sv15$meanf[,1], Sv20$meanf[,1]))
Mmr <- as.data.frame(cbind(seq(0, tsteps, 1),matrix(4,tsteps+1,1), Sv05$meanf[,2], Sv10$meanf[,2],Sv15$meanf[,2], Sv20$meanf[,2]))
Mhigh <- as.data.frame(cbind(seq(0, tsteps, 1),matrix(5,tsteps+1,1), Sv05$M[5,5,], Sv10$M[5,5,], Sv15$M[5,5,], Sv20$M[5,5,]))
Musr <- rbind(Mu, Mb,Mmf,Mmr,Mhigh)
colnames(Musr) <- c("t", "pop",4,3,2,1)
Msrm <- melt(Musr, id.vars = c("t","pop"))
facet_names <- c(`0` = "U", `1` = "B", `3` = "mean fit", `4` = "mean res", `5` = "Highest fit/res")
pmb<-ggplot(Msrm, aes(x=t, y = value, colour = variable)) + geom_line() + facet_wrap(~pop,labeller = as_labeller(facet_names), scales = "free")
ggsave(paste(pref,"norsd_TimeSeries_output_06.pdf",sep=""))
# number in highest fitness changes but mean r and f don't? 

# How is the bacterial load affected? want to control.
colnames(Mb) <- c("t", "pop",4,3,2,1)
Mbm <- melt(Mb, id.vars = c("t","pop"))
ggplot(Mbm, aes(x=t, y = value, colour = variable)) + geom_line() + scale_x_continuous(lim = c(0,1200)) + scale_y_continuous(lim = c(0,100))
ggsave(paste(pref,"norsd_TimeSeries_Boutput_06.pdf",sep=""))

# plots only proportion with "resistance" - not background initial resistance
M20n <- Sv20$M;M10n <- Sv10$M;M15n <- Sv15$M;M05n <- Sv05$M
M20n[susr,susf,] <- 0; M10n[susr,susf,] <- 0; M05n[susr,susf,] <- 0; M15n[susr,susf,] <- 0; 
for (i in 1:dim(Sv20$M)[3]){ M20n[,,i] <- M20n[,,i]/sum(M20n[,,i]) # re-normalise
M15n[,,i] <- M15n[,,i]/sum(M15n[,,i]) # re-normalise
M10n[,,i] <- M10n[,,i]/sum(M10n[,,i]) # re-normalise
M05n[,,i] <- M05n[,,i]/sum(M05n[,,i])} # re-normalise

# What happens? 
mm20<-c() ; mm10<-c() ; mm05<-c() ; mm15<-c() 
ll<-dim(Sv20$M)[3];
ss<-seq(0,ll,1/dt) # Don't want to grab all 
for(i in 2:length(ss)){
  mm220<-as.data.frame(melt(M20n[,,ss[i]])); 
  mm215<-as.data.frame(melt(M15n[,,ss[i]])); 
  mm210<-as.data.frame(melt(M10n[,,ss[i]])); 
  mm205<-as.data.frame(melt(M05n[,,ss[i]])); 
  mm220$tstep=ss[i]*dt; mm215$tstep=ss[i]*dt; mm210$tstep=ss[i]*dt; mm205$tstep=ss[i]*dt # To go to generations
  mm20<-rbind(mm20,mm220);mm15<-rbind(mm15,mm215); mm10<-rbind(mm10,mm210) ; mm05<-rbind(mm05,mm205) 
} 
colnames(mm20)<-c("x","y","z","t"); colnames(mm15)<-c("x","y","z","t"); colnames(mm10)<-c("x","y","z","t") ;colnames(mm05)<-c("x","y","z","t")

# plots
p1<-ggplot(mm20[w,],aes(x,y,fill=z))  + facet_wrap( ~ t, ncol=3) + ggtitle(paste("w = ", 1,sep=""))
p1<-p1 + scale_fill_gradient("Proportion", limits=c(0,1),low="white", high="red",guide = FALSE)
p1<-p1 + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p1<-p1 + geom_tile() + scale_y_continuous(breaks=c(1,nfit),"Relative fitness levels",labels=c("Least","Most")) + scale_x_continuous(breaks=c(mres,1),"Resistance levels",labels=c("Most","Least"))
p1
ggsave(paste(pref,"norsdR_Array_w=",omega1,"_06.pdf",sep=""))
p2<-ggplot(mm15[w,],aes(x,y,fill=z))  + facet_wrap( ~ t, ncol=3) + ggtitle(paste("w = ", 2,sep=""))
p2<-p2 + scale_fill_gradient("Proportion", limits=c(0,1),low="white", high="red",guide=FALSE)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p2<-p2 + geom_tile() + scale_y_continuous(breaks=c(1,nfit),"Relative fitness levels",labels=c("Least","Most")) + scale_x_continuous(breaks=c(mres,1),"Resistance levels",labels=c("Most","Least"))
p2
ggsave(paste(pref,"norsdR_Array_w=",omega2,"_06.pdf",sep=""))
p3<-ggplot(mm10[w,],aes(x,y,fill=z))  + facet_wrap( ~ t, ncol=3) + ggtitle(paste("w = ", 3,sep=""))
p3<-p3 + scale_fill_gradient("Proportion", limits=c(0,1),low="white", high="red",guide = FALSE)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p3<-p3 + geom_tile() + scale_y_continuous(breaks=c(1,nfit),"Relative fitness levels",labels=c("Least","Most")) + scale_x_continuous(breaks=c(mres,1),"Resistance levels",labels=c("Most","Least"))
p3
ggsave(paste(pref,"norsdR_Array_w=",omega3,"_06.pdf",sep=""))
p4<-ggplot(mm05[w,],aes(x,y,fill=z))  + facet_wrap( ~ t, ncol=3) + ggtitle(paste("w = ", 4,sep=""))
p4<-p4 + scale_fill_gradient("Proportion", limits=c(0,1),low="white", high="red",guide=FALSE)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p4<-p4 + geom_tile() + scale_y_continuous(breaks=c(1,nfit),"Relative fitness levels",labels=c("Least","Most")) + scale_x_continuous(breaks=c(mres,1),"Resistance levels",labels=c("Most","Least"))
p4
ggsave(paste(pref,"norsdR_Array_w=",omega4,"_06.pdf",sep=""))

setwd(plots)
pdf(paste(pref,"norsdR_Array_w_all.pdf",sep=""),width=18,height=18)
multiplot(p1,p2,p4,p3,cols=2)
dev.off()

pdf(paste(pref,"all_norsd_Array_w_all.pdf",sep = ""),width=18,height=18)
multiplot(go,p1,p2,pmb,p4,p3,cols=2)
dev.off()



#############********************************************** HIGH SHORT OR LOW LONG *********************************************########
## prefix = HSoLL
pref <- "HSoLL_"
dt=0.1
tsteps<-500*(1/dt)
omega1 <- c(matrix(30,1,10*(1/dt)),matrix(0,1,tsteps - 10*(1/dt)))
omega2 <- c(matrix(30,1,20*(1/dt)),matrix(0,1,tsteps - 20*(1/dt)))
omega3 <- c(matrix(0.4,1,50*(1/dt)),matrix(0,1,tsteps - 50*(1/dt)))
omega4 <- c(matrix(0.4,1,100*(1/dt)),matrix(0,1,tsteps - 100*(1/dt)))
omegam <- as.data.frame(cbind(seq(1,tsteps),omega1,omega2,omega3, omega4)); colnames(omegam) <- c("time","1","2","3","4"); momegam <- melt(omegam,id.vars = "time")
go<-ggplot(momegam,aes(x=time,y=value,colour=variable)) + geom_line(size=2) + scale_x_continuous(lim=c(0,200), breaks = seq(0,200,50), labels = dt*seq(0,200,50)) + scale_y_continuous("Omega") + scale_color_discrete("")
ggsave(paste(pref,"norsd_omegaconcs_06.pdf",sep=""))

kk <- 500
Sv20<-ec_srs_funcf_mean_varsr(tsteps,home, c(kk,omega1),iniv,M0,acqdistn,dt)
Sv15<-ec_srs_funcf_mean_varsr(tsteps,home, c(kk,omega2),iniv,M0,acqdistn,dt)
Sv10<-ec_srs_funcf_mean_varsr(tsteps,home, c(kk,omega3),iniv,M0,acqdistn,dt)
Sv05<-ec_srs_funcf_mean_varsr(tsteps,home, c(kk,omega4),iniv,M0,acqdistn,dt)
## NEED TO SPEED IT UP?? Fast for 5 x 5... ~6 sec on laptop 

# What happens? 
mm20<-c() ; mm10<-c() ; mm05<-c() ; mm15<-c() 
ll<-dim(Sv20$M)[3];
ss<-seq(0,ll,1/dt) # Don't want to grab all 
for(i in 2:length(ss)){
  mm220<-as.data.frame(melt(Sv20$M[,,ss[i]])); 
  mm215<-as.data.frame(melt(Sv15$M[,,ss[i]])); 
  mm210<-as.data.frame(melt(Sv10$M[,,ss[i]])); 
  mm205<-as.data.frame(melt(Sv05$M[,,ss[i]])); 
  mm220$tstep=ss[i]*dt; mm215$tstep=ss[i]*dt; mm210$tstep=ss[i]*dt; mm205$tstep=ss[i]*dt # To go to generations
  mm20<-rbind(mm20,mm220);mm15<-rbind(mm15,mm215); mm10<-rbind(mm10,mm210) ; mm05<-rbind(mm05,mm205) 
} 
colnames(mm20)<-c("x","y","z","t"); colnames(mm15)<-c("x","y","z","t"); colnames(mm10)<-c("x","y","z","t") ;colnames(mm05)<-c("x","y","z","t")
#mm20$x<-seq(mres,1,-1); mm10$x<-seq(mres,1,-1); mm05$x<-seq(mres,1,-1)
setwd(plots)
# Grab a subset
#sub<-c(1/dt,250,500,750,1000,1500,2000,2500,seq(3000,4001,1000),4500,tsteps)*dt
sub<-c(1/dt,100,250,500,750,1000,1500,2000,2500,seq(3000,4001,500))*dt
w<-which(mm20[,"t"] %in% sub)
# plots
p1<-ggplot(mm20[w,],aes(x,y,fill=z))  + facet_wrap( ~ t, ncol=3) + ggtitle(paste("w = ", 1,sep=""))
p1<-p1 + scale_fill_gradient("Proportion", limits=c(0,1),low="white", high="red",guide = FALSE)
p1<-p1 + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p1<-p1 + geom_tile() + scale_y_continuous(breaks=c(1,nfit),"Relative fitness levels",labels=c("Least","Most")) + scale_x_continuous(breaks=c(mres,1),"Resistance levels",labels=c("Most","Least"))
p1
ggsave(paste(pref,"norsd_Array_w=",omega1,"_06.pdf",sep=""))
p2<-ggplot(mm15[w,],aes(x,y,fill=z))  + facet_wrap( ~ t, ncol=3) + ggtitle(paste("w = ", 2,sep=""))
p2<-p2 + scale_fill_gradient("Proportion", limits=c(0,1),low="white", high="red",guide=FALSE)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p2<-p2 + geom_tile() + scale_y_continuous(breaks=c(1,nfit),"Relative fitness levels",labels=c("Least","Most")) + scale_x_continuous(breaks=c(mres,1),"Resistance levels",labels=c("Most","Least"))
p2
ggsave(paste(pref,"norsd_Array_w=",omega2,"_06.pdf",sep=""))
p3<-ggplot(mm10[w,],aes(x,y,fill=z))  + facet_wrap( ~ t, ncol=3) + ggtitle(paste("w = ", 3,sep=""))
p3<-p3 + scale_fill_gradient("Proportion", limits=c(0,1),low="white", high="red",guide = FALSE)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p3<-p3 + geom_tile() + scale_y_continuous(breaks=c(1,nfit),"Relative fitness levels",labels=c("Least","Most")) + scale_x_continuous(breaks=c(mres,1),"Resistance levels",labels=c("Most","Least"))
p3
ggsave(paste(pref,"norsd_Array_w=",omega3,"_06.pdf",sep=""))
p4<-ggplot(mm05[w,],aes(x,y,fill=z))  + facet_wrap( ~ t, ncol=3) + ggtitle(paste("w = ", 4,sep=""))
p4<-p4 + scale_fill_gradient("Proportion", limits=c(0,1),low="white", high="red",guide=FALSE)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p4<-p4 + geom_tile() + scale_y_continuous(breaks=c(1,nfit),"Relative fitness levels",labels=c("Least","Most")) + scale_x_continuous(breaks=c(mres,1),"Resistance levels",labels=c("Most","Least"))
p4
ggsave(paste(pref,"norsd_Array_w=",omega4,"_06.pdf",sep=""))


setwd(plots)
pdf(paste(pref,"norsd_Array_w_all.pdf",sep=""),width=18,height=18)
multiplot(p1,p2,p4,p3,cols=2)
dev.off()

## plot U, S & R over time
Mu <- as.data.frame(cbind(seq(0, tsteps, 1),matrix(0,tsteps+1,1), Sv05$U, Sv10$U,Sv15$U, Sv20$U))
Mb <- as.data.frame(cbind(seq(0, tsteps, 1),matrix(1,tsteps+1,1), Sv05$B, Sv10$B,Sv15$B, Sv20$B))
Mmf <- as.data.frame(cbind(seq(0, tsteps, 1),matrix(3,tsteps+1,1), Sv05$meanf[,1], Sv10$meanf[,1],Sv15$meanf[,1], Sv20$meanf[,1]))
Mmr <- as.data.frame(cbind(seq(0, tsteps, 1),matrix(4,tsteps+1,1), Sv05$meanf[,2], Sv10$meanf[,2],Sv15$meanf[,2], Sv20$meanf[,2]))
Mhigh <- as.data.frame(cbind(seq(0, tsteps, 1),matrix(5,tsteps+1,1), Sv05$M[5,5,], Sv10$M[5,5,], Sv15$M[5,5,], Sv20$M[5,5,]))
Musr <- rbind(Mu, Mb,Mmf,Mmr,Mhigh)
colnames(Musr) <- c("t", "pop",4,3,2,1)
Msrm <- melt(Musr, id.vars = c("t","pop"))
facet_names <- c(`0` = "U", `1` = "B", `3` = "mean fit", `4` = "mean res", `5` = "Highest fit/res")
ggplot(Msrm, aes(x=t, y = value, colour = variable)) + geom_line() + facet_wrap(~pop,labeller = as_labeller(facet_names), scales = "free")
ggsave(paste(pref,"norsd_TimeSeries_output_06.pdf",sep=""))
# number in highest fitness changes but mean r and f don't? 

# How is the bacterial load affected? want to control.
colnames(Mb) <- c("t", "pop",4,3,2,1)
Mbm <- melt(Mb, id.vars = c("t","pop"))
ggplot(Mbm, aes(x=t, y = value, colour = variable)) + geom_line() + scale_x_continuous(lim = c(0,1200)) + scale_y_continuous(lim = c(0,100))
ggsave(paste(pref,"norsd_TimeSeries_Boutput_06.pdf",sep=""))

# plots only proportion with "resistance" - not background initial resistance
M20n <- Sv20$M;M10n <- Sv10$M;M15n <- Sv15$M;M05n <- Sv05$M
M20n[susr,susf,] <- 0; M10n[susr,susf,] <- 0; M05n[susr,susf,] <- 0; M15n[susr,susf,] <- 0; 
for (i in 1:dim(Sv20$M)[3]){ M20n[,,i] <- M20n[,,i]/sum(M20n[,,i]) # re-normalise
M15n[,,i] <- M15n[,,i]/sum(M15n[,,i]) # re-normalise
M10n[,,i] <- M10n[,,i]/sum(M10n[,,i]) # re-normalise
M05n[,,i] <- M05n[,,i]/sum(M05n[,,i])} # re-normalise

# What happens? 
mm20<-c() ; mm10<-c() ; mm05<-c() ; mm15<-c() 
ll<-dim(Sv20$M)[3];
ss<-seq(0,ll,1/dt) # Don't want to grab all 
for(i in 2:length(ss)){
  mm220<-as.data.frame(melt(M20n[,,ss[i]])); 
  mm215<-as.data.frame(melt(M15n[,,ss[i]])); 
  mm210<-as.data.frame(melt(M10n[,,ss[i]])); 
  mm205<-as.data.frame(melt(M05n[,,ss[i]])); 
  mm220$tstep=ss[i]*dt; mm215$tstep=ss[i]*dt; mm210$tstep=ss[i]*dt; mm205$tstep=ss[i]*dt # To go to generations
  mm20<-rbind(mm20,mm220);mm15<-rbind(mm15,mm215); mm10<-rbind(mm10,mm210) ; mm05<-rbind(mm05,mm205) 
} 
colnames(mm20)<-c("x","y","z","t"); colnames(mm15)<-c("x","y","z","t"); colnames(mm10)<-c("x","y","z","t") ;colnames(mm05)<-c("x","y","z","t")

# plots
p1<-ggplot(mm20[w,],aes(x,y,fill=z))  + facet_wrap( ~ t, ncol=3) + ggtitle(paste("w = ", 1,sep=""))
p1<-p1 + scale_fill_gradient("Proportion", limits=c(0,1),low="white", high="red",guide = FALSE)
p1<-p1 + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p1<-p1 + geom_tile() + scale_y_continuous(breaks=c(1,nfit),"Relative fitness levels",labels=c("Least","Most")) + scale_x_continuous(breaks=c(mres,1),"Resistance levels",labels=c("Most","Least"))
p1
ggsave(paste(pref,"norsdR_Array_w=",omega1,"_06.pdf",sep=""))
p2<-ggplot(mm15[w,],aes(x,y,fill=z))  + facet_wrap( ~ t, ncol=3) + ggtitle(paste("w = ", 2,sep=""))
p2<-p2 + scale_fill_gradient("Proportion", limits=c(0,1),low="white", high="red",guide=FALSE)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p2<-p2 + geom_tile() + scale_y_continuous(breaks=c(1,nfit),"Relative fitness levels",labels=c("Least","Most")) + scale_x_continuous(breaks=c(mres,1),"Resistance levels",labels=c("Most","Least"))
p2
ggsave(paste(pref,"norsdR_Array_w=",omega2,"_06.pdf",sep=""))
p3<-ggplot(mm10[w,],aes(x,y,fill=z))  + facet_wrap( ~ t, ncol=3) + ggtitle(paste("w = ", 3,sep=""))
p3<-p3 + scale_fill_gradient("Proportion", limits=c(0,1),low="white", high="red",guide = FALSE)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p3<-p3 + geom_tile() + scale_y_continuous(breaks=c(1,nfit),"Relative fitness levels",labels=c("Least","Most")) + scale_x_continuous(breaks=c(mres,1),"Resistance levels",labels=c("Most","Least"))
p3
ggsave(paste(pref,"norsdR_Array_w=",omega3,"_06.pdf",sep=""))
p4<-ggplot(mm05[w,],aes(x,y,fill=z))  + facet_wrap( ~ t, ncol=3) + ggtitle(paste("w = ", 4,sep=""))
p4<-p4 + scale_fill_gradient("Proportion", limits=c(0,1),low="white", high="red",guide=FALSE)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p4<-p4 + geom_tile() + scale_y_continuous(breaks=c(1,nfit),"Relative fitness levels",labels=c("Least","Most")) + scale_x_continuous(breaks=c(mres,1),"Resistance levels",labels=c("Most","Least"))
p4
ggsave(paste(pref,"norsdR_Array_w=",omega4,"_06.pdf",sep=""))

setwd(plots)
pdf(paste(pref,"norsdR_Array_w_all.pdf",sep=""),width=18,height=18)
multiplot(p1,p2,p4,p3,cols=2)
dev.off()

pdf(paste(pref,"all_norsd_Array_w_all.pdf",sep = ""),width=18,height=18)
multiplot(go,p1,p2,pmb,p4,p3,cols=2)
dev.off()
