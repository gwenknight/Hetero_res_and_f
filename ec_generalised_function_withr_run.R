### Run simple e coli model

##### Code to run generalised function and examples

##*** Libraries needed
library(mvtnorm);library(plyr); library(ggplot2);library(reshape2);library(deSolve);library(grid);library(gtools); library(directlabels); library(mvtnorm)
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


# Initial conditions
iniv<-c(98,1,1)

#############********************************************** LOAD UP TO HERE ********************************************** 
dt=0.1
tsteps<-500*(1/dt)
omega1 <- 40
omega2 <- 20
omega3 <- 5
Sv20<-ec_funcf_mean_varsr(tsteps,home, c(omega1),iniv,M0,acqdistn,dt,500)
Sv10<-ec_funcf_mean_varsr(tsteps,home, c(omega2),iniv,M0,acqdistn,dt,500)
Sv05<-ec_funcf_mean_varsr(tsteps,home, c(omega3),iniv,M0,acqdistn,dt,500)
## NEED TO SPEED IT UP?? Fast for 5 x 5... ~6 sec on laptop 

# What happens? 
mm20<-c() ; mm10<-c() ; mm05<-c() 
ll<-dim(Sv20$M)[3];
ss<-seq(0,ll,1/dt) # Don't want to grab all 
for(i in 2:length(ss)){
  mm220<-as.data.frame(melt(Sv20$M[,,ss[i]])); 
  mm210<-as.data.frame(melt(Sv10$M[,,ss[i]])); 
  mm205<-as.data.frame(melt(Sv05$M[,,ss[i]])); 
  mm220$tstep=ss[i]*dt; mm210$tstep=ss[i]*dt; mm205$tstep=ss[i]*dt # To go to generations
  mm20<-rbind(mm20,mm220); mm10<-rbind(mm10,mm210) ; mm05<-rbind(mm05,mm205) 
} 
colnames(mm20)<-c("x","y","z","t"); colnames(mm10)<-c("x","y","z","t") ;colnames(mm05)<-c("x","y","z","t")
#mm20$x<-seq(mres,1,-1); mm10$x<-seq(mres,1,-1); mm05$x<-seq(mres,1,-1)
setwd(plots)
# Grab a subset
#sub<-c(1/dt,250,500,750,1000,1500,2000,2500,seq(3000,4001,1000),4500,tsteps)*dt
sub<-c(1/dt,100,250,500,750,1000,1500,2000,2500,seq(3000,4001,500))*dt
w<-which(mm20[,"t"] %in% sub)
# plots
p1<-ggplot(mm20[w,],aes(x,y,fill=z))  + facet_wrap( ~ t, ncol=3) + ggtitle(paste("w = ", omega1,sep=""))
p1<-p1 + scale_fill_gradient("Proportion", limits=c(0,1),low="white", high="red",guide = FALSE)
p1<-p1 + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p1<-p1 + geom_tile() + scale_y_continuous(breaks=c(1,nfit),"Relative fitness levels",labels=c("Least","Most")) + scale_x_continuous(breaks=c(mres,1),"Resistance levels",labels=c("Most","Least"))
p1
ggsave(paste("Array_w=",omega1,"_06.pdf",sep=""))
p9<-ggplot(mm10[w,],aes(x,y,fill=z))  + facet_wrap( ~ t, ncol=3) + ggtitle(paste("w = ", omega2,sep=""))
p9<-p9 + scale_fill_gradient("Proportion", limits=c(0,1),low="white", high="red",guide = FALSE)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p9<-p9 + geom_tile() + scale_y_continuous(breaks=c(1,nfit),"Relative fitness levels",labels=c("Least","Most")) + scale_x_continuous(breaks=c(mres,1),"Resistance levels",labels=c("Most","Least"))
p9
ggsave(paste("Array_w=",omega2,"_06.pdf",sep=""))
p4<-ggplot(mm05[w,],aes(x,y,fill=z))  + facet_wrap( ~ t, ncol=3) + ggtitle(paste("w = ", omega3,sep=""))
p4<-p4 + scale_fill_gradient("Proportion", limits=c(0,1),low="white", high="red",guide=FALSE)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p4<-p4 + geom_tile() + scale_y_continuous(breaks=c(1,nfit),"Relative fitness levels",labels=c("Least","Most")) + scale_x_continuous(breaks=c(mres,1),"Resistance levels",labels=c("Most","Least"))
p4
ggsave(paste("Array_w=",omega3,"_06.pdf",sep=""))

setwd(plots)
pdf("Array_w_all.pdf",width=12,height=8)
multiplot(p1,p9,cols=2)
dev.off()







### Look at proportion in each of the 30 levels over time for each - facet = level
mm20$omega = omega1; mm10$omega = omega2; mm05$omega = omega3
mega<-as.data.frame(rbind(mm20,mm10,mm05)); colnames(mega)<-c("x","y","z","time","omega")
mega$level = c(seq(21,25,1),seq(16,20,1),seq(11,15,1),seq(6,10,1),seq(1,5,1))
g<-ggplot(mega,aes(x=time,y=z,colour=factor(omega))) + geom_line(size=2) + facet_wrap( ~ level, ncol=5) +  scale_colour_manual(values=cbPalette,"Abx\nLevel",breaks=c(omega1,omega2, omega3),labels=c("0.2","0.1","0.05"))
g<-g + scale_x_continuous("Generations",breaks=c(0,200,400)) + scale_y_continuous("Proportion at this level",breaks=c(0.25,0.75))+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
g
setwd(plots)
ggsave("mega.pdf",width=12,height=8)

theme_set(theme_bw(base_size = 14))
g<-g+ facet_wrap( ~ level,scales = "free", ncol=5)+ scale_y_continuous("Proportion at this level")
ggsave("mega_freescale.pdf")

### Look at change in R & S over time
theme_set(theme_bw(base_size = 34))
rrmr<-as.data.frame(cbind(seq(1,ll,1)*dt,Sv05$R,Sv10$R,Sv20$R,1)); colnames(rrmr)<-c("time","05","10","20","type")
rrms<-as.data.frame(cbind(seq(1,ll,1)*dt,Sv05$S,Sv10$S,Sv20$S,2)); colnames(rrms)<-c("time","05","10","20","type")
rrm<-as.data.frame(rbind(rrmr,rrms))
rrm2<-melt(rrm,id.vars=c("time","type")); rrm2[which(rrm2$type==1),"type"]<-"Resistant"; rrm2[which(rrm2$type==2),"type"]<-"Susceptible"
g<-ggplot(rrm2,aes(x=time,y=value,colour=factor(variable))) + geom_line(size=2) + scale_x_continuous("Generations") + scale_y_continuous("Percentage with R")
g<-g + scale_colour_manual("Abx\nLevel",breaks=c("20","10","05"),labels=c(omega1, omega2,omega3),values = cbPalette) + facet_wrap(~type)
g
ggsave("r&s_overtime.pdf",width=12,height=8)

# Time to dominance...
#t05<-min(intersect(intersect(which(rrmr[,2]>79.99),which(rrmr[,2]< 80.007)),which(floor(rrmr[,2])==80))*dt)
#t10<-min(intersect(intersect(which(rrmr[,3]>79.99),which(rrmr[,3]< 80.002)),which(floor(rrmr[,3])==80))*dt)
#t20<-min(intersect(intersect(which(rrmr[,4]>79.99),which(rrmr[,4]< 80.007)),which(floor(rrmr[,4])==80))*dt)

t05 <- rrmr[min(which(rrmr[,2] > 50)),"time"]
t10 <- rrmr[min(which(rrmr[,3] > 50)),"time"]
t20 <- rrmr[min(which(rrmr[,4] > 50)),"time"]

mm20_2<-mm20;mm10_2<-mm10;mm05_2<-mm05
mm20_2$t<-mm20_2$t/t20;mm10_2$t<-mm10_2$t/t10;# mm05_2$t<-mm05_2$t/t05 no t05 at moment
mega_2<-as.data.frame(rbind(mm20_2,mm10_2)) #,mm05_2)); 
colnames(mega_2)<-c("x","y","z","time","omega")
mega_2$level = c(seq(21,25,1),seq(16,20,1),seq(11,15,1),seq(6,10,1),seq(1,5,1))
g<-ggplot(mega_2,aes(x=time,y=z,colour=factor(omega))) + geom_line(size=2) + facet_wrap( ~ level, ncol=5, scales = "free") + scale_colour_manual("Abx\nLevel",breaks=c(20,10,5),labels=c("0.2","0.1","0.05"),values = cbPalette)
g<-g + scale_x_continuous("Time to full resistance") + scale_y_continuous("Proportion at this level")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
g
setwd(plots)
ggsave("mega_normtodom.pdf")

# Time to full resistance and fitness
w<-which(mega$level==5); 
t05a<-min(intersect(which(mega[w,"omega"]==omega1),which(mega[w,"z"]>0.5)))
t10a<-min(intersect(which(mega[w,"omega"]==omega2),which(mega[w,"z"]>0.4)))
t20a<-min(intersect(which(mega[w,"omega"]==omega3),which(mega[w,"z"]>0.4)))

#t05<-mega[w[t05a],"time"]; 
t10<-mega[w[t10a],"time"]; t20<-mega[w[t20a],"time"]

mm20_2<-mm20;mm10_2<-mm10;#mm05_2<-mm05
mm20_2$t<-mm20_2$t/t20;mm10_2$t<-mm10_2$t/t10;#mm05_2$t<-mm05_2$t/t05
mega_2<-as.data.frame(rbind(mm20_2,mm10_2,mm05_2)); colnames(mega_2)<-c("x","y","z","time","omega")
mega_2$level = c(seq(21,25,1),seq(16,20,1),seq(11,15,1),seq(6,10,1),seq(1,5,1))
theme_set(theme_bw(base_size = 12)); 
g<-ggplot(mega_2,aes(x=time,y=z,colour=factor(omega))) + geom_line(size=2) + facet_wrap( ~ level, ncol=5, scales = "free") + scale_colour_manual("Abx\nLevel",breaks=c(omega1, omega2, omega3),values = cbPalette)
g<-g + scale_x_continuous("Time to full resistance") + scale_y_continuous("Proportion at this level")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
g
setwd(plots)
ggsave("mega_normtofullR.pdf",width=12,height=8)

# Plot proportions in each fitness / resistance level over time
pp<-c();
ll<-dim(Sv20$M)[3];
ss<-seq(0,ll,1/dt) # Don't want to grab all 
for(i in 2:length(ss)){
  pp220<-c(ss[i]*dt,colSums(Sv20$M[,,ss[i]]), rowSums(Sv20$M[,,ss[i]]),20)
  pp210<-c(ss[i]*dt,colSums(Sv10$M[,,ss[i]]), rowSums(Sv10$M[,,ss[i]]),10)
  pp205<-c(ss[i]*dt,colSums(Sv05$M[,,ss[i]]), rowSums(Sv05$M[,,ss[i]]),5)
  pp<-rbind(pp,pp220,pp210,pp205);
} 
pp<-as.data.frame(pp);colnames(pp)<-c("t","Fitness level 1\n(low)","Fitness level 2","Fitness level 3","Fitness level 4","Fitness level 5\n(high)","Res. level 1\n(low)","Res. level 2","Res. level 3","Res. level 4","Res. level 5\n(high)","w"); 
pp2<-melt(pp,id.vars = c("t","w"))
theme_set(theme_bw(base_size = 34)); 
g<-ggplot(pp2,aes(x=t,y=value,colour=factor(w))) + facet_wrap(~variable,ncol=5) + geom_line(size=2) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
g<-g + scale_x_continuous("Generation") + scale_y_continuous("Proportion") + scale_colour_manual(values=cbPalette,"Abx\nlevel",labels=c(0.05,0.1,0.2))
g # Suggests that although v similar proportions in the most fit fewer are in the higher resistance levels with low level antibiotics use. In fact with this model
# the same rate of selection for no cost mutations is seen whether there is high or low anitbiotic use 
setwd(plots)
ggsave("f&r_overtime.pdf",width=18,height=12)

#### Compare with and without fitness and resistance levels. 
### Range of omega
setwd(home)
omegav<-c(0.2,0.1,0.05) * dt
#omegav<-c(0.3,0.4,0.2) * dt
para<-read.csv("para_ecoli.csv",header=TRUE,check.names=F,stringsAsFactors = FALSE)[,1:2]
for(i in 1:length(para[,1])){assign(para[i,1],para[i,2])}
# Correct for timestep
mu<-mu*dt;beta<-beta*dt;eps<-eps*dt
assign("ks",1) # Prescribe effect treatment success levels
m<-dim(acqdistn)[1]; vs<-seq(1/m,1,1/m); vf=vs
assign("kr",1-sum(colSums(acqdistn)*vs)) # Set the levels to be those of the acqdistn
assign("f",sum(colSums(acqdistn)*vf))
## SAME as writeup_ecoli
#kr = 0.4; f = 0.6 ## 40% cost to both
bigall<-c(); lambdasv<-c(); lambdarv<-c(); 
endp<-200*1/dt
U<-matrix(0,1,endp); S<-matrix(0,1,endp); R<-matrix(0,1,endp);
U[1]<-iniv[1]; S[1]<-iniv[2]; R[1]<-iniv[3];
lambdasv<-matrix(0,1,endp);lambdarv<-matrix(0,1,endp);
lambdasv[1] = beta * S[1]/sum(iniv);   lambdarv[1] = sum(colSums(acqdistn*seq(1/nfit,1,1/nfit))) * beta * R[1]/sum(iniv) # function outputs just meanfit when all popns 0
setwd(home) # might have para in different place for different models 
for(j in 1:length(omegav)){
  assign("omega",omegav[j])
  for(i in 1:endp){
    lambdas=lambdasv[i];lambdar=lambdarv[i];
    # Dynamics
    U[i+1] =  U[i] + mu*(S[i]+R[i]) - (lambdas+lambdar)*U[i] + omega*ks*S[i] + omega*kr*R[i]
    S[i+1] =  S[i] + lambdas*U[i] - (mu + omega*ks)*S[i] - eps * S[i]
    R[i+1] =  R[i] + lambdar*U[i] - (mu + omega*kr)*R[i] + eps * S[i] 
    lambdasv[i+1] = beta * S[i+1] / 100; 
    lambdarv[i+1] = f * beta * R[i+1] / 100;   #X$meanfit * beta * R[i+1]/N; 
  } 
  all<-as.data.frame(cbind(seq(0,endp,1)*dt,U,S,R,omega/dt)); colnames(all)<-c("time","U","Susceptible","Resistant","w")
  bigall<-rbind(bigall,all)
}
allm<-melt(bigall[,c("time","Susceptible","Resistant","w")], id.vars=c("w","time"))
allm$nw = allm$w
theme_set(theme_bw(base_size = 34))
colnames(rrm2)<-c("time","variable","w","value")
rrm2n<-rrm2[,c("w","time","variable","value")]; 
rrm2n$with<-1; rrm2n$nw<-0
rrm2n[which(rrm2n$w=="05"),"nw"]=allm[12000,"nw"]; 
rrm2n[which(rrm2n$w=="10"),"nw"]<-allm[4000,"nw"]; 
rrm2n[which(rrm2n$w=="20"),"nw"]<-allm[1,"nw"]
allm$with<-0; allm$nw<-allm$w; 
allmn<-rbind(allm,rrm2n)

p<-ggplot(allmn,aes(x=time,y=value,colour=variable,linetype=factor(with)))+geom_line(size=2) + scale_x_continuous("Generations",lim=c(0,endp*dt))
p<-p+scale_colour_discrete("Sub-\npopulation",breaks=c("Susceptible","Resistant")) + scale_y_continuous("Percentage of population") + facet_wrap( ~ nw)  + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p<-p+ scale_linetype_discrete("With\ndiversity",breaks=c(0,1),labels=c("None","With diversity")) 
p
setwd(plots)
ggsave("WithnWithoutdiversity.pdf",width=16,height=8)

p + scale_y_continuous("Percentage of population",lim=c(0,10)) 
ggsave("WithnWithoutdiversity_zoom.pdf",width=16,height=8)

## Plot like in Gulberg
pp2n<-pp2[7501:15000,]
w<-intersect(which(pp2n$w==5),c(which(pp2n$t==100),which(pp2n$t==200),which(pp2n$t==300),which(pp2n$t==400),which(pp2n$t==500)))
ggplot(pp2n[w,], aes(x=t, y= value,colour=factor(variable))) + geom_point(aes(shape = factor(variable)),size=5)

library(data.table)
drrm2n <- data.table(pp2n[w,])
factors <- unique(drrm2n[, variable])
nfactors <- length(factors)
width = 25
for (id in seq_along(factors)) {
  drrm2n[variable == factors[id], adj.time := t - width + (id - 1) * (2*width) / (nfactors - 1)]
}
g<-ggplot(drrm2n, aes(x=adj.time, y= value,colour=factor(variable))) + geom_point(aes(shape = factor(variable)),size=5)+ theme(legend.position="top")
g<-g+ scale_x_continuous("Generations") + scale_y_continuous("Proportion") + scale_colour_discrete("")+ scale_shape_discrete("")
setwd(plots)
ggsave("Gulberg5d.pdf",width=16,height=8)

#############******************************************************************************************** 
#############******************************************************************************************** 
#### change the initial conditions
#############******************************************************************************************** 
#############******************************************************************************************** 
# 2: normalise but mean at 0.2
f <- function(x, y) dmvnorm(cbind(x, y), mean = c(0.2, 0.2),sigma = diag(2)/20)
z <- outer(x, y, FUN = f); z <- z/sum(z) # Generate and normalise
dev.off(); persp(x, y, z, theta = -30, phi = 30, ticktype = "detailed")
acqdistn<-z

plot_diff_acd_output(acqdistn,plots,2)

# 3: flatten to all low resistance but high fitness
f <- function(x, y) dmvnorm(cbind(x, y), mean = c(0.01, 0.01),sigma = diag(2)/20)
z <- outer(x, y, FUN = f); z <- z/sum(z) # Generate and normalise
dev.off(); persp(x, y, z, theta = -30, phi = 30, ticktype = "detailed")
zn<-z
zn[1:5,]<-colSums(z) # Doesn't matter if col or row Sums. 
acqdistn<-zn/sum(zn)

plot_diff_acd_output(acqdistn,plots,3)

# 4: flatten to all low resistance but high fitness
f <- function(x, y) dmvnorm(cbind(x, y), mean = c(0.01, 0.01),sigma = diag(2)/20)
z <- outer(x, y, FUN = f); z <- z/sum(z) # Generate and normalise
dev.off(); persp(x, y, z, theta = -30, phi = 30, ticktype = "detailed")
zn<-z
zn[,1:5]<-colSums(z) # Doesn't matter if col or row Sums. 
acqdistn<-t(zn)/sum(zn)

plot_diff_acd_output(acqdistn,plots,4)

# 5: EG: start at 0.8
x=y=seq(0.2,1,0.2)
f <- function(x, y) dmvnorm(cbind(x, y), mean = c(0.2, 0.2),sigma = diag(2)/20)
z <- outer(x, y, FUN = f); z <- z/sum(z) # Generate and normalise
dev.off(); persp(x, y, z, theta = -30, phi = 30, ticktype = "detailed")

nfit = 10; mres = 10;
# Array of distribution of fitness and resistance c[resistance, fitness]
M0 <- array(0,c(mres,nfit,10))
aa<-M0[,,1]
aa[5,1:5]<-z[1,];aa[4,1:5]<-z[2,];aa[3,1:5]<-z[3,];aa[2,1:5]<-z[4,];aa[1,1:5]<-z[5,]
acqdistn<-aa
x=y=seq(0.1,1,0.1)
dev.off(); persp(x, y, aa, theta = -30, phi = 30, ticktype = "detailed")

plot_diff_acd_output(acqdistn,plots,5)

# 6: use standard initial distribution but start at 0.8
x=y=seq(0.2,1,0.2)
f <- function(x, y) dmvnorm(cbind(x, y), mean = c(0.6, 0.6),sigma = diag(2)/20)
z <- outer(x, y, FUN = f); z <- z/sum(z) # Generate and normalise
dev.off(); persp(x, y, z, theta = -30, phi = 30, ticktype = "detailed")

nfit = 25; mres = 25;
# Array of distribution of fitness and resistance c[resistance, fitness]
M0 <- array(0,c(mres,nfit,10))
aa<-M0[,,1]
aa[20,21:25]<-z[1,];aa[21,21:25]<-z[2,];aa[22,21:25]<-z[3,];aa[23,21:25]<-z[4,];aa[24,21:25]<-z[5,]
acqdistn<-aa
x=y=seq(1/25,1,1/25)
dev.off(); persp(x, y, aa, theta = -30, phi = 30, ticktype = "detailed")

plot_diff_acd_output(acqdistn,plots,6)

