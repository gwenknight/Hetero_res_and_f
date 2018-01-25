###### E coli model. 
## From basics of no diversity / treatment etc

# Home
home<-"~/Documents/Hetero_res_and_f/"
plots<-paste(home,"plots",sep="")
setwd(home)
library(mvtnorm);library(plyr); library(ggplot2);library(reshape2);library(deSolve);library(grid);library(gtools); library(directlabels); library(mvtnorm)
#source("ec_norsdiff_generalised_function.R")
theme_set(theme_bw(base_size = 34))
## Key parameters
dt<-0.1 # 1/10th of a generation
endp<-200*1/dt

#************************************************************************************************************************
## Step 1. Fit to get, in the absence of antibiotics, 90% of the population is bacteria at equilibrium
# in absence of antibiotics
# bugs have high relative fitness
para<-read.csv("data/para_ecoli.csv",header=TRUE,check.names=F,stringsAsFactors = FALSE)[,1:2]
for(i in 1:length(para[,1])){assign(para[i,1],para[i,2])}
# Correct for timestep
mu<-mu*dt;beta<-beta*dt;

# Initial conditions
# 99% of the culture is space -> 1% is bug
initial<-c(99,1); N=sum(initial)
U<-matrix(0,1,endp); B<-matrix(0,1,endp); # JUST ONE TYPE OF BUG
U[1]<-initial[1]; B[1]<-initial[2];

lambdav<-matrix(0,1,endp);
lambdav[1]=beta  # no quite right but ok for now? (equilibrium interested in)

assign("omega",0) # No treatment

for(i in 1:endp){
  lambda=lambdav[i];
  
  # NEW Dynamics
  U[i+1] =  U[i] + mu*(B[i]) - lambda*(U[i]/(U[i] + kk))
  B[i+1] =  B[i] + lambda*(U[i]/(U[i] + kk)) - mu*B[i] 
  
  #U[i+1] =  U[i] + mu*(S[i]) - (lambdas)*(U[i]/(U[i] + kk)) 
  #S[i+1] =  S[i] + lambdas*(U[i]/(U[i] + kk)) - mu*S[i]
  
  # Dynamics
  #U[i+1] =  U[i] + mu*(S[i]) - (lambdas)*U[i] + omega*ks*S[i] 
  #S[i+1] =  S[i] + lambdas*U[i] - (mu + omega*ks)*S[i] #- eps * S[i]
  #lambdasv[i+1] = beta * S[i+1] / 100; 
  lambdav[i+1] = max(0,(1-omega)/1) * beta * B[i+1] 
}  

all<-as.data.frame(cbind(seq(0,200,dt),U,B)); colnames(all)<-c("time","U","B")
allm<-melt(all, id.vars="time")
all[20/dt+1,]
p<-ggplot(allm,aes(x=time,y=value,colour=variable))+geom_line(size=2) + scale_x_continuous("Generations",lim=c(0,30)) +
  scale_y_continuous("Percentage of population",limits = c(0,100)) + scale_colour_discrete("",labels = c("Space","Total bug"))
p # Fitted beta and mu to give 80% susceptible at equilibrium. 
setwd(plots)
ggsave("initial.pdf",width=12,height=8)

#************************************************************************************************************************
### Second set some at beginning - should get even level + little extra via eps?
# Same fitness
initial<-c(98,1,1); N=sum(initial)
U<-matrix(0,1,endp); S<-matrix(0,1,endp); R<-matrix(0,1,endp);
U[1]<-initial[1]; S[1]<-initial[2]; R[1]<-initial[3];

lambdasv<-matrix(0,1,endp);lambdarv<-matrix(0,1,endp);
lambdasv[1]=beta * S[1]/N; lambdarv[1]= 1 * beta * R[1]/N ; # Same fitness

for(i in 1:endp){
  lambdas=lambdasv[i];lambdar=lambdarv[i];# krr = krv[i]
  
  # NEW Dynamics
  U[i+1] =  U[i] + mu*(S[i]+R[i]) - (lambdas+lambdar)*(U[i]/(U[i] + kk)) 
  S[i+1] =  S[i] + lambdas*(U[i]/(U[i] + kk)) - mu*S[i] - eps * S[i]
  R[i+1] =  R[i] + lambdar*(U[i]/(U[i] + kk)) - mu*R[i] + eps * S[i] 
  
  lambdasv[i+1] = max(0,(1-omega)/1) * beta * S[i+1] / ( S[i+1] + R[i+1] );
  lambdarv[i+1] = f * max(0,(6-omega)/6) * beta * R[i+1] / ( S[i+1] + R[i+1] );   # resistant strain has an MIC of 6
} 

all<-as.data.frame(cbind(seq(0,endp,1),U,S,R)); colnames(all)<-c("time","U","S","R")
allm<-melt(all, id.vars="time")
p<-ggplot(allm,aes(x=time,y=value,colour=variable))+geom_line(size=4) + scale_x_continuous("Time",lim=c(0,endp))
p # Get very slight more resistance as only additional eps generation
# Checked: When set eps to zero get exactly the same. 

### Second set some at beginning - should get even level + little extra via eps?
omegav<-c(0.1,0.3,0.4,0.42,0.45,0.5,0.8,1,1.5,2,4,6,8) 
assign("f",0.6)
bigall<-c()
for(j in 1:length(omegav)){
  assign("omega",omegav[j])
  for(i in 1:endp){
    lambdas=lambdasv[i];lambdar=lambdarv[i];# krr = krv[i]
    
    # NEW Dynamics
    U[i+1] =  U[i] + mu*(S[i]+R[i]) - (lambdas+lambdar)*(U[i]/(U[i] + kk)) 
    S[i+1] =  S[i] + lambdas*(U[i]/(U[i] + kk)) - mu*S[i] - eps * S[i]
    R[i+1] =  R[i] + lambdar*(U[i]/(U[i] + kk)) - mu*R[i] + eps * S[i] 
    
    lambdasv[i+1] = max(0,(1-omega)/1) * beta * S[i+1] / ( S[i+1] + R[i+1] );
    lambdarv[i+1] = f * max(0,(6-omega)/6) * beta * R[i+1] / ( S[i+1] + R[i+1] );   #
  } 
  all<-as.data.frame(cbind(seq(0,endp,1)*dt,U,S,R,omega)); colnames(all)<-c("time","U","S","R","w")
  bigall<-rbind(bigall,all)
}

allm<-melt(bigall, id.vars=c("w","time"))
p<-ggplot(allm,aes(x=time,y=value,colour=variable))+geom_line(size=4) + scale_x_continuous("Generations",lim=c(0,endp*dt)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p<-p+scale_colour_manual("Sub-\npopulation", values=cbPaletteSR)  + scale_y_continuous("Percentage of population") + facet_wrap( ~ w)
p 
setwd(plots)
ggsave("initial_omega_change.pdf",width=12,height=8) # Change omega above for changing output

# Increase f to 0.9
assign("f",0.9)
bigall<-c()
for(j in 1:length(omegav)){
  assign("omega",omegav[j])
  for(i in 1:endp){
    lambdas=lambdasv[i];lambdar=lambdarv[i];# krr = krv[i]
    
    # NEW Dynamics
    U[i+1] =  U[i] + mu*(S[i]+R[i]) - (lambdas+lambdar)*(U[i]/(U[i] + kk))
    S[i+1] =  S[i] + lambdas*(U[i]/(U[i] + kk)) - mu*S[i] - eps * S[i]
    R[i+1] =  R[i] + lambdar*(U[i]/(U[i] + kk)) - mu*R[i] + eps * S[i] 
    
    lambdasv[i+1] = max(0,(1-omega)/1) * beta * S[i+1] / ( S[i+1] + R[i+1] );
    lambdarv[i+1] = f * max(0,(6-omega)/6) * beta * R[i+1] / ( S[i+1] + R[i+1] );   #
  } 
  all<-as.data.frame(cbind(seq(0,endp,1)*dt,U,S,R,omega)); colnames(all)<-c("time","U","S","R","w")
  bigall<-rbind(bigall,all)
}

allm<-melt(bigall, id.vars=c("w","time"))
p<-ggplot(allm,aes(x=time,y=value,colour=variable))+geom_line(size = 4) + scale_x_continuous("Generations",lim=c(0,endp*dt))
p<-p+scale_colour_manual("Sub-\npopulation", values=cbPaletteSR) + scale_y_continuous("Percentage of population") + facet_wrap( ~ w)
p 
setwd(plots)
ggsave("initial_omega_change_f1.pdf",width=12,height=8) # Change omega above for changing output

#### ALLOW RESOURCE TO INCREASE - continuous feed
bigall<-c()
for(j in 1:length(omegav)){
  assign("omega",omegav[j])
  for(i in 1:endp){
    lambdas=lambdasv[i];lambdar=lambdarv[i];# krr = krv[i]
    
    # NEW Dynamics
    U[i+1] =  U[i] + mu*(S[i]+R[i]) - (lambdas+lambdar)*(U[i]/(U[i] + kk)) + 0.5
    S[i+1] =  S[i] + lambdas*(U[i]/(U[i] + kk)) - mu*S[i] - eps * S[i]
    R[i+1] =  R[i] + lambdar*(U[i]/(U[i] + kk)) - mu*R[i] + eps * S[i] 
    
    lambdasv[i+1] = max(0,(1-omega)/1) * beta * S[i+1] / ( S[i+1] + R[i+1] );
    lambdarv[i+1] = f * max(0,(6-omega)/6) * beta * R[i+1] / ( S[i+1] + R[i+1] );   #
  } 
  all<-as.data.frame(cbind(seq(0,endp,1)*dt,U,S,R,omega)); colnames(all)<-c("time","U","S","R","w")
  bigall<-rbind(bigall,all)
}

allm<-melt(bigall, id.vars=c("w","time"))
p<-ggplot(allm,aes(x=time,y=value,colour=variable))+geom_line(size = 4) + scale_x_continuous("Generations",lim=c(0,endp*dt))
p<-p+scale_colour_manual("Sub-\npopulation", values=cbPaletteSR) + scale_y_continuous("Percentage of population") + facet_wrap( ~ w)
p 
setwd(plots)
ggsave("initial_omega_f1_9_continuous_resource5.pdf",width=12,height=8) # Change omega above for changing output
#### ALLOW RESOURCE TO INCREASE - continuous feed
bigall<-c()
for(j in 1:length(omegav)){
  assign("omega",omegav[j])
  for(i in 1:endp){
    lambdas=lambdasv[i];lambdar=lambdarv[i];# krr = krv[i]
    
    # NEW Dynamics
    U[i+1] =  U[i] + mu*(S[i]+R[i]) - (lambdas+lambdar)*(U[i]/(U[i] + kk)) + 2
    S[i+1] =  S[i] + lambdas*(U[i]/(U[i] + kk)) - mu*S[i] - eps * S[i]
    R[i+1] =  R[i] + lambdar*(U[i]/(U[i] + kk)) - mu*R[i] + eps * S[i] 
    
    lambdasv[i+1] = max(0,(1-omega)/1) * beta * S[i+1] / ( S[i+1] + R[i+1] );
    lambdarv[i+1] = f * max(0,(6-omega)/6) * beta * R[i+1] / ( S[i+1] + R[i+1] );   #
  } 
  all<-as.data.frame(cbind(seq(0,endp,1)*dt,U,S,R,omega)); colnames(all)<-c("time","U","S","R","w")
  bigall<-rbind(bigall,all)
}

allm<-melt(bigall, id.vars=c("w","time"))
p<-ggplot(allm,aes(x=time,y=value,colour=variable))+geom_line(size = 4) + scale_x_continuous("Generations",lim=c(0,endp*dt))
p<-p+scale_colour_manual("Sub-\npopulation", values=cbPaletteSR) + scale_y_continuous("Percentage of population") + facet_wrap( ~ w)
p 
setwd(plots)
ggsave("initial_omega_f1_9_continuous_resource2.pdf",width=12,height=8) # Change omega above for changing output

#************************************************************************************************************************
###?? NOT UPDATED
# Step three match ratio diagrams from Gulberg.
rat_all<-c(); ratbig<-c()
wv<-c(0,seq(0.001,0.04,0.001))
fv=c(0.6,0.8)
for(k in 1:length(fv)){
  assign("f",fv[k])
  rat_all<-c()
  for(j in 1:length(wv)){
    assign("omega",wv[j])
    for(i in 1:endp){
      lambdas=lambdasv[i];lambdar=lambdarv[i];# krr = krv[i]
      
      # Dynamics
      U[i+1] =  U[i] + mu*(S[i]+R[i]) - (lambdas+lambdar)*U[i] + omega*ks*S[i] + omega*kr*R[i]
      S[i+1] =  S[i] + lambdas*U[i] - (mu + omega*ks)*S[i] - eps * S[i]
      R[i+1] =  R[i] + lambdar*U[i] - (mu + omega*kr)*R[i] + eps * S[i] 
      
      lambdasv[i+1] = beta * S[i+1] / 100; 
      lambdarv[i+1] = f * beta * R[i+1] / 100;   #X$meanfit * beta * R[i+1]/N; 
    }  
    rat<-as.data.frame(cbind(seq(0,endp,1)*dt,U,S,R,R/S,omega/dt)); # Remember to convert back timestep
    rat_all<-rbind(rat_all,rat)
  }
  rat_all2<-rat_all
  rat_all2$f<-fv[k]
  ratbig<-as.data.frame(rbind(ratbig,rat_all2))
}
colnames(ratbig)<-c("time","U","S","R","ratio","w","f")
pg<-ggplot(ratbig,aes(x=time,y=ratio,colour=factor(w)))+geom_line() + scale_colour_discrete(guide = guide_legend(reverse=TRUE)) + scale_x_continuous("Generations",lim=c(0,endp*dt))
pg<-pg + geom_hline(y=1) + scale_y_log10("Percentage of population") 
pg

pg<-ggplot(ratbig,aes(x=time,y=ratio,colour=factor(w)))+geom_line() + scale_colour_discrete("Antibiotic\nconcentration",guide = guide_legend(reverse=TRUE),breaks=c(0,0.01,0.05,0.1,0.15,0.2,0.3,0.4)) + scale_x_continuous("Generations",lim=c(0,40))
pg<-pg + geom_hline(y=1) + scale_y_log10("Ratio resistant : susceptible",lim=c(0.00001,10)) + facet_wrap( ~ f)
pg
ggsave("ratio_match.pdf",width=12,height=8)

wp<-c(which(ratbig$time==10),which(ratbig$time==20),which(ratbig$time==30),which(ratbig$time==40))
wp2<-intersect(wp,c(which(ratbig$w==0),which(ratbig$w==0.05),which(ratbig$w==0.15),which(ratbig$w==0.25),which(ratbig$w==0.3)))
pg<-ggplot(ratbig[wp2,],aes(x=time,y=ratio,colour=factor(w)))+geom_point(size=4) + geom_line() + scale_colour_discrete("Antibiotic\nconcentration",guide = guide_legend(reverse=TRUE)) + scale_x_continuous("Generations",lim=c(0,40))
pg<-pg + geom_hline(y=1) + scale_y_log10("Ratio resistant : susceptible",lim=c(0.00001,100)) + facet_wrap( ~ f)
pg
ggsave("ratio_match_points.pdf",width=12,height=8)




#************************************************************************************************************************
# Step four - add in variation


