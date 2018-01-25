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
  lambdav[i+1] = max(0,(1-omega)/1) * beta * B[i] 
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
# Step two: 


