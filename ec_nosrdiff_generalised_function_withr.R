##### Code for generalised fitness function 
# With the inclusion of different resistance levels
# Resistance = increased ability to survive. 

######*********************************************************************** Main generalised function code  *******#####
## Function to update the mean fitness & MEAN SUSCEPTIBILITY at each timestep
ec_srs_meanfit_varsr=function(bigM,acq,trans,trate,nA,n,m,acqdist){
  # Needs bigM: Array of distribution  [susceptible rows, fitness columns, time]
  #       acq: number of acquisitions of resistant strains
  #       trans: number of transmissions of resistant strains
  #       trate: antibiotic exposure level
  #       nA: number with bacteria
  #       n: number of fitness levels
  #       m: number of resistance levels
  #       acqdist: 2D distribution of new cases
  
  # Vector of fitness values
  vf<-seq(1/n,1,1/n)
  
  # VECTOR of MIC LEVELs => Susceptible has MIC of 1.  
  # ?? 0.1 => 10% resistant, 0.9 = 90% resistant.
  # vs<-seq(1/m,1,1/m)
  # vector of MIC values
  # MAX MIC = 32
  # MIN MIC = 0.1
  vs_mic <- seq(0.1,32,length.out = m) 
  # Vector of resistance levels
  vs_rel <- pmax((vs_mic - trate),0)/vs_mic # constant 
  
  #print(c("baseline","vs_rel",vs_rel,"vf",vf))
  
  # Which column of M? last that is non-zero plus one for this timestep
  sumM<-colSums(colSums(bigM,dims=1)) # Gives a vector of the sums over the array for each timestep
  if(length(which(sumM > 0))){tt<-tail(which(sumM > 0), n = 1) + 1}else{tt=2} # which timestep.
  
  #************************************************ At initial time point
  if(acq==0 && trans==0){meanfit = 0; return(meanfit); meanres = 0; return(meanres);
  #************************************************ Otherwise look at distribution  
  }else{
    M <- bigM[,,tt-1] # Grab the appropriate matrix for this timestep
    M_new <- M;  # New matrices to update as go through
    
    #**** Update M
    ##### 1) New cases - how distribute acquisitions across fitness and resistance? 
    
    #*** How are the acquisitions distributed? Input
    new_a <- acqdist
    
    #*** How are the new growth bugs distributed?
    M_temp<-M_new
    # if none resistant then none grow... 
    Mr<-rowSums(M_temp) # = rowSums(M_new) same # Proportions at each RESISTANCE level
    if(trans>0){ pastmean=sum( Mr*vs_rel ) }else{ pastmean = 0; }
    if(pastmean > 0){ # if resistant then do the following... otherwise no new bugs
      
      # now growth rate dependent on fitness and resistance 
      # FITNESS
      Mf<-colSums(M_new) # Proportions at each FITNESS level
      if(trans>0){ pastmean=sum( Mf*vf ) }else{ pastmean = 1 }
      M_temp<-M_new
      #for(i in 1:length(Mf)){M_temp[,i] = M_temp[,i] * vf[i] / pastmean } # Updated matrix: colSums(M_new) = Mf_new
      M_temp <- t(t(M_temp) * vf / pastmean)
      
      # RESISTANCE
      Mr<-rowSums(M_temp) # = rowSums(M_new) same # Proportions at each RESISTANCE level
      if(trans>0){ pastmean=sum( Mr*vs_rel ) }else{ pastmean = 1 }
      # update M_temp with fitness then resistance
      #for(i in 1:length(Mf)){M_temp[i,] = M_temp[i,] * vs_rel[i] / pastmean } # Updated matrix: colSums(M_new) = Mf_new
      if(pastmean > 0){M_temp <- vs_rel*M_temp/pastmean} 
      # if some of the 
      # As some may be zero (if antibiotic use too high - remember this is new growth)
      # M_temp <- M_temp/sum(M_temp) ## don't need as sum(vs_rel) now = 1.
      #?? does it matter to do fitness then resistance? don't think it makes a difference. 
      
      # if no new transmission then 
      new_b = M_temp 
    }else{new_b = 0; trans = 0} # if all resistant then no growth
    
    
    #*** Assign and update M
    # nA = number of actives left after treatment removed (affects distribution) and death (doesn't affect)
    #if(nA>0){M[,tt] = M[,tt-1]*( nA/(nA+acq+trans+react) ) + new * (acq+trans+react)/(nA+acq+trans+react)  ## PREVIOUS WITHOUT SUS
    # if(nA>0){M_new = M_new * nA / (nA+acq+trans) + new * (acq+trans)/(nA+acq+trans)  ### previous with S&R
    if(nA>0){M_new = M_new * nA / (nA+acq+trans) + 
      new_a * acq/(nA+acq+trans) + new_b * trans/(nA+acq+trans) 
    }else{M_new =  new_a * acq/(nA+acq+trans) + new_b * trans/(nA+acq+trans)  }
    
    #*** New matrix
    bigM[,,tt] <- M_new # Grab the appropriate matrix for this timestep
    
    
    #*** Checks
    if(any(colSums(M_new)>1.001)){print(c(sum(M[,which(colSums(M_new>1))]),"colsums",colSums(M_new,"ERROR in colsums of M")));break}
    if(any(rowSums(M_new)>1.001)){print(c(sum(M[,which(rowSums(M_new>1))]),"rowsums",rowSums(M_new,"ERROR in rowsums of M")));break}
    
    #*** Calculate mean fitness
    meanfit = sum(colSums(bigM[,,tt])*vf)
    
    #*** Calculate mean susceptibility
    # This is the percentage of those at each fitness level that die * 
    #? To match previous work this should be 0.6/0.95 to convert (1-ks) [0.95] to (1-kr) [0.6]
    # 0.1 is 10% resistant: gives a number, (1- this number) is the amount of treatment success
    meanres = sum(rowSums(bigM[,,tt])*vs_rel)
  }
  
  #print(c("meanfit",meanfit))
  #*** Output mean fitness, new distributions of active and latent 
  return(list(meanfit=meanfit,vs_rel=vs_rel,bigM=bigM,meanres = meanres, vf=vf))
}


######***************************************************************************** Difference equation set up taking in above *******#####
## Function which calculates population dynamics over time
ec_srs_funcf_mean_varsr=function(endp,home,vary,initial,M0,acqdist,dt,kk){
  # Needs endp: length of simulation
  #       home: location
  #       vary: key parameters can change (epsilon, f and treatment differential (krdiff))
  #       initial: initial population (150/100K incidence etc) 
  #       M0,L0: initial distribution arrays, gives information on number of fitness levels 
  #       type: type for acquisition distribution (options: "orig" "normal")
  
  #*** Parameter assignment
  # Might not need if globally assigned
  setwd(home) # might have para in different place for different models 
  para<-read.csv("data/para_ecoli.csv",header=TRUE,check.names=F,stringsAsFactors = FALSE)[,1:2]
  for(i in 1:length(para[,1])){assign(para[i,1],para[i,2])}
  # Assign vary after to overwrite (? not checked overruled global assignment for meanfit function)
  # Parameters assigned by length of vary so for baseline just eps = 0
  vary_n = c("omega") # NB only ks important (kr not used). This is proportion of cases treated successfully if susceptible, mean resistance multiplies this
  if(length(vary)>0){for(i in 1:length(vary)){assign(vary_n[i],vary[i])}}
  
  # Correct for timestep
  mu<-mu*dt;beta<-beta*dt;eps<-eps*dt
  
  #*** Build and intialise population
  U<-matrix(0,1,endp); B<-matrix(0,1,endp);
  U[1]<-initial[1]; B[1]<-initial[2]
  
  #*** Fitness levels: Add columns to M and L if not enough time steps
  if(dim(M0)[3]<endp){M<-array(0,c(dim(M0)[1],dim(M0)[2],endp+1)); M[,,1]<-M0[,,1]; rownames(M)<-rownames(M0);}else{M<-M0}
  
  # If B present then how distributed? 
  if(B[1] > 0 ){if(sum(M0[,,1]) < 1){print("ERROR M0 distribution")}}
  
  # How many fitness levels?
  nfit<-dim(M)[2]
  # How many resistance levels?
  mres<-dim(M)[1]
  
  #*** Initial foi
  lambda<-matrix(0,1,endp)
  
  vs_mic1 <- seq(0.1,32,length.out = mres); 
  vs_rel1 <- pmax((vs_mic1 - omega),0)/vs_mic1 # constant 
  lambda[1] = sum(colSums(acqdistn)*seq(1/nfit,1,1/nfit)) * sum(rowSums(acqdistn)*vs_rel1) * beta * B[1] # function outputs just meanfit when all popns 0
  meanf<-c(0,0);
  
  print(c("mu",mu,"beta",beta,"eps",eps,"kk",kk,"omega",omega))
  
  #*** Main model dynamics
  for(i in 1:endp){
    lambda=lambda[i]
    
    #print(c("ks,kr",ks,kr,X$meanfit,lambdas))
    # Dynamics
    U[i+1] =  U[i] + mu*(B[i]) - lambda*(U[i]/(U[i] + kk))
    B[i+1] =  B[i] + lambda*(U[i]/(U[i] + kk)) - mu*B[i] 
    
    # Mean fitness update and foi
    X<-ec_srs_meanfit_varsr(M,eps*lambda*(U[i]/(U[i] + kk)),
                            (1-eps)*lambda*(U[i]/(U[i] + kk)),omega,B[i]*(1 - mu),
                            nfit,mres,acqdist)
   # print(c(i,"X",X))
    # MIC Susceptible = 1
    lambda[i+1] = X$meanfit * X$meanres * beta * B[i+1]   
    
    M<-X$bigM
    meanf<-rbind(meanf,c(X$meanfit,X$meanres))
  
  }
  
  #*** Storing and output
  All<-c();D<-c()
  All<-as.data.frame(cbind(seq(1,endp+1,1),U,B,meanf))
  colnames(All)<-c("time","U","B","mf","mr")
  # In form for plotting
  D1<-melt(All,id.vars="time")
  # What to output (storage of all vector, plus individual levels, at 5 and 50, foi s, foi r, fitness distributions and mean relative fitness over time)
  return(list(D1=D1,U=U,B=B,lambda=lambda,M=M,meanf=meanf))
}



######****************************************************************************************** For multiple plots *******#####
## Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}






######**********************************************************************************  *******#####
## Plot all output for different acqdistns
# new
plot_diff_acd_output <- function(acqdistn,plots,num) {
  # acqdistn = matrix of acqdistn
  # plots = address
  # num = identifier
  
  ### Plot necessities
  theme_set(theme_bw(base_size = 34))
  cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  ### Plot acqdistn
  z<-as.data.frame(acqdistn)
  rownames(z) <- seq(1/mres,1,1/mres);colnames(z) <- seq(1,nfit,1);
  z2<-as.data.frame(melt(z)); z2$res<-seq(1/mres,1,1/mres); colnames(z2)<-c("fitness","value","res")
  p<-ggplot(z2, aes(x=res, y=value, fill=factor(fitness))) + geom_bar(stat="identity",colour="black") + facet_grid(~fitness) 
  p<-p + scale_x_continuous("Resistance level",breaks=c(0,0.2,0.4,0.6,0.8,1)) + scale_y_continuous("Proportion") + scale_fill_brewer("Fitness \nlevel",palette="Reds") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p
  setwd(plots)
  ggsave(paste(num,"_acqdistn_05.pdf",sep=""),width=14,height=10)
  
  omega1 <- 25
  omega2 <- 2
  omega3 <- 0.4
  Sv20<-ec_funcf_mean_varsr(tsteps,home, c(omega1),iniv,M0,acqdistn,dt,500)
  Sv10<-ec_funcf_mean_varsr(tsteps,home, c(omega2),iniv,M0,acqdistn,dt,500)
  Sv05<-ec_funcf_mean_varsr(tsteps,home, c(omega3),iniv,M0,acqdistn,dt,500)
  
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
  ggsave(paste(num,"_Array_w=",omega1,"_06.pdf",sep=""))
  p9<-ggplot(mm10[w,],aes(x,y,fill=z))  + facet_wrap( ~ t, ncol=3) + ggtitle(paste("w = ", omega2,sep=""))
  p9<-p9 + scale_fill_gradient("Proportion", limits=c(0,1),low="white", high="red",guide = FALSE)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p9<-p9 + geom_tile() + scale_y_continuous(breaks=c(1,nfit),"Relative fitness levels",labels=c("Least","Most")) + scale_x_continuous(breaks=c(mres,1),"Resistance levels",labels=c("Most","Least"))
  p9
  ggsave(paste(num,"_Array_w=",omega2,"_06.pdf",sep=""))
  p4<-ggplot(mm05[w,],aes(x,y,fill=z))  + facet_wrap( ~ t, ncol=3) + ggtitle(paste("w = ", omega3,sep=""))
  p4<-p4 + scale_fill_gradient("Proportion", limits=c(0,1),low="white", high="red",guide=FALSE)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p4<-p4 + geom_tile() + scale_y_continuous(breaks=c(1,nfit),"Relative fitness levels",labels=c("Least","Most")) + scale_x_continuous(breaks=c(mres,1),"Resistance levels",labels=c("Most","Least"))
  p4
  ggsave(paste(num,"_Array_w=",omega3,"_06.pdf",sep=""))
  
  setwd(plots)
  pdf(paste(num,"_Array_w_all.pdf",sep=""),width=18,height=8)
  multiplot(p1,p9,p4,cols=3)
  dev.off()
  
  ## plot U, S & R over time
  Mu <- as.data.frame(cbind(seq(0, tsteps, 1),matrix(0,tsteps+1,1), Sv05$U, Sv10$U, Sv20$U))
  Ms <- as.data.frame(cbind(seq(0, tsteps, 1),matrix(1,tsteps+1,1), Sv05$S, Sv10$S, Sv20$S))
  Mr <- as.data.frame(cbind(seq(0, tsteps, 1),matrix(2,tsteps+1,1), Sv05$R, Sv10$R, Sv20$R))
  Mmf <- as.data.frame(cbind(seq(0, tsteps, 1),matrix(3,tsteps+1,1), Sv05$meanf[,1], Sv10$meanf[,1], Sv20$meanf[,1]))
  Mmr <- as.data.frame(cbind(seq(0, tsteps, 1),matrix(4,tsteps+1,1), Sv05$meanf[,2], Sv10$meanf[,2], Sv20$meanf[,2]))
  Mhigh <- as.data.frame(cbind(seq(0, tsteps, 1),matrix(5,tsteps+1,1), Sv05$M[5,5,], Sv10$M[5,5,],Sv20$M[5,5,]))
  Musr <- rbind(Mu, Ms, Mr,Mmf,Mmr,Mhigh)
  colnames(Musr) <- c("t", "pop",omega3, omega2, omega1)
  Msrm <- melt(Musr, id.vars = c("t","pop"))
  facet_names <- c(`0` = "U", `1` = "S", `2` = "R", `3` = "mean fit", `4` = "mean res", `5` = "Highest fit/res")
  ggplot(Msrm, aes(x=t, y = value, colour = variable)) + geom_line() + facet_wrap(~pop,labeller = as_labeller(facet_names), scales = "free")
  ggsave(paste(num,"_TimeSeries_output_06.pdf",sep=""))
  # number in highest fitness changes but mean r and f don't? 
  
  
  # plots actual numbers with resistance
  for(i in 1:length(w)){mm20[w[i],"zR"] <- Sv20$R[(1/dt)*mm20[w[i],"t"]]*mm20[w[i],"z"]}
  for(i in 1:length(w)){mm10[w[i],"zR"] <- Sv10$R[(1/dt)*mm10[w[i],"t"]]*mm10[w[i],"z"]}
  for(i in 1:length(w)){mm05[w[i],"zR"] <- Sv05$R[(1/dt)*mm05[w[i],"t"]]*mm05[w[i],"z"]}
  
  p1<-ggplot(mm20[w,],aes(x,y,fill=zR))  + facet_wrap( ~ t, ncol=3) + ggtitle(paste("w = ", omega1,sep=""))
  p1<-p1 + scale_fill_gradient("Proportion", limits=c(0,100),low="white", high="red",guide = FALSE)
  p1<-p1 + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p1<-p1 + geom_tile() + scale_y_continuous(breaks=c(1,nfit),"Relative fitness levels",labels=c("Least","Most")) + scale_x_continuous(breaks=c(mres,1),"Resistance levels",labels=c("Most","Least"))
  p1
  ggsave(paste(num,"_Array_w=",omega1,"_zR_06.pdf",sep=""))
  p9<-ggplot(mm10[w,],aes(x,y,fill=zR))  + facet_wrap( ~ t, ncol=3) + ggtitle(paste("w = ", omega2,sep=""))
  p9<-p9 + scale_fill_gradient("Proportion", limits=c(0,100),low="white", high="red",guide = FALSE)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p9<-p9 + geom_tile() + scale_y_continuous(breaks=c(1,nfit),"Relative fitness levels",labels=c("Least","Most")) + scale_x_continuous(breaks=c(mres,1),"Resistance levels",labels=c("Most","Least"))
  p9
  ggsave(paste(num,"_Array_w=",omega2,"_zR_06.pdf",sep=""))
  p4<-ggplot(mm05[w,],aes(x,y,fill=zR))  + facet_wrap( ~ t, ncol=3) + ggtitle(paste("w = ", omega3,sep=""))
  p4<-p4 + scale_fill_gradient("Proportion", limits=c(0,100),low="white", high="red",guide=FALSE)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p4<-p4 + geom_tile() + scale_y_continuous(breaks=c(1,nfit),"Relative fitness levels",labels=c("Least","Most")) + scale_x_continuous(breaks=c(mres,1),"Resistance levels",labels=c("Most","Least"))
  p4
  ggsave(paste(num,"_Array_w=",omega3,"_zR_06.pdf",sep=""))
  
  setwd(plots)
  pdf(paste(num,"_Array_w_all_zR.pdf",sep=""),width=18,height=8)
  multiplot(p4,p9,p1,cols=3)
  dev.off()
  # shows that not all potential 100 units are R
  
  ### Look at proportion in each of the 30 levels over time for each - facet = level
  mm20$omega = omega1; mm10$omega = omega2; mm05$omega = omega3
  mega<-as.data.frame(rbind(mm20,mm10,mm05)); colnames(mega)<-c("x","y","z","time","omega")
  mega$level = c(seq(21,25,1),seq(16,20,1),seq(11,15,1),seq(6,10,1),seq(1,5,1))
  g<-ggplot(mega,aes(x=time,y=z,colour=factor(omega))) + geom_line(size=2) + facet_wrap( ~ level, ncol=5) +  scale_colour_manual(values=cbPalette,"Abx\nLevel",breaks=c(omega1,omega2, omega3))
  g<-g + scale_x_continuous("Generations",breaks=c(0,200,400)) + scale_y_continuous("Proportion at this level",breaks=c(0.25,0.75))+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
  g
  setwd(plots)
  ggsave(paste(num,"_mega.pdf",sep=""),width=12,height=8)
  
  theme_set(theme_bw(base_size = 14))
  g<-g+ facet_wrap( ~ level,scales = "free", ncol=5)+ scale_y_continuous("Proportion at this level")
  ggsave(paste(num,"_mega_freescale.pdf",sep=""))
  
  ### Look at change in R & S over time
  theme_set(theme_bw(base_size = 34))
  rrmr<-as.data.frame(cbind(seq(1,ll,1)*dt,Sv05$R,Sv10$R,Sv20$R,1)); colnames(rrmr)<-c("time","05","10","20","type")
  rrms<-as.data.frame(cbind(seq(1,ll,1)*dt,Sv05$S,Sv10$S,Sv20$S,2)); colnames(rrms)<-c("time","05","10","20","type")
  rrm<-as.data.frame(rbind(rrmr,rrms))
  rrm2<-melt(rrm,id.vars=c("time","type")); rrm2[which(rrm2$type==1),"type"]<-"Resistant"; rrm2[which(rrm2$type==2),"type"]<-"Susceptible"
  g<-ggplot(rrm2,aes(x=time,y=value,colour=factor(variable))) + geom_line(size=2) + scale_x_continuous("Generations") + scale_y_continuous("Percentage with R")
  g<-g + scale_colour_manual("Abx\nLevel",breaks=c("20","10","05"),labels=c(omega1, omega2,omega3),values = cbPalette) + facet_wrap(~type)
  g
  ggsave(paste(num,"_r&s_overtime.pdf",sep=""),width=12,height=8)
  
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
  ggsave(paste(num,"_mega_normtodom.pdf",sep=""))
  
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
  ggsave(paste(num,"_mega_normtofullR.pdf",sep=""),width=12,height=8)
  
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
  ggsave(paste(num,"_f&r_overtime.pdf",sep=""),width=18,height=12)
  
  #### Compare with and without fitness and resistance levels. 
  ### Range of omega
  setwd(home)
  omegav <- c(omega1,omega2,omega3)
  para<-read.csv("data/para_ecoli.csv",header=TRUE,check.names=F,stringsAsFactors = FALSE)[,1:2]
  for(i in 1:length(para[,1])){assign(para[i,1],para[i,2])}
  # Correct for timestep
  mu<-mu*dt;beta<-beta*dt;eps<-eps*dt
  m<-dim(acqdistn)[1]; vs<-seq(1/m,1,1/m); 
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
      # NEW Dynamics
      U[i+1] =  U[i] + mu*(S[i]+R[i]) - (lambdas+lambdar)*(U[i]/(U[i] + kk)) 
      S[i+1] =  S[i] + lambdas*(U[i]/(U[i] + kk)) - mu*S[i] - eps * S[i]
      R[i+1] =  R[i] + lambdar*(U[i]/(U[i] + kk)) - mu*R[i] + eps * S[i] 
      
      lambdasv[i+1] =     max(0,(1-omega)/1) * beta * S[i+1] / ( S[i+1] + R[i+1] );
      lambdarv[i+1] = f * max(0,(20-omega)/20) * beta * R[i+1] / ( S[i+1] + R[i+1] );   # resistant strain has an MIC of 6  
    } 
    all<-as.data.frame(cbind(seq(0,endp,1)*dt,U,S,R,omega)); colnames(all)<-c("time","U","Susceptible","Resistant","w")
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
  
  w<-which(allmn$with == 0)
  p<-ggplot(allmn[w,],aes(x=time,y=value,colour=variable,linetype=factor(with)))+geom_line(size=2) + 
    scale_x_continuous("Time steps",lim=c(0,endp*dt))
  p<-p+scale_colour_manual("Sub-\npopulation",breaks=c("Susceptible","Resistant"), values = c("blue","red")) + 
    scale_y_continuous("Percentage of population", limits = c(0,100)) + facet_wrap( ~ nw)  + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p<-p + scale_linetype_discrete("With\ndiversity",breaks=c(0,1),labels=c("None","With diversity")) + theme(legend.position="none")
  p
  setwd(plots)
  ggsave(paste(num,"_Withoutdiversity.pdf",sep=""),width=12,height=7)
  
  p<-ggplot(allmn,aes(x=time,y=value,colour=variable,linetype=factor(with)))+geom_line(size=2) + 
    scale_x_continuous("Time steps",lim=c(0,endp*dt))
  p<-p+scale_colour_manual("Sub-\npopulation",breaks=c("Susceptible","Resistant"), values = c("blue","red")) + 
    scale_y_continuous("Percentage of population", limits = c(0,100)) + facet_wrap( ~ nw)  + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p<-p + scale_linetype_discrete("With\ndiversity",breaks=c(0,1),labels=c("None","With diversity")) 
  p
  setwd(plots)
  ggsave(paste(num,"_WithnWithoutdiversity.pdf",sep=""),width=12,height=7)
  p + theme(legend.position="none")
  ggsave(paste(num,"_WithnWithoutdiversity_nolegend.pdf",sep=""),width=12,height=7)
  
  p + scale_y_continuous("Percentage of population",lim=c(0,10)) 
  ggsave(paste(num,"_WithnWithoutdiversity_zoom.pdf",sep=""),width=12,height=7)
  
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
  ggsave(paste(num,"_Gulberg5d.pdf",sep=""),width=16,height=8)
  
}


######**********************************************************************************  *******#####
## Plot all output for different acqdistns
# old

plot_diff_acd_output <- function(acqdistn,plots,num) {
  # acqdistn = matrix of acqdistn
  # plots = address
  # num = identifier
  
  ### Plot necessities
  theme_set(theme_bw(base_size = 34))
  cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  ### Plot acqdistn
  z<-as.data.frame(acqdistn)
  rownames(z) <- seq(1/mres,1,1/mres);colnames(z) <- seq(1,nfit,1);
  z2<-as.data.frame(melt(z)); z2$res<-seq(1/mres,1,1/mres); colnames(z2)<-c("fitness","value","res")
  p<-ggplot(z2, aes(x=res, y=value, fill=factor(fitness))) + geom_bar(stat="identity",colour="black") + facet_grid(~fitness) 
  p<-p + scale_x_continuous("Resistance level",breaks=c(0,0.2,0.4,0.6,0.8,1)) + scale_y_continuous("Proportion") + scale_fill_brewer("Fitness \nlevel",palette="Reds") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p
  setwd(plots)
  ggsave(paste(num,"_acqdistn_05.pdf",sep=""),width=14,height=10)
  
  # Initial conditions
  Sv20<-ec_funcf_mean_varsr(tsteps,home, c(0.2),iniv,M0,acqdistn,dt)
  Sv10<-ec_funcf_mean_varsr(tsteps,home, c(0.1),iniv,M0,acqdistn,dt)
  Sv05<-ec_funcf_mean_varsr(tsteps,home, c(0.05),iniv,M0,acqdistn,dt)
  
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
  setwd(plots)
  # Grab a subset
  sub<-c(1/dt,100,250,500,750,1000,1500,2000,2500,seq(3000,4001,500))*dt
  w<-which(mm20[,"t"] %in% sub)
  # plots
  p1<-ggplot(mm20[w,],aes(x,y,fill=z))  + facet_wrap( ~ t, ncol=3) + ggtitle("w = 0.2")
  p1<-p1 + scale_fill_gradient("Proportion", limits=c(0,1),low="white", high="red",guide = FALSE)
  p1<-p1 + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p1<-p1 + geom_tile() + scale_y_continuous(breaks=c(1,nfit),"Relative fitness levels",labels=c("Least","Most")) + scale_x_continuous(breaks=c(mres,1),"Resistance levels",labels=c("Most","Least"))
  p1
  ggsave(paste(num,"_Array_w=20.pdf",sep=""))
  p9<-ggplot(mm10[w,],aes(x,y,fill=z))  + facet_wrap( ~ t, ncol=3) + ggtitle("w = 0.1")
  p9<-p9 + scale_fill_gradient("Proportion", limits=c(0,1),low="white", high="red",guide = FALSE)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p9<-p9 + geom_tile() + scale_y_continuous(breaks=c(1,nfit),"Relative fitness levels",labels=c("Least","Most")) + scale_x_continuous(breaks=c(mres,1),"Resistance levels",labels=c("Most","Least"))
  p9
  ggsave(paste(num,"_Array_w=10.pdf",sep=""))
  p4<-ggplot(mm05[w,],aes(x,y,fill=z))  + facet_wrap( ~ t, ncol=3) + ggtitle("w = 0.05")
  p4<-p4 + scale_fill_gradient("Proportion", limits=c(0,1),low="white", high="red",guide=FALSE)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p4<-p4 + geom_tile() + scale_y_continuous(breaks=c(1,nfit),"Relative fitness levels",labels=c("Least","Most")) + scale_x_continuous(breaks=c(mres,1),"Resistance levels",labels=c("Most","Least"))
  p4
  ggsave(paste(num,"_Array_w=05.pdf",sep=""))
  
  setwd(plots)
  pdf(paste(num,"_Array_w_all.pdf",sep=""),width=20,height=18)
  multiplot(p1,p4,cols=2)
  dev.off()
  
  ### Look at proportion in each of the 30 levels over time for each - facet = level
  mm20$omega = 20; mm10$omega = 10; mm05$omega = 5
  mega<-as.data.frame(rbind(mm20,mm10,mm05)); colnames(mega)<-c("x","y","z","time","omega")
  mega$level = c(seq(21,25,1),seq(16,20,1),seq(11,15,1),seq(6,10,1),seq(1,5,1))
  g<-ggplot(mega,aes(x=time,y=z,colour=factor(omega))) + geom_line(size=2) + facet_wrap( ~ level, ncol=5) + scale_colour_manual("Abx\nLevel",breaks=c(20,10,5),labels=c("0.2","0.1","0.05"),values = cbPalette)
  g<-g + scale_x_continuous("Generations",breaks=c(0,200,400)) + scale_y_continuous("Proportion at this level",breaks=c(0.25,0.75))+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
  g
  setwd(plots)
  ggsave(paste(num,"_mega.pdf",sep=""))
  
  theme_set(theme_bw(base_size = 14))
  g<-g+ facet_wrap( ~ level,scales = "free", ncol=5)+ scale_y_continuous("Proportion at this level")
  ggsave(paste(num,"_mega_freescale.pdf",sep=""))
  
  ### Look at change in R & S over time
  rrmr<-as.data.frame(cbind(seq(1,ll,1)*dt,Sv05$R,Sv10$R,Sv20$R,1)); colnames(rrmr)<-c("time","05","10","20","type")
  rrms<-as.data.frame(cbind(seq(1,ll,1)*dt,Sv05$S,Sv10$S,Sv20$S,2)); colnames(rrms)<-c("time","05","10","20","type")
  rrm<-as.data.frame(rbind(rrmr,rrms))
  rrm2<-melt(rrm,id.vars=c("time","type")); rrm2[which(rrm2$type==1),"type"]<-"Resistant"; rrm2[which(rrm2$type==2),"type"]<-"Susceptible"
  g<-ggplot(rrm2,aes(x=time,y=value,colour=factor(variable))) + geom_line(size=2) + scale_x_continuous("Generations") + scale_y_continuous("Percentage with R")
  g<-g + scale_colour_manual("Abx\nLevel",breaks=c("20","10","05"),labels=c("0.2","0.1","0.05"),values = cbPalette) + facet_wrap(~type)
  g
  ggsave(paste(num,"_r&s_overtime.pdf",sep=""),width=12,height=8)
  
  # Time to dominance...
  t05<-min(intersect(intersect(which(rrmr[,2]>79.99),which(rrmr[,2]< 80.9)),which(floor(rrmr[,2])==80))*dt)
  t10<-min(intersect(intersect(which(rrmr[,3]>79.99),which(rrmr[,3]< 80.001)),which(floor(rrmr[,3])==80))*dt)
  t20<-min(intersect(intersect(which(rrmr[,4]>79.99),which(rrmr[,4]< 80.007)),which(floor(rrmr[,4])==80))*dt)
  
  mm20_2<-mm20;mm10_2<-mm10;mm05_2<-mm05
  mm20_2$t<-mm20_2$t/t20;mm10_2$t<-mm10_2$t/t10;mm05_2$t<-mm05_2$t/t05
  mega_2<-as.data.frame(rbind(mm20_2,mm10_2,mm05_2)); colnames(mega_2)<-c("x","y","z","time","omega")
  mega_2$level = c(seq(21,25,1),seq(16,20,1),seq(11,15,1),seq(6,10,1),seq(1,5,1))
  g<-ggplot(mega_2,aes(x=time,y=z,colour=factor(omega))) + geom_line(size=2) + facet_wrap( ~ level, ncol=5, scales = "free") + scale_colour_manual("Abx\nLevel",breaks=c(20,10,5),labels=c("0.2","0.1","0.05"),values = cbPalette)
  g<-g + scale_x_continuous("Time to full resistance") + scale_y_continuous("Proportion at this level")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
  g
  setwd(plots)
  ggsave(paste(num,"_mega_normtodom.pdf",sep=""))
  
  # Time to full resistance and fitness
  w<-which(mega$level==5); 
  t05a<-min(intersect(which(mega[w,"omega"]==5),which(mega[w,"z"]>0.75)))
  t10a<-min(intersect(which(mega[w,"omega"]==10),which(mega[w,"z"]>0.75)))
  t20a<-min(intersect(which(mega[w,"omega"]==20),which(mega[w,"z"]>0.75)))
  
  t05<-mega[w[t05a],"time"]; t10<-mega[w[t10a],"time"]; t20<-mega[w[t20a],"time"]
  
  mm20_2<-mm20;mm10_2<-mm10;mm05_2<-mm05
  mm20_2$t<-mm20_2$t/t20;mm10_2$t<-mm10_2$t/t10;mm05_2$t<-mm05_2$t/t05
  mega_2<-as.data.frame(rbind(mm20_2,mm10_2,mm05_2)); colnames(mega_2)<-c("x","y","z","time","omega")
  mega_2$level = c(seq(21,25,1),seq(16,20,1),seq(11,15,1),seq(6,10,1),seq(1,5,1))
  g<-ggplot(mega_2,aes(x=time,y=z,colour=factor(omega))) + geom_line(size=2) + facet_wrap( ~ level, ncol=5, scales = "free") + scale_colour_manual("Abx\nLevel",breaks=c(20,10,5),labels=c("0.2","0.1","0.05"),values = cbPalette)
  g<-g + scale_x_continuous("Time to full resistance") + scale_y_continuous("Proportion at this level")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
  g
  setwd(plots)
  ggsave(paste(num,"_mega_normtofullR.pdf",sep=""))
  
  #### Compare with and without fitness and resistance levels. 
  ### Range of omega
  setwd(home)
  print(getwd())
  
  omegav<-c(0.2,0.1,0.05) * dt
  para<-read.csv("data/para_ecoli.csv",header=TRUE,check.names=F,stringsAsFactors = FALSE)[,1:2]
  for(i in 1:length(para[,1])){assign(para[i,1],para[i,2])}
  # Correct for timestep
  mu<-mu*dt;beta<-beta*dt;eps<-eps*dt
  assign("ks",1) # Prescribe effect treatment success levels
  assign("kr",1-sum(rowSums(acqdistn)*vs)) # Set the levels to be those of the acqdistn
  assign("f",sum(colSums(acqdistn)*vf))
  bigall<-c(); lambdasv<-c(); lambdarv<-c(); 
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
  colnames(rrm2)<-c("time","variable","w","value")
  rrm2n<-rrm2[,c("w","time","variable","value")]; 
  rrm2n$with<-1; rrm2n$nw<-0
  allm$with<-0; allm$nw<-allm$w; 
  rrm2n[which(rrm2n$w=="05"),"nw"]=allm[12000,"nw"]; 
  rrm2n[which(rrm2n$w=="10"),"nw"]<-allm[4000,"nw"]; 
  rrm2n[which(rrm2n$w=="20"),"nw"]<-allm[1,"nw"]
  allmn<-rbind(allm,rrm2n)
  theme_set(theme_bw(base_size = 34))
  p<-ggplot(allmn,aes(x=time,y=value,colour=variable,linetype=factor(with)))+geom_line(size=2) + scale_x_continuous("Generations",lim=c(0,endp*dt))
  p<-p+scale_colour_discrete("Sub-\npopulation",breaks=c("Susceptible","Resistant")) + scale_y_continuous("Percentage of population") + facet_wrap( ~ nw)  + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p<-p+ scale_linetype_discrete("With\ndiversity",breaks=c(0,1),labels=c("None","With diversity")) 
  p
  setwd(plots)
  ggsave(paste(num,"_WithnWithoutdiversity.pdf",sep=""),width=16,height=8)
  
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
  theme_set(theme_bw(base_size = 34))
  pp<-as.data.frame(pp);colnames(pp)<-c("t","Fitness level 1\n(low)","Fitness level 2","Fitness level 3","Fitness level 4","Fitness level 5\n(high)","Res. level 1\n(low)","Res. level 2","Res. level 3","Res. level 4","Res. level 5\n(high)","w"); 
  pp2<-melt(pp,id.vars = c("t","w"))
  g<-ggplot(pp2,aes(x=t,y=value,colour=factor(w))) + facet_wrap(~variable,ncol=5) + geom_line(size=2)
  g<-g + scale_x_continuous("Generation") + scale_y_continuous("Proportion") + scale_colour_manual("Abx\nlevel",labels=c(0.05,0.1,0.2),values=cbPalette)
  g # slightly slower to the most fit levels but v similar. Bigger delay for low level resistance
  setwd(plots)
  ggsave(paste(num,"_f&r_overtime.pdf",sep=""),width=18,height=12)
}
