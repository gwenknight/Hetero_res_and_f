##### Code for generalised fitness function 
# With the inclusion of different resistance levels
# Resistance = increased ability to survive. 

######*********************************************************************** Main generalised function code  *******#####
## Function to update the mean fitness & MEAN SUSCEPTIBILITY at each timestep
stoch_ec_srs_meanfit_varsr=function(bigM,acq,trans,trate,nA,n,m,acqdist,submic){
  # Needs bigM: Array of distribution  [susceptible rows, fitness columns, time]
  #       acq: number of acquisitions of new strains
  #       trans: number of new bacterial strains
  #       trate: antibiotic exposure level
  #       nA: number with bacteria
  #       n: number of fitness levels
  #       m: number of resistance levels
  #       acqdist: 2D distribution of new cases
  #       submic: pattern below MIC
  
  # Vector of fitness values
  vf<-seq(1/n,1,1/n)
  
  # VECTOR of MIC LEVELs => ?Susceptible has MIC of 1.  
  # vector of MIC values
  # MAX MIC = 32
  # MIN MIC = 0.1
  vs_mic <- seq(0.1,32,length.out = m) 
  
  ## How does growth rate change sub-mic? 
  if(submic == 1){
    # Vector of resistance levels in face of antibiotic treatment trate
    vs_rel <- pmax((vs_mic - trate),0)/vs_mic} # constant 
  if(submic == 2){
    # # FLAT AND THEN DROPS LINEARLY AT 9/10 of MIC to 0
    vs_rel_multiply <- dirac_att(vs_mic - trate) # new function (see below): 0 if vs_mic < trate, 1 if greater
    for(i in 1:length(vs_mic)){ if((0.9*vs_mic[i]) > trate){ # if below 9/10ths of MIC then relative growth flat
      vs_rel[i] = 1} else {
        vs_rel[i] = 10*(vs_mic[i] - trate) / vs_mic[i]} } # else linear decline to 0 at MIC
    vs_rel <- vs_rel_multiply * vs_rel}
  if(submic == 3){
    # FLAT AND THEN DROPS LINEARLY AT 9/10 of MIC to 0.5
    vs_rel_multiply <- dirac_att(vs_mic - trate) # new function (see below): 0 if vs_mic < trate, 1 if greater
    for(i in 1:length(vs_mic)){ if((0.9*vs_mic[i]) > trate){ # if below 9/10ths of MIC then relative growth flat
      vs_rel[i] = 1} else {
        vs_rel[i] = (11*vs_mic[i] - 10*trate) / (2*vs_mic[i]) } } # else linear decline to 0 at MIC
    vs_rel <- vs_rel_multiply * vs_rel
    vs_rel <- dirac_att_5(vs_rel) # Make 0.5 instead of 0 above MIC
  }
  
  # Which column of M? last that is non-zero plus one for this timestep
  sumM<-colSums(colSums(bigM,dims=1)) # Gives a vector of the sums over the array for each timestep
  if(length(which(sumM > 0))){tt<-tail(which(sumM > 0), n = 1) + 1}else{tt=2} # which timestep.
  
  M_new <- bigM[,,tt-1] # Grab the appropriate matrix for this timestep, use as ew matrix to update as go through
  
  # If there are acq and transmissions then M_new needs to be updated
  if(acq != 0 | trans !=0){
    #**** Update 
    ##### 1) New cases - how distribute acquisitions across fitness and resistance? 
    
    #*** How are the acquisitions distributed? Input STOCHASTIC
    # sample using acqdist
    new_a <- acqdist
    
    #*** How are the new growth bugs distributed?
    M_temp<-M_new
    # if none resistant then none grow... 
    res_levels = sum( vs_rel )
    if(trans > 0){ # if no growth then no need to do
      if(res_levels > 0){ # if resistant then do the following... otherwise no new bugs
        
        # now growth rate dependent on fitness and resistance 
        # update M_temp with fitness then resistance
        
        # FITNESS
        Mf<-colSums(M_temp) # Proportions at each FITNESS level NOT STOCHASTIC - NUMBER STOCHASTIC BUT NOT HOW FITNESS AFFECTS
        pastmean=sum( Mf*vf )
        M_temp <- t(t(M_temp) * vf / pastmean)
        
        # RESISTANCE
        Mr<-rowSums(M_temp) # = rowSums(M_new) same # Proportions at each RESISTANCE level
        pastmean=sum( Mr*vs_rel )
        if(pastmean > 0){ M_temp <- vs_rel*M_temp/pastmean } # Otherwise leave M_temp as is. 
        
        # new transmission are then distributed
        new_b = M_temp 
        
      }else{new_b = 0; trans = 0} # if all resistant then no growth
    }else{new_b = 0} # no distribution of new cases
    
    #*** Assign and update M
    # nA = number of bacteria left after treatment removed (affects distribution) and death (doesn't affect)
    #M_new = M_new * nA / (nA+acq+trans) + new_a * acq / (nA+acq+trans) + new_b * trans / (nA+acq+trans) 
    M_new <- M_new * nA + rMydist_a(acq, new_a) + rMydist_a(trans, new_b) # use numbers of events to select
    M_new <- M_new / (nA+acq+trans) 
  }
  
  #*** New matrix or same as before
  bigM[,,tt] <- M_new # Grab the appropriate matrix for this timestep
  
  #*** Checks
  if(any(colSums(M_new)>1.001)){print(c(sum(M_new[,which(colSums(M_new)>1.001)]),"colsums",colSums(M_new),"ERROR in colsums of M (too high)"));break}
  if(any(rowSums(M_new)>1.001)){print(c(sum(M_new[,which(rowSums(M_new)>1.001)]),"rowsums",rowSums(M_new),"ERROR in rowsums of M (too high)"));break}
  if(sum(M_new)<0.999){print(c(sum(M_new),"colsums",colSums(M_new),"rowsums",rowSums(M_new),"ERROR in sum of M (too small)"));break}
  
  #*** Calculate mean fitness
  meanfit = sum(colSums(bigM[,,tt])*vf)
  #print(c("meanres = sum of", rowSums(bigM[,,tt]),"vs_rel",vs_rel,"sum",sum(bigM[,,tt])))
  #print(c("cc",colSums(bigM[,,tt]),"rr",rowSums(bigM[,,tt])))
  #*** Calculate mean susceptibility
  meanres = sum(rowSums(bigM[,,tt])*vs_rel)
  
  
  #*** Output mean fitness, new distributions of active and latent 
  return(list(meanfit=meanfit, vs_rel=vs_rel, bigM=bigM, meanres = meanres, vf=vf))
}

######***************************************************************************** Difference equation set up taking in above *******#####
## Function which calculates population dynamics over time
stoch_ec_srs_funcf_mean_varsr=function(endp,home,vary,initial,M0,acqdist,dt,submic){
  # Needs endp: length of simulation
  #       home: location
  #       vary: key parameters can change (OMEGA [abx], kk)
  #       initial: initial population
  #       M0: initial distribution arrays, gives information on number of fitness levels 
  #       acqdist: incoming distribution of new mutations
  #       dt = timestep
  #       submic: pattern below MIC
  
  #*** Parameter assignment
  # Might not need if globally assigned
  setwd(home) # might have para in different place for different models 
  para<-read.csv("data/para_ecoli.csv",header=TRUE,check.names=F,stringsAsFactors = FALSE)[,1:2]
  for(i in 1:length(para[,1])){assign(para[i,1],para[i,2])}
  vary_n = c("kk","omega") 
  
  if(length(vary) > 0){
    if(length(vary) < 3){for(i in 1:length(vary)){assign(vary_n[i],vary[i])}
    }else{ assign("kk",vary[1]); assign("omega",vary[2:length(vary)])}
  }
  # make a matrix of omega values
  if(length(omega)==1){omega <- matrix(omega, 1, endp)}
  
  # Correct for timestep
  mu<-mu*dt;  beta<-beta*dt;
  
  #*** Build and intialise population
  U<-matrix(0,1,endp);      B<-matrix(0,1,endp);
  U[1]<-initial[1];         B[1]<-initial[2]
  
  #*** Fitness levels: Add columns to M and L if not enough time steps
  if(dim(M0)[3]<endp){M<-array(0,c(dim(M0)[1],dim(M0)[2],endp+1)); M[,,1]<-M0[,,1]; rownames(M)<-rownames(M0);}else{M<-M0}
  
  # If B present then how distributed? 
  if(B[1] > 0 ){if(sum(M0[,,1]) < 1){print("ERROR M0 distribution")}}
  
  # How many fitness levels?
  nfit<-dim(M)[2]
  # How many resistance levels?
  mres<-dim(M)[1]
  
  #*** Initial foi
  lambda_v<-matrix(0,1,endp)
  
  vs_mic1 <- seq(0.1,32,length.out = mres); 
  vs_rel1 <- pmax((vs_mic1 - omega[1]),0)/vs_mic1 # constant 
  lambda_v[1] = sum(colSums(acqdistn)*seq(1/nfit,1,1/nfit)) * sum(rowSums(acqdistn)*vs_rel1) * beta  # function outputs just meanfit when all popns 0
  
  # Store place
  meanf<-c(0,0);
  print(c("mu",mu,"beta",beta,"eps",eps,"kk",kk,"mean_omega",mean(omega)))
  
  #*** Main model dynamics
  for(i in 1:endp){
    #print(i)
    
    # Deaths - uniform sample. If less than mu then die (Monte Carlo?) 
    # mu is the probability in this time step that the bacteria dies
    mu_n <- runif(B[i],min = 0, max = 1)
    die_b <- length(which(mu_n < mu))
    #if(die_b > 0){print(c("death",i,die_b))}
    U[i+1] = U[i] + die_b
    B[i+1] = B[i] - die_b
    
    # Births - uniform sample. 
    # birth is the probability in this time step that the bacteria multiplies
    bi_n <- runif(B[i+1],min = 0, max = 1) # already killed some - those that were going to die wouldn't have multiplied
    # also need to scale by U resource left:
    bi_b <- length(which(bi_n < (U[i+1]/(U[i+1] + kk)) * lambda_v[i])) # X$meanfit * X$meanres * beta => past fitness
    
    #if(bi_b > 0){print(c("birth",i,bi_b))}
    # Dynamics 
    if(  U[i] < bi_b ){ # can't have negative 
      U[i+1] = 0
      B[i+1] = B[i] + U[i] 
    } else {
      U[i+1] =  U[i] - bi_b
      B[i+1] =  B[i] + bi_b }
    
    # mutations
    mut_a <- runif(bi_b, min = 0, max = 1)
    mut_e <- length(which(mut_a < mu))
    #if(mut_e > 0){print(c("mut",i,mut_e))}
    
    #print(c("num with mut", mut_e,"bi_b",bi_b, "B[i]", B[i],"lambda",lambda_v[i]))
    # Mean fitness update and foi
    
    X<-stoch_ec_srs_meanfit_varsr(M,mut_e,bi_b - mut_e,
                                  omega[i],B[i] - die_b,
                                  nfit,mres,acqdist,submic)
    
    # Update fitness
    if(typeof(X) != "double"){
      #print(c(X$meanfit , X$meanres , beta))
      lambda_v[i+1] = X$meanfit * X$meanres * beta  # no B[i] here anymore
      
      # Store
      M<-X$bigM
      meanf<-rbind(meanf,c(X$meanfit,X$meanres))
    }else{meanf<-rbind(meanf,c(0,0))}
    
    if(B[i+1] == 0){ break }
    
  }
  
  
  
  #*** Storing and output
  All<-c();D<-c()
  All<-as.data.frame(cbind(seq(1,endp+1,1),U,B,meanf))
  colnames(All)<-c("time","U","B","mf","mr")
  # In form for plotting
  D1<-melt(All,id.vars="time")
  # What to output (storage of all vector, plus individual levels, at 5 and 50, foi s, foi r, fitness distributions and mean relative fitness over time)
  return(list(D1=D1,U=U,B=B,lambda_v=lambda_v,M=M,meanf=meanf))
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
plot_diff_acd_output_stoch <- function(acqdistn,plots,num, omega_M, submic_M, wildtype,pref) {
  # acqdistn = matrix of acqdistn
  # plots = address
  # num = identifier of acquisition distribution 
  # omega_M = omega values (max 4 in vector or 4 rows in matrix)
  # submic_M = submic indicator functions (max 4 in vector)
  # c(susr,susf) = Wildtype = level to remove in renormalising
  # label = for files (may not be needed on top of giving a plot home directory)
  
  ### Plot acqdistn
  z<-as.data.frame(acqdistn)
  mres <- dim (acqdistn)[1]; nfit <- dim (acqdistn)[2]
  rownames(z) <- seq(1/mres,1,1/mres);colnames(z) <- seq(1,nfit,1);
  z2<-as.data.frame(melt(z)); z2$res<-seq(1/mres,1,1/mres); colnames(z2)<-c("fitness","value","res")
  p<-ggplot(z2, aes(x=res, y=value, fill=factor(fitness))) + geom_bar(stat="identity",colour="black") + facet_grid(~fitness) 
  p<-p + scale_x_continuous("Resistance level",breaks=c(0,0.2,0.4,0.6,0.8,1)) + scale_y_continuous("Proportion") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) #+scale_fill_brewer("Fitness \nlevel",palette="Reds")
  p
  setwd(plots)
  ggsave(paste(num,"_",pref,"_acqdistn_",num,".pdf",sep=""),width=14,height=10)
  
  ### What to run
  omegam <- as.data.frame(matrix(0,tsteps,5))
  omegam[,1] <- seq(1,tsteps); colnames(omegam) <- c("time","1","2","3","4"); 
  if(is.vector(omega_M)){
    for(iij in 1:length(omega_M)){
      omega_value <- omega_M[iij]
      submic_value <- submic_M[iij]
      assign(paste("omega",iij,sep=""),omega_M[iij])
      assign(paste("Sv",iij,sep=""),stoch_ec_srs_funcf_mean_varsr(tsteps,home, c(kk,omega_value),iniv,M0,acqdistn,dt,submic_value))
      omegam[,(iij+1)]<-omega_value # store
    }
    
  } else {
    for(iij in 1:dim(omega_M)[1]){ # each row a different omega
      omega_value <- omega_M[iij,]
      submic_value <- submic_M[iij]
      assign(paste("omega",iij,sep=""),omega_M[iij,1])
      assign(paste("Sv",iij,sep=""),stoch_ec_srs_funcf_mean_varsr(tsteps,home, c(kk,omega_value),iniv,M0,acqdistn,dt,submic_value))
      omegam[,(iij+1)]<-omega_value
    }
  }
  
  ### What does omega do? 
  momegam <- melt(omegam,id.vars = "time")
  go<-ggplot(momegam,aes(x=time,y=value,colour=variable)) + geom_line(size=2) + scale_x_continuous(breaks = seq(0,tsteps,tsteps / 10), labels = dt*seq(0,tsteps,tsteps/10)) + scale_y_continuous("Omega") + scale_color_discrete("")
  setwd(plots)
  ggsave(paste(num,"_",pref,"norsd_omegaconcs.pdf",sep=""))
  
  # What happens? 
  mm1<-c() ; mm2<-c() ; mm3<-c() ; mm4<-c() 
  ll<-dim(Sv1$M)[3];
  ss<-seq(0,ll,1/dt) # Don't want to grab all 
  for(i in 2:length(ss)){
    mm01<-as.data.frame(melt(Sv1$M[,,ss[i]])); 
    mm02<-as.data.frame(melt(Sv2$M[,,ss[i]])); 
    mm03<-as.data.frame(melt(Sv3$M[,,ss[i]])); 
    mm04<-as.data.frame(melt(Sv4$M[,,ss[i]])); 
    mm01$tstep=ss[i]*dt; mm02$tstep=ss[i]*dt; mm03$tstep=ss[i]*dt; mm04$tstep=ss[i]*dt # To go to generations
    mm1<-rbind(mm1,mm01);mm2<-rbind(mm2,mm02); mm3<-rbind(mm3,mm03) ; mm4<-rbind(mm4,mm04) 
  } 
  colnames(mm1)<-c("x","y","z","t"); colnames(mm2)<-c("x","y","z","t"); colnames(mm3)<-c("x","y","z","t") ;colnames(mm4)<-c("x","y","z","t")
  setwd(plots)
  # Grab a subset
  sub<-c(1/dt,100,250,500,750,1000,1500,2000,2500,seq(3000,4001,500))*dt
  w<-which(mm1[,"t"] %in% sub)
  # plots
  p1<-ggplot(mm1[w,],aes(x,y,fill=z))  + facet_wrap( ~ t, ncol=3) + ggtitle(paste("w = ", 1,sep=""))
  p1<-p1 + scale_fill_gradient("Proportion", limits=c(0,1),low="white", high="red",guide = FALSE)
  p1<-p1 + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p1<-p1 + geom_tile() + scale_y_continuous(breaks=c(1,nfit),"Relative fitness levels",labels=c("Least","Most")) + scale_x_continuous(breaks=c(mres,1),"Resistance levels",labels=c("Most","Least"))
  p1
  ggsave(paste(num,"_",pref,"norsd_Array_w=", omega1,"_06.pdf",sep=""))
  p2<-ggplot(mm2[w,],aes(x,y,fill=z))  + facet_wrap( ~ t, ncol=3) + ggtitle(paste("w = ", 2,sep=""))
  p2<-p2 + scale_fill_gradient("Proportion", limits=c(0,1),low="white", high="red",guide=FALSE)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p2<-p2 + geom_tile() + scale_y_continuous(breaks=c(1,nfit),"Relative fitness levels",labels=c("Least","Most")) + scale_x_continuous(breaks=c(mres,1),"Resistance levels",labels=c("Most","Least"))
  p2
  ggsave(paste(num,"_",pref,"norsd_Array_w=",omega2,"_06.pdf",sep=""))
  p3<-ggplot(mm3[w,],aes(x,y,fill=z))  + facet_wrap( ~ t, ncol=3) + ggtitle(paste("w = ", 3,sep=""))
  p3<-p3 + scale_fill_gradient("Proportion", limits=c(0,1),low="white", high="red",guide = FALSE)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p3<-p3 + geom_tile() + scale_y_continuous(breaks=c(1,nfit),"Relative fitness levels",labels=c("Least","Most")) + scale_x_continuous(breaks=c(mres,1),"Resistance levels",labels=c("Most","Least"))
  p3
  ggsave(paste(num,"_",pref,"norsd_Array_w=",omega3,"_06.pdf",sep=""))
  p4<-ggplot(mm4[w,],aes(x,y,fill=z))  + facet_wrap( ~ t, ncol=3) + ggtitle(paste("w = ", 4,sep=""))
  p4<-p4 + scale_fill_gradient("Proportion", limits=c(0,1),low="white", high="red",guide=FALSE)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p4<-p4 + geom_tile() + scale_y_continuous(breaks=c(1,nfit),"Relative fitness levels",labels=c("Least","Most")) + scale_x_continuous(breaks=c(mres,1),"Resistance levels",labels=c("Most","Least"))
  p4
  ggsave(paste(num,"_",pref,"norsd_Array_w=",omega4,"_06.pdf",sep=""))
  
  
  setwd(plots)
  pdf(paste(num,"_",pref,"norsd_Array_w_all.pdf",sep=""),width=18,height=18)
  multiplot(p1,p2,p4,p3,cols=2)
  dev.off()
  
  ## plot U, S & R over time
  Mu <- as.data.frame(cbind(seq(0, tsteps, 1),matrix(0,tsteps+1,1), Sv1$U, Sv2$U,Sv3$U, Sv4$U))
  Mb <- as.data.frame(cbind(seq(0, tsteps, 1),matrix(1,tsteps+1,1), Sv1$B, Sv2$B,Sv3$B, Sv4$B))
  Mmf <- as.data.frame(cbind(seq(0, tsteps, 1),matrix(3,tsteps+1,1), Sv1$meanf[,1], Sv2$meanf[,1],Sv3$meanf[,1], Sv4$meanf[,1]))
  Mmr <- as.data.frame(cbind(seq(0, tsteps, 1),matrix(4,tsteps+1,1), Sv1$meanf[,2], Sv2$meanf[,2],Sv3$meanf[,2], Sv4$meanf[,2]))
  Mhigh <- as.data.frame(cbind(seq(0, tsteps, 1),matrix(5,tsteps+1,1), Sv1$M[5,5,], Sv2$M[5,5,], Sv3$M[5,5,], Sv4$M[5,5,]))
  Musr <- rbind(Mu, Mb,Mmf,Mmr,Mhigh)
  colnames(Musr) <- c("t", "pop",1,2,3,4)
  Msrm <- melt(Musr, id.vars = c("t","pop"))
  facet_names <- c(`0` = "U", `1` = "B", `3` = "mean fit", `4` = "mean res", `5` = "Highest fit/res")
  ggplot(Msrm, aes(x=t, y = value, colour = variable)) + geom_line() + facet_wrap(~pop,labeller = as_labeller(facet_names), scales = "free")
  ggsave(paste(num,"_",pref,"norsd_TimeSeries_output_06.pdf",sep=""))
  # number in highest fitness changes but mean r and f don't? 
  
  # How is the bacterial load affected? want to control.
  colnames(Mb) <- c("t", "pop",1,2,3,4)
  Mbm <- melt(Mb, id.vars = c("t","pop"))
  pmb<-ggplot(Mbm, aes(x=t, y = value, colour = variable)) + geom_line() + scale_x_continuous(breaks = seq(0,tsteps,tsteps / 10), labels = dt*seq(0,tsteps,tsteps/10)) + 
    scale_y_continuous(lim = c(0,1000))
  ggsave(paste(num,"_",pref,"norsd_TimeSeries_Boutput_06.pdf",sep=""))
  
  # plots only proportion with "resistance" - not background initial resistance
  M10n <- Sv1$M;M20n <- Sv2$M;M30n <- Sv3$M;M40n <- Sv4$M
  susr <- wildtype[1]; susf <- wildtype[2]
  M10n[susr,susf,] <- 0; M20n[susr,susf,] <- 0; M30n[susr,susf,] <- 0; M40n[susr,susf,] <- 0; 
  for (i in 1:dim(Sv4$M)[3]){ M10n[,,i] <- M10n[,,i]/sum(M10n[,,i]) # re-normalise
  M20n[,,i] <- M20n[,,i]/sum(M20n[,,i]) # re-normalise
  M30n[,,i] <- M30n[,,i]/sum(M30n[,,i]) # re-normalise
  M40n[,,i] <- M40n[,,i]/sum(M40n[,,i])} # re-normalise
  
  # What happens? 
  mm1r<-c() ; mm2r<-c() ; mm3r<-c() ; mm4r<-c() 
  ll<-dim(Sv1$M)[3];
  ss<-seq(0,ll,1/dt) # Don't want to grab all 
  for(i in 2:length(ss)){
    mm01r<-as.data.frame(melt(M10n[,,ss[i]])); 
    mm02r<-as.data.frame(melt(M20n[,,ss[i]])); 
    mm03r<-as.data.frame(melt(M30n[,,ss[i]])); 
    mm04r<-as.data.frame(melt(M40n[,,ss[i]])); 
    mm01r$tstep=ss[i]*dt; mm02r$tstep=ss[i]*dt; mm03r$tstep=ss[i]*dt; mm04r$tstep=ss[i]*dt # To go to generations
    mm1r<-rbind(mm1r,mm01r);mm2r<-rbind(mm2r,mm02r); mm3r<-rbind(mm3r,mm03r) ; mm4r<-rbind(mm4r,mm04r) 
  } 
  colnames(mm1r)<-c("x","y","z","t"); colnames(mm2r)<-c("x","y","z","t"); colnames(mm3r)<-c("x","y","z","t") ;colnames(mm4r)<-c("x","y","z","t")
  setwd(plots)
  
  # plots
  p1<-ggplot(mm1r[w,],aes(x,y,fill=z))  + facet_wrap( ~ t, ncol=3) + ggtitle(paste("w = ", 1,sep=""))
  p1<-p1 + scale_fill_gradient("Proportion", limits=c(0,1),low="white", high="red",guide = FALSE)
  p1<-p1 + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p1<-p1 + geom_tile() + scale_y_continuous(breaks=c(1,nfit),"Relative fitness levels",labels=c("Least","Most")) + scale_x_continuous(breaks=c(mres,1),"Resistance levels",labels=c("Most","Least"))
  p1
  ggsave(paste(num,"_",pref,"norsdR_Array_w=",omega1,"_06.pdf",sep=""))
  p2<-ggplot(mm2r[w,],aes(x,y,fill=z))  + facet_wrap( ~ t, ncol=3) + ggtitle(paste("w = ", 2,sep=""))
  p2<-p2 + scale_fill_gradient("Proportion", limits=c(0,1),low="white", high="red",guide=FALSE)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p2<-p2 + geom_tile() + scale_y_continuous(breaks=c(1,nfit),"Relative fitness levels",labels=c("Least","Most")) + scale_x_continuous(breaks=c(mres,1),"Resistance levels",labels=c("Most","Least"))
  p2
  ggsave(paste(num,"_",pref,"norsdR_Array_w=",omega2,"_06.pdf",sep=""))
  p3<-ggplot(mm3r[w,],aes(x,y,fill=z))  + facet_wrap( ~ t, ncol=3) + ggtitle(paste("w = ", 3,sep=""))
  p3<-p3 + scale_fill_gradient("Proportion", limits=c(0,1),low="white", high="red",guide = FALSE)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p3<-p3 + geom_tile() + scale_y_continuous(breaks=c(1,nfit),"Relative fitness levels",labels=c("Least","Most")) + scale_x_continuous(breaks=c(mres,1),"Resistance levels",labels=c("Most","Least"))
  p3
  ggsave(paste(num,"_",pref,"norsdR_Array_w=",omega3,"_06.pdf",sep=""))
  p4<-ggplot(mm4r[w,],aes(x,y,fill=z))  + facet_wrap( ~ t, ncol=3) + ggtitle(paste("w = ", 4,sep=""))
  p4<-p4 + scale_fill_gradient("Proportion", limits=c(0,1),low="white", high="red",guide=FALSE)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p4<-p4 + geom_tile() + scale_y_continuous(breaks=c(1,nfit),"Relative fitness levels",labels=c("Least","Most")) + scale_x_continuous(breaks=c(mres,1),"Resistance levels",labels=c("Most","Least"))
  p4
  ggsave(paste(num,"_",pref,"norsdR_Array_w=",omega4,"_06.pdf",sep=""))
  
  setwd(plots)
  pdf(paste(num,"_",pref,"norsdR_Array_w_all.pdf",sep=""),width=18,height=18)
  multiplot(p1,p2,p4,p3,cols=2)
  dev.off()
  
  pdf(paste(num,"_",pref,"all_norsd_Array_w_all.pdf",sep = ""),width=18,height=18)
  multiplot(go,p1,p2,pmb,p4,p3,cols=2)
  dev.off()
  
  # Plot proportions in each fitness / resistance level over time
  pp<-c();
  ll<-dim(Sv1$M)[3];
  ss<-seq(0,ll,1/dt) # Don't want to grab all 
  for(i in 2:length(ss)){
    pp1<-c(ss[i]*dt,colSums(Sv1$M[,,ss[i]]), rowSums(Sv1$M[,,ss[i]]),omega1)
    pp2<-c(ss[i]*dt,colSums(Sv2$M[,,ss[i]]), rowSums(Sv2$M[,,ss[i]]),omega2)
    pp3<-c(ss[i]*dt,colSums(Sv3$M[,,ss[i]]), rowSums(Sv3$M[,,ss[i]]),omega3)
    pp4<-c(ss[i]*dt,colSums(Sv4$M[,,ss[i]]), rowSums(Sv4$M[,,ss[i]]),omega4)
    pp<-rbind(pp,pp1,pp2,pp3,pp4);
  } 
  pp<-as.data.frame(pp);
  if (mres == 5){
    colnames(pp)<-c("t","Fitness level 1\n(low)","Fitness level 2","Fitness level 3","Fitness level 4","Fitness level 5\n(high)",
                    "Res. level 1\n(low)","Res. level 2","Res. level 3","Res. level 4","Res. level 5\n(high)","w") 
  }
  if (mres == 10){
    colnames(pp)<-c("t","Fitness level 1\n(low)","Fitness level 2","Fitness level 3","Fitness level 4","Fitness level 5","Fitness level 6","Fitness level 7","Fitness level 8","Fitness level 9","Fitness level 10\n(high)",
                    "Res. level 1\n(low)","Res. level 2","Res. level 3","Res. level 4","Res. level 5","Res. level 6","Res. level 7","Res. level 8","Res. level 9","Res. level 10\n(high)","w"); }
  
  pp2<-melt(pp,id.vars = c("t","w"))
  g<-ggplot(pp2,aes(x=t,y=value,colour=factor(w))) + facet_wrap(~variable,ncol=mres) + geom_line(size=2) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  g<-g + scale_x_continuous("Generation") + scale_y_continuous("Proportion") #+ scale_colour_manual(values=cbPalette,"Abx\nlevel",labels=c(0.05,0.1,0.2))
  g
  ggsave(paste(num,"_",pref,"f&r_overtime.pdf",sep = ""),width=18,height=12)
  
}

### Function to set to zero if negative one otherwise
dirac_att   <- function(x) { x[x<0] <- 0  ;x[x>0] <- 1; x }
dirac_att_5 <- function(x) { x[x==0] <- 0.5; x }


######**********************************************************************************  *******#####
#sample from own distribution - decide who leaves. From http://stackoverflow.com/questions/12848736/how-to-declare-a-user-defined-distribution-in-r
rMydist_a <- function(n,acqdistn) {
  if(n==0){aa <- 0}else{
    nn <- dim(acqdistn)[1];
    mm <- dim(acqdistn)[2];
    # transform acqdistn into vector
    dim(acqdistn) <- NULL 
    ss <- sample(x = seq(1,length(acqdistn),1), size = n, prob = acqdistn, replace=T) # Sample with given distribution
    aa <- tabulate(ss,nbins=length(acqdistn)) 
    dim(aa) <- c(nn,mm)}
  return(aa) # Return matrix of distribution of where new are
}

# rMydist_a(10,rbind(c(0.1,0.3),c(0,0.6))) # e.g.
