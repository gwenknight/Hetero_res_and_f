trans = 3 
trate = 3
m = 5
n = 5

vf<-seq(1/n,1,1/n)

vs_mic <- seq(0.1,32,length.out = m) 
vs_rel <- pmax((vs_mic - trate),0)/vs_mic # constant ### SHOULD SUM TO 1? 

M_new <- matrix(0.1,5,5)
M_new[,1] <- 0.1
M_new[,2] <- 0.2
M_new[,3] <- 0.3
M_new[,4] <- 0.4
M_new[,5] <- 0.5
M_new <- M_new/sum(M_new)

# Varying fitness levels but resistance at start all same
colSums(M_new)
rowSums(M_new)
vf
vs_rel

Mf<-colSums(M_new) # Proportions at each FITNESS level
if(trans>0){ pastmean=sum( Mf*vf ) }else{ pastmean = 1 }
M_temp<-M_new
#for(i in 1:length(Mf)){M_temp[,i] = M_temp[,i] * vf[i] / pastmean } # Updated matrix: colSums(M_new) = Mf_new
M_temp <- t(t(M_temp) * vf / pastmean)

colSums(M_new) # should have changed
colSums(M_temp)
rowSums(M_new) # should be same
rowSums(M_temp)


# RESISTANCE
Mr<-rowSums(M_temp) # = rowSums(M_new) same # Proportions at each RESISTANCE level
if(trans>0){ pastmean=sum( Mr*vs_rel ) }else{ pastmean = 1 }
# update M_temp with fitness then resistance
#for(i in 1:length(Mf)){M_temp[i,] = M_temp[i,] * vs_rel[i] / pastmean } # Updated matrix: colSums(M_new) = Mf_new
M_temp <- vs_rel*M_temp/pastmean

colSums(M_new) # should have changed
colSums(M_temp)
rowSums(M_new) # should now have changed
rowSums(M_temp)


#####****
M_new <- matrix(0.1,5,5)
M_new[1,] <- 0.1
M_new[2,] <- 0.2
M_new[3,] <- 0.3
M_new[4,] <- 0.4
M_new[5,] <- 0.5
M_new <- M_new/sum(M_new)

# Varying fitness levels but resistance at start all same
colSums(M_new)
rowSums(M_new)
vf
vs_rel

Mf<-colSums(M_new) # Proportions at each FITNESS level
if(trans>0){ pastmean=sum( Mf*vf ) }else{ pastmean = 1 }
M_temp<-M_new
#for(i in 1:length(Mf)){M_temp[,i] = M_temp[,i] * vf[i] / pastmean } # Updated matrix: colSums(M_new) = Mf_new
M_temp <- t(t(M_temp) * vf / pastmean)

colSums(M_new) # should have changed
colSums(M_temp)
rowSums(M_new) # should be same
rowSums(M_temp)


# RESISTANCE
Mr<-rowSums(M_temp) # = rowSums(M_new) same # Proportions at each RESISTANCE level
if(trans>0){ pastmean=sum( Mr*vs_rel ) }else{ pastmean = 1 }
# update M_temp with fitness then resistance
#for(i in 1:length(Mf)){M_temp[i,] = M_temp[i,] * vs_rel[i] / pastmean } # Updated matrix: colSums(M_new) = Mf_new
M_temp <- vs_rel*M_temp/pastmean

colSums(M_new) # should have changed
colSums(M_temp)
rowSums(M_new) # should now have changed
rowSums(M_temp)
