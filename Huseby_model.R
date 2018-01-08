### Model from Huseby paper 2017

# Baseline parameters
k = 500 #mg / L = concentration of resource at which max f is half value (Monod constant)
pop_fit_absence <- c(0.6,0.7,0.8,0.9,1)
pop_mic         <- c(0.25,0.2,0.15,0.1,0.05)
en <- 2.5 * 10^(-8) # mu/cell growth 

# Growth rate
phi <- function(r,f){
  phi_est <- f * (r / (r + k))
  return(phi_est)
}

# Relative fitness of bacteria
fitness <- function(i,ab){
  f_absence <- pop_fit_absence[i]
  mic <- pop_mic[i]
  if(ab < mic){fit_est <- f_absence * ((mic - ab)/mic)
  }else{fit_est <- 0}
  return(fit_est)
} 

# outputs
mtstep = 100
n = matrix(0,5,mtstep)
r = matrix(0,1,mtstep)

# initial conditions
n[,1] = matrix(0.2,5,1)*10^6
r[1] = 100

###### Difference equations #######
for(t in 2:mtstep){
  ab = 0.04
  
  for(i in 1:5){
  n[i,t] = n[i,t-1] * exp(fitness(i,ab)*R[t-1]/(R[t-1]+k))}
  
  R[t] = R[t-1] * (1 - sum(n[,t-1])*en)
  
}
par(mfrow = c(1, 2))
plot(seq(1,mtstep,1),n[3,],type = "l")
lines(seq(1,mtstep,1),n[2,],type = "l",col="blue")
lines(seq(1,mtstep,1),n[3,],type = "l",col="green")
lines(seq(1,mtstep,1),n[4,],type = "l",col="pink")
lines(seq(1,mtstep,1),n[5,],type = "l",col="red")
lines(seq(1,mtstep,1),n[1,],type = "l")
plot(seq(1,mtstep,1),R, type = "l")


# ##### ODE
# ## Integrals too difficult for R... 
# for(t in 2:(mtstep)){
#   ab = 0.004
#   # density of population i
#   for(i in 1:5){
#     n[i,t] = exp(t*phi(r[(t-1)], fitness(i,ab))) + (n[i,1] - 1)
#   }
#   
#   # amount of resource...
#   # bexp(r)r = -aexp(At)... 
#   sum_r_use = 0
#   for(i in 1:5){ 
#     a =  - n[i,t]*en
#     sum_r_use = sum_r_use - a*r[t] - log(r[t] + k) + (r[1] + n[i,1]*en*r[1] + log(r[1]+k))
#   }
#   r[t] = sum_r_use
# }
# 
# plot(seq(1,length(n[1,]),1),n[1,])
# lines(seq(1,length(n[1,]),1),n[2,],col='blue')
# lines(seq(1,length(n[1,]),1),n[3,],col='green')
# lines(seq(1,length(n[1,]),1),n[4,],col='pink')
# lines(seq(1,length(n[1,]),1),n[5,],col='red')
