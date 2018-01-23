### Original EC model dynamics
# do they match those of huseby without any change in structure? 

# original:
#U[i+1] =  U[i] + mu*(S[i]+R[i]) - (lambdas+lambdar)*U[i] + omega*ks*S[i] + omega*kr*R[i]
#S[i+1] =  S[i] + lambdas*U[i] - (mu + omega*ks)*S[i] - eps * S[i]
#R[i+1] =  R[i] + lambdar*U[i] - (mu + omega*kr)*R[i] + eps * S[i] 

##### (1) is the growth rate / relative fitness a "linear decreasing function of antibiotic
# concentration? 
### EC has 
#R[i+1] =  R[i] (1 + X$meanfit * beta/N * U[i] - (mu + omega*kr)) + eps * S[i] 
# growth rate = X$meanfit * beta/N * U[i] - (mu + omega*kr)
# call meanfit = mean fitness in absence of antibiotics 
# U = nutrient 
# omega = antibiotic concentration
# kr = MIC dependent
# then yes gr = linear decline as omega increases

##### (2) bacterial growth rate a monotonically increasing function of nutrient?
### EC has U = space = nutrient... 
# need to add a limiting factor - less space, harder to reach 
# lambda * U/ (U + k)?????


gr = matrix(0,1,100)
R[1] = 100

meanfit = 0.9
beta = 20
N = 1
U = 1
S = 0
mu = 0
omega = 0.6
kr = 1
eps = 0.00008
for(i in 2:100){
  
  R[i] =  R[i-1] * (1 + meanfit * beta/N * U - (mu + omega*kr)) + eps * S
  
}
plot(seq(0,100,1),R)
