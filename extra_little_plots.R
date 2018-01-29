#### CHECKS 

### SUB-MIC PLOTS
##### LINEAR MIC WHAT IS THE SELECTION? 

m1 <- pmax(0,(10 - seq(1,30,1))/10)
m2 <- pmax(0,(20 - seq(1,30,1))/20)
m3 <- pmax(0,(30 - seq(1,30,1))/30)

plot(m1, type = "l")
lines(m2,col="red")
lines(m3,col="blue")

# Difference - what selection acts on
plot(m3 - m1, type = "l")
lines(m2 - m1,col="red")

##### STEPPED MIC WHAT IS THE SELECTION? 

m1 <- c(matrix(1,9,1),matrix(0,21,1))
m2 <- c(matrix(1,19,1),matrix(0,11,1))
m3 <- c(matrix(1,29,1),matrix(0,1,1))

plot(m1, type = "l")
lines(m2,col="red")
lines(m3,col="blue")

# Difference - what selection acts on
plot(m3 - m1, type = "l")
lines(m2 - m1,col="red")
