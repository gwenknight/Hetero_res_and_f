#### CHECKS 

### SUB-MIC PLOTS
##### LINEAR MIC WHAT IS THE SELECTION? 

m1 <- pmax(0,(10 - seq(0,30,1))/10)
m2 <- pmax(0,(20 - seq(0,30,1))/20)
m3 <- pmax(0,(30 - seq(0,30,1))/30)

setwd("~/Documents/Hetero_res_and_f/plots/")
pdf("MIC_linear_decline.pdf")
plot(m3, type = "l",ylim = c(0,1),xlab = "[Antibiotic]",ylab = "Relative fitness")
lines(m2,col="red")
lines(m1,col="blue")
abline( v = 30, lty = 2)
abline( v = 20, col="red",lty = 2)
abline( v = 10, col="blue",lty = 2)
dev.off()

# Difference - what selection acts on
plot(m3 - m1, type = "l")
lines(m2 - m1,col="red")

##### STEPPED MIC WHAT IS THE SELECTION? 

m1 <- c(matrix(1,9,1),matrix(0,21,1))
m2 <- c(matrix(1,19,1),matrix(0,11,1))
m3 <- c(matrix(1,29,1),matrix(0,1,1))

pdf("MIC_stepped_decline.pdf")
plot(m3, type = "l",ylim = c(0,1),xlab = "[Antibiotic]",ylab = "Relative fitness")
lines(m2,col="red")
lines(m1,col="blue")
abline( v = 30, lty = 2)
abline( v = 20, col="red",lty = 2)
abline( v = 10, col="blue",lty = 2)
dev.off()

# Difference - what selection acts on
plot(m3 - m1, type = "l")
lines(m2 - m1,col="red")

##### HALF AT MIC

m1 <- c(matrix(1,9,1),matrix(0.5,25,1))
m2 <- c(matrix(1,19,1),matrix(0.5,15,1))
m3 <- c(matrix(1,29,1),matrix(0.5,5,1))

pdf("MIC_half_decline.pdf")
plot(m3, type = "l",ylim = c(0,1),xlab = "[Antibiotic]",ylab = "Relative fitness")
lines(m2,col="red")
lines(m1,col="blue")
abline( v = 30, lty = 2)
abline( v = 20, col="red",lty = 2)
abline( v = 10, col="blue",lty = 2)
dev.off()
# Difference - what selection acts on
plot(m3 - m1, type = "l")
lines(m2 - m1,col="red")

