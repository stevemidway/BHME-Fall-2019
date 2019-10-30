# rm(list=ls())

library(car)
library(R2jags)
library(lattice)
library(MCMCpack)
library(coda)


## Data Generation
y <- c(1,8,10,9,8,7.5,7,6.75,6.5,6.4)
x <- c(1,2,3,4,5,6,7,8,9,10)

y <- c(y+rnorm(1, 0, 0.25), 
       y + 3+rnorm(1,0, 0.25) ,
       y + 6+rnorm(1,0, 0.25) )
x <- c(x, x, x)
g <- c(rep(1,times=10),
       rep(2,times=10),
       rep(3,times=10))
df <- cbind(y, x, g)
df <- as.data.frame(df)
plot(y ~ x, col = g, pch=16)

## nls() appraoch
h4 <- nls(y ~ (a*x^2)/(b + c*x + x^2),start = list(a = 5, b = 1, c = 1)) 
summary(h4)

x2 <- seq(1,10,length.out = 100)
y2 <- (coef(h4)[1] * x2^2)/(coef(h4)[2] + coef(h4)[3] * x2 + x2^2)

points(y2 ~ x2, type="l", col="blue", lwd=3)
text(7, 15, "Hollings Type IV Response: y = a*x^2/b + c*x + x^2")
abline(v = (-2*(coef(h4)[2]) / (coef(h4)[3])), lty=2, lwd=3, col="darkgray")
text(4,5,"Peak = (-2*b)/c")

# -2b/c = peak

###################### 
sink("holling4.txt")
cat("
model {
  # Likelihood: 
  # Level-1 of the model
  for (i in 1:n){ 
    y[i] ~ dnorm(y.hat[i], tau.y)               
    y.hat[i] <- a[group[i]] * x[i]^2 / (b[group[i]] + c[group[i]]*x[i] + x[i]^2) 
  } 
  
  tau.y <- pow(sigma.y,-2)
  sigma.y ~ dunif(0,100)
  
  # Level-2 of the model
  for (j in 1:J){
    a[j] <- BB[j,1]
    b[j] <- BB[j,2]
    c[j] <- BB[j,3]

    BB[j,1:3] ~ dmnorm(BB.hat[j,], Tau.B[,])

    BB.hat[j,1] <- mu.a 
    BB.hat[j,2] <- mu.b 
    BB.hat[j,3] <- mu.c

    peak[j] <- -2*b[j]/c[j]
  }

  # Priors and derived quantities
  mu.a ~ dnorm(0, 0.0001)
  mu.b ~ dnorm(0, 0.0001)
  mu.c ~ dnorm(0, 0.0001)
  
  ### Model variance-covariance
  Tau.B[1:K,1:K] ~ dwish(W[,], df)
  df <- K+1
  Sigma.B[1:K,1:K] <- inverse(Tau.B[,])
  for (k in 1:K){
    for (k.prime in 1:K){
      rho.B[k,k.prime] <- Sigma.B[k,k.prime]/
        sqrt(Sigma.B[k,k]*Sigma.B[k.prime,k.prime])
    }
    sigma.B[k] <- sqrt(Sigma.B[k,k])
  }
  
} # end model
",fill = TRUE)
sink()

n <- nrow(df) #Number of observations

J <- length(unique(df$g))

# Create identity matrix for Wishart dist'n
# Number of parameters to estimate (K)
K <- 3

W <- diag(K)


# load data
data <- list(y = df$y, 
             x = df$x, 
             group = as.numeric(df$g),
             n = n,
             J = J, 
             K = K, 
             W = W)


# Initial values
inits <- function (){
  list (mu.a = rnorm(1), mu.b=rnorm(1), mu.c=rnorm(1),
        sigma.y=runif(1),
        BB=matrix(rnorm(J*K),nrow=J,ncol=K),Tau.B=rwish(K+1,diag(K)) )
}


# Parameters monitored
params1 <- c("BB","mu.a","mu.b","mu.c",
             "sigma.B","Tau.B", "peak")

# MCMC settings
ni <- 100000
nt <- 4
nb <- 50000
nc <- 3


out1 <- jags(data, inits, params1, "holling4.txt", n.chains = nc, 
             n.thin = nt, n.iter = ni, n.burnin = nb)

system("say Happy Birthday Shane Flinn. I fit this awesome non-linear model for you. Lets cross our fingers for some convergence.")

# Summarize posteriors
print(out1, dig = 3)

plot(out1)
# Create MCMC object so you can use the CODA package to manipulate and examine data
out.mcmc <- as.mcmc(out1)
# Create traceplots
xyplot(out.mcmc[,1:6])
# Look at posterior density plots
densityplot(out.mcmc[,1:6])


## Plot JAGS results by group
x.fake <- seq(1,10,length.out = 100)
y.g1 <- (out1$BUGSoutput$mean$BB[1,1] * x.fake^2) / (out1$BUGSoutput$mean$BB[1,2] + out1$BUGSoutput$mean$BB[1,3]*x.fake + x.fake^2) 
y.g2 <- (out1$BUGSoutput$mean$BB[2,1] * x.fake^2) / (out1$BUGSoutput$mean$BB[2,2] + out1$BUGSoutput$mean$BB[2,3]*x.fake + x.fake^2) 
y.g3 <- (out1$BUGSoutput$mean$BB[3,1] * x.fake^2) / (out1$BUGSoutput$mean$BB[3,2] + out1$BUGSoutput$mean$BB[3,3]*x.fake + x.fake^2) 

# Add curves
plot(y.g1 ~ x.fake, pch=16, ylim=c(0,20), type="l", lwd=4, las=1, ylab="y", col="blue")
points(y.g2 ~ x.fake, pch=16, ylim=c(0,20), type="l", lwd=4, col="red")
points(y.g3 ~ x.fake, pch=16, ylim=c(0,20), type="l", lwd=4, col="green")

# Add peaks with uncertainty
segments(out1$BUGSoutput$summary[23,3],max(y.g1),out1$BUGSoutput$summary[23,7],max(y.g1), lwd=6)
segments(out1$BUGSoutput$summary[24,3],max(y.g2),out1$BUGSoutput$summary[24,7],max(y.g2), lwd=6)
segments(out1$BUGSoutput$summary[25,3],max(y.g3),out1$BUGSoutput$summary[25,7],max(y.g3), lwd=6)
points(out1$BUGSoutput$mean$peak[1],max(y.g1), cex=1, pch=16, col="gray")
points(out1$BUGSoutput$mean$peak[2],max(y.g2), cex=1, pch=16, col="gray")
points(out1$BUGSoutput$mean$peak[3],max(y.g3), cex=1, pch=16, col="gray")
