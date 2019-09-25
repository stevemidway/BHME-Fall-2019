# rm(list=ls())
# Install packages
library(R2jags)
library(arm)
library(lattice)

# Read in data
dat <- read.table('PLD.txt',na.strings='NA',header=T)

# Look at structure of the data


# Check first 6 rows of data


# Check dimensions of data


# Obtain summary of data


# Plot PLD temp relationship
plot(pld ~ temp, data=dat, xlab='Temperature (C)', ylab='PLD (days)',las=1)
abline(lm(pld ~ temp,data=dat))

# Fit linear model and check diagnostics
m1 <- lm(pld ~ temp, data=dat)
summary(m1)

# Diagnostic plots
par(mfrow=c(2,2))
plot(m1)

# Histogram of residuals


# Assumption of homoscedasticity is violated, we will work with log(temp) and log(pld)


###########################################
# Log-transform the pld and temp data
###########################################
m2 <- lm(log(pld) ~ log(temp), data=dat)
par(mfrow=c(2,2))
plot(m2)
hist(resid(m2),breaks=20)

# Residual plots look better

# Summary of model output
summary(m2)

# Quick plot data using log(pld) and log(temp)
plot(log(pld)~log(temp),data=dat,ylab='log(PLD)',xlab='log(Temperature (C))',pch=16)
abline(lm(log(pld)~log(temp),data=dat))


####### FIT MODEL IN JAGS################
#########################################

# Define the model in the BUGS language and write a text file
sink("linReg.txt")
cat("
model {

# Likelihood: 
for (i in 1:n){                           # loop over observation i to n
                  # 1. Distribution for random part
                  # 2. Linear predictor

 # Calculate residuals
   resid[i] <- y[i] - mu[i]

   }  # end i loop

# Derived quantities
tau <- pow(sigma,-2)
sigma2 <- pow(sigma,2)

# Priors



} # end model
",fill = TRUE)
sink()



############################
# Bundle data
data <- list(y= , 
             x= , 
             n= )

# Initial values

inits <- function(){list(sigma = runif(1, 0, 10), alpha=rnorm(1), beta=rnorm(1))}

# Parameters to estimate/keep track of
parameters <- c("alpha","beta","sigma","resid", "mu")

# MCMC settings
niter <- 8000
nthin <- 3
nburn <- 3000
nchains <- 3


# Do the MCMC stuff calling JAGS from R
out <- jags(data = data, inits = inits, parameters.to.save = parameters, 
model.file = "linReg.txt", n.chains = nchains, n.thin = nthin, n.iter = niter, 
n.burnin = nburn)


# Summarize the result
print(out, digits = 3)

# Look at structure of object out
# str(out)


######## Diagnostic Plots #############


# Grab all MCMC samples for plotting
out.array <- out$BUGSoutput$sims.array
output <- rbind(out.array[,1,],out.array[,2,],out.array[,3,]) #convert to matrix; the matrix directly from out output (out$BUGSoutput$sims.matrix) does not sort by chain

# Number of samples retained
n.keep <-out$BUGSoutput$n.keep

# Trace plots


# Density plots




######## GRAPH RESULTS #############





