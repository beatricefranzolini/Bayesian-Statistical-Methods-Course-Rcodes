rm(list = ls())

## SIMULATE A DATASET
set.seed(0)
n = 100
z_true = rbinom(n, 1, 0.5) #to simulate from a bernoulli we use the bin with param 1 
beta_true = c(1, 4)
y = rep(NA, n)
y[(z_true==0)] = rnorm(sum(z_true==0), beta_true[1], 1)
y[(z_true==1)] = rnorm(sum(z_true==1), beta_true[2], 1)
hist(y)


## SET PRIOR HYPERPARAMETERS
meanb = mean(y); varb = 1

# BUILD EMPTY MATRIX OF POSTERIOR SAMPLES
S = 10000
posterior.samples.Z = matrix(NA, nrow = S, ncol = n)
posterior.samples.beta = matrix(NA, nrow = S, ncol = 2)
posterior.samples.pi = matrix(NA, nrow = S, ncol = 1)

# initialize unknown parameters z, beta and pi
z = rbinom(n, 1, 0.5) #cluster them randomly to start
beta = rnorm(2, meanb, sqrt(varb))
pi = runif(1)

#save initialization
t = 1
posterior.samples.Z[t,] = z
posterior.samples.beta[t,] = beta
posterior.samples.pi[t,] = pi

#the following lines create a progress bar 
pb <- txtProgressBar(min = 1,      # Minimum value of the progress bar
                     max = S, # Maximum value of the progress bar
                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
                     width = 50,   # Progress bar width. Defaults to getOption("width")
                     char = "=")   # Character used to create the bar

for (t in 2:S){
  
  #sample z #####################################################
  
  for(i in 1:n){
    
    #compute parameter of full conditional of z 
    p1 = pi * dnorm(y[i], beta[1], 1)
    p2 = (1 - pi) * dnorm(y[i], beta[2], 1)
    
    #sample z from full conditional
    z[i] = sample(c(0,1), 1, prob = c(p1, p2))
    
  }
  
  #sample beta #################################################
  
  #compute parameter of full conditional of beta0
  ybar0 = ifelse(is.na(mean(y[z==0])), 0, mean(y[z==0]))
  mean_beta0 = (1 / varb / (1 / varb + sum(z==0))) * meanb + 
      ( sum(z==0) / (1 / varb + sum(z==0)) ) * ybar0
  var_beta0 = 1 / (1 /varb + sum(z==0))
  
  #sample beta0 from full conditional 
  beta[1] = rnorm(1, mean_beta0, sqrt(var_beta0))
  
  #compute parameter of full conditional of beta1
  ybar1 = ifelse(is.na(mean(y[z==1])), 0, mean(y[z==1]))
  mean_beta1 = (1 / (1 + sum(z==1))) * meanb + 
      ( sum(z==1) / (1 +  sum(z==1)) ) * ybar1
  var_beta1 = 1 / (1 + sum(z==1))
  
  #sample beta1 from full conditional 
  beta[2] = rnorm(1, mean_beta1, sqrt(var_beta1))
  
  #sample pi #################################################
  
  #compute parameter of full conditional of pi
  a = 1 + sum(z)
  b = 1 + n - sum(z)
  pi = rbeta(1, a, b)
  
  #save chain value 
  posterior.samples.Z[t,] = z
  posterior.samples.beta[t,] = beta
  posterior.samples.pi[t,] = pi
  
  setTxtProgressBar(pb, t)
}
close(pb)

burnin = 2000

#plot one realization of z
hist(y)
stripchart(y,
           method = "jitter",
           pch = 23,
           bg = posterior.samples.Z[7000,],
           add = TRUE)

#compute posterior coclustering matrix 
library(T4cluster) #library to compute the posterior coclustering matrix
library(plot.matrix) #library to plot a matrix
PSM_estimate = psm(posterior.samples.Z[(burnin+1):S,])
plot(PSM_estimate)
PSM_true = psm(matrix(z_true, nrow = 1, ncol = n))
plot(PSM_true)

plot(posterior.samples.beta[(burnin+1):S, 1])
plot(posterior.samples.beta[(burnin+1):S, 2])
plot(posterior.samples.pi[(burnin+1):S, 1])

#to derive posterior point estimates we have to deal with label switching(!)... 
#upcoming! 
