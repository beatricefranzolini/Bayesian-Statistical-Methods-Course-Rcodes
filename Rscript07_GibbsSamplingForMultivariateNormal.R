## IN THIS CODE WE PERFORM GIBBS SAMPLING ON A MULTIVARIATE NORMAL 
## MODEL WITH UNKNOWN MEAN AND VARIANCE-COVARIANCE MATRIX SEE ALSO CHAPTER 7 
## OF HOFF'S BOOK FOR MORE DETAILS ON DERIVATION OF FULL CONDITIONALS AND R CODE

#install.packages("mvtnorm")
require(mvtnorm)
##  PROVIDE THE MATRIX OF OBSERVED DATAPOINTS (IN THIS 
## CASE WE ARE USING "SYNTETHIC" DATA). WE HAVE n OBSERVATIONS OF DIMENSION d
n<-100
d<-2
x<-matrix(rnorm(d*n,mean = 3.2,sd = sqrt(1.6)),nrow = n,ncol = d)
mean.x<-colMeans(x)

##  DEFINE THE PRIOR HYPERPARAMETERS
mu0<-rep(0,d)
L0<-diag(100,nrow = d,ncol = d)
nu0<-1
S0<-diag(1,nrow = d,ncol = d)

## BUILD AN EMPTY MATRIX OF POSTERIOR SAMPLES
s<-1000
posterior.samples<-matrix(NA,nrow = s,ncol = d+d^2)
if(d==2){colnames(posterior.samples)<-c("mu1","mu2","Sigma11","Sigma12","Sigma21","Sigma22")}

## INITIALIZE THE MARKOV CHAIN (mu,Sigma) AT SOME ARBITRARY VALUES IN THE PARAMETRIC SPACE
mu<-rep(1,d)
Sigma<-diag(1,nrow = d,ncol = d)
t<-1
posterior.samples[t,]<-c(mu,Sigma)

for (t in 2:s){# SIMULATE THE MARKOV CHAIN (mu,Sigma) FROM t=1 to t=s
  
  # COMPUTE THE PARAMETERS OF THE FULL CONDITIONAL OF mu
  Ln<-solve( solve(L0) + n*solve(Sigma) )
  mun<-Ln%*%( solve(L0)%*%mu0 + n*solve(Sigma)%*%mean.x )
  # SIMULATE A NEW mu FROM ITS FULL CONDITIONAL
  mu<-rmvnorm(n = 1,mean = mun,sigma = Ln)  
  
  # COMPUTE THE PARAMETERS OF THE FULL CONDITIONAL OF Sigma
  ###update Sigma
  nun<-nu0+n
  Sn<- S0 + ( t(x)-c(mu) )%*%t( t(x)-c(mu) ) 
  # SIMULATE A NEW Sigma FROM ITS FULL CONDITIONAL
  Sigma<-solve( rWishart(n = 1, df = nun, Sigma = solve(Sn))[,,1] )
  
  # STORE THE NEW VALUE OF THE MARKOV CHAIN IN THE MATRIX OF POSTERIOR SAMPLES
  posterior.samples[t,]<-c(mu,c(Sigma))
}

print("Posterior expected values:")
print(colMeans(posterior.samples))
print("Posterior standard deviations:")
print(apply(posterior.samples,MARGIN = 2,sd))

par(mfrow=c(2,1))
plot(posterior.samples[,1],type="l",main="trajectory of mu_1^(t)")
plot(posterior.samples[,2],type="l",main="trajectory of mu_2^(t)")
par(mfrow=c(3,1))
plot(posterior.samples[,3],type="l",main="trajectory of Sigma_11^(t)")
plot(posterior.samples[,6],type="l",main="trajectory of Sigma_22^(t)")
plot(posterior.samples[,4],type="l",main="trajectory of Sigma_12^(t)")
par(mfrow=c(1,1))

plot(posterior.samples[-c(1:10),1],posterior.samples[-c(1:10),2],
     main="joint distribution of (mu_1,mu_2)")



### DO ANY INFERENCE OF INTEREST USING THE SAMPLES ....
require(coda)
samples.mcmc<-as.mcmc(posterior.samples[-c(1:10),-4])
summary(samples.mcmc)
plot(samples.mcmc)
crosscorr(samples.mcmc)
