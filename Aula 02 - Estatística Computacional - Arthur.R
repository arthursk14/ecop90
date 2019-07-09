#### Exercise 1 ###
### Superior limit estimation for a uniform distrbution ###

# Sets the seed of R random generator, for reproduction of this example 
set.seed(123)

# Install the necessary package (EnvStat)
install.packages("EnvStats")
library(EnvStats)

bml = c(1:10000)
bmm = c(1:10000)

for (i in 1:10e3) {
  # Creates the uniform distribution sample vector with n observations using b = 1
  n = 5
  s = runif(n, 0, 1)
  
  # Estimates using maximum likelihood and saves estimation of b in bml vector 
  ml = eunif(s, "mle")
  bml[i] = ml$parameters[2]
  
  # Estimates using moments and saves estimation of b in bmm vector
  mm = eunif(s, "mme")
  bmm[i] = mm$parameters[2]

}

# Boxplot of estimation
boxplot(bml-1, bmm-1)

#
#
#

### Professor Solution ###

nr = 10e3
n = c(5,10e0,10e1,10e3)
nn = length(n)

b = 1
eml = matrix(rep(0,nn*nr),ncol=nn)
emm = eml

set.seed(123)

for (i in 1:nn) {
  for (j in 1:nr){
    x = runif(n[i],0,b)
    eml[j,i] = max(x)
    emm[j,i] = 2*mean(x)
  }
  
}

eeml = eml - b
eemm = emm - b

boxplot(eeml[,1], eemm[,1], eeml[,2], eemm[,2], eeml[,3], eemm[,3], eeml[,4], eemm[,4])

ebml = apply(eeml,2,"mean")
ebmm = apply(eemm,2,"mean")

mseeml = apply(eeml^2,2,"mean")
mseemm = apply(eemm^2,2,"mean")

#
#
#

### Exercise 2 ###
### Variance estimation for a normal distribution and for a student's T distribution ###

# Sets the seed of R random generator, for reproduction of this example 
set.seed(123)

# Number of repetitions
nr = 10e3

# Sample sizes
n = c(5,20,100)
nn = length(n)

# Set paramaters of the normal distribution
m = c(0,10)
v = c(1,9)
l = length(m)
ll = l - 1

# Creates the matrices of results, i.e estimated variance for different sample sizes
# and using the defined different paramaters for the normal distribution
evs = matrix(rep(0,nn*nr*l),ncol=nn*l)
evm = evs

# Loop for estimation
for (i in 0:ll){
  for (j in 1:nn){
    for (k in 1:nr){
      
      x = rnorm(n[j],m[i+1],v[i+1]^(1/2))
      
      # Sample variance estimator
      es = (1/(n[j]-1))*sum((x-mean(x))^2)
      
      ji = j + 4*i
      
      evs[k,ji] = es
      
      # Method of moments
      evm[k,ji] = ((n[j]-1)/n[j])*es
      
    }
    
  }

}

# Calculate the error of the estimation

eevs = evs
eevm = evm

for (i in 0:ll){
  for (j in 1:nn){
    ji = j + 4*i
    eevs[,ji] = evs[,ji] - v[i+1]
    eevm[,ji] = evm[,ji] - v[i+1]
  
  }
}


# Boxplot of estimation error
boxplot(eevs[,1], eevm[,1], eevs[,2], eevm[,2], eevs[,3], eevm[,3], eevs[,4], eevm[,4],
        eevs[,5], eevm[,5], eevs[,6], eevm[,6], eevs[,7], eevm[,7], eevs[,8], eevm[,8], outline = F)

#
#
#

### Professor Solution (for the student's T distribution with v degrees of freedom) ###

nr = 10e3
n = c(5,20,100)
nn = length(n)
v = 4

s2 = matrix(rep(0,nn*nr),ncol=nn)
m2 = s2

set.seed(123)

for (i in 1:nn){
  for (j in 1:nr){
    x = rt(n[i],v)
    s2[j,i] = var(x)
    m2[j,i] = var(x)*((n[i]-1)/n[i])
    
  }
  
}

es2 = s2 - (v/(v-2))
em2 = m2 - (v/(v-2))

boxplot(es2[,1], em2[,1], es2[,2], em2[,2], es2[,3], em2[,3], outline = F)
abline(h=0)

bes2 = apply(es2,2,mean)
medes2 = apply(es2,2,median)
bem2 = apply(em2,2,mean)
medem2 = apply(em2,2,median)

del = density(es2[,1])
plot(del)

#
#
#