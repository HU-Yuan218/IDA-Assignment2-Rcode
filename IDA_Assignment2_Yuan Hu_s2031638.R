### IDA Assignment 2
### Yuan Hu s2031638

### Q2
## (b)
install.packages("maxLik")
#load("dataex2.Rdata")
load("~/Desktop/Semester 1/Incomplete Data Analysis (IDA)/Assignment2/dataex2.Rdata")

# sigma is known and equal to 1.5
log_like_normal = function(param, data){
  x = data$X
  r = data$R
  mu = param
  sum(r*dnorm(x, mean = mu, sd = 1.5, log = TRUE) + (1-r)*pnorm(x, mean = mu, sd = 1.5, lower.tail = TRUE, log.p = TRUE))
}

library(maxLik)
require(maxLik)
# Use the 'maxLik' function, using as starting values 5 (mu).
mle = maxLik(logLik = log_like_normal, data = dataex2, start = c(5))
summary(mle)

# Alternatively, use the 'optim' function, choosing method = "Brent", and using as starting values 10 (mu).
mleoptim_brent <- optim(par = c(10), fn = log_like_normal, data = dataex2, method = "Brent", lower = -10, upper = 10, control = list("fnscale"=-1), hessian = TRUE)
mleoptim_brent



### Q4
# load("dataex4.Rdata")
load("~/Desktop/Semester 1/Incomplete Data Analysis (IDA)/Assignment2/dataex4.Rdata")
em.bernoulli = function(data, beta, eps){
  x = data$X
  y = data$Y
  ind_yobs = which(is.na(y) == FALSE)
  y_obs = y[ind_yobs]
  x_obs = x[ind_yobs]
  x_mis = x[-ind_yobs]
  diff = 1
  while(diff > eps){
    beta_t.old = beta
    beta0_t = beta_t.old[1]
    beta1_t = beta_t.old[2]
    # E-step
    log_like_Q = function(beta, data){
      beta0 = beta[1]
      beta1 = beta[2]
      log_like = -sum(log(1 + exp(beta0 + x*beta1))) + sum(y_obs*(beta0 + x_obs*beta1)) + sum((beta0 + x_mis*beta1)*(exp(beta0_t + x_mis*beta1_t)/(1 + exp(beta0_t + x_mis*beta1_t))))
      return(log_like)
    } 
    # M-step
    beta_t.1 = maxLik(logLik = log_like_Q, data = data, start = c(beta))
    beta = c(beta_t.1$estimate[1], beta_t.1$estimate[2])
    diff = sum(abs(beta - beta_t.old))
  }
  return(beta)
}
em.bernoulli(data = dataex4, beta = c(1,1), eps = 0.00001)



### Q5
## (b)
em.mixture = function(y, theta0, eps){
  n = length(y)
  theta = theta0
  p = theta[1]
  mu = theta[2]
  sigma = theta[3]
  lambda = theta[4]
  diff = 1
  while(diff > eps){
    theta.old = theta
    # E-step
    ptilde1 = p*dlnorm(y, meanlog = mu, sdlog = sigma)
    ptilde2 = (1 - p)*dexp(y, rate = lambda)
    ptilde = ptilde1/(ptilde1 + ptilde2)
    # M-step
    p = mean(ptilde)
    mu = sum(log(y)*ptilde)/sum(ptilde)
    sigma = sqrt(sum(((log(y) - mu)^2)*ptilde)/sum(ptilde))
    lambda = sum(1 - ptilde)/sum(y*(1 - ptilde))
    theta = c(p, mu, sigma, lambda)
    diff = sum(abs(theta - theta.old))
  }
  return(theta)
}

# load("dataex5.Rdata")
load("~/Desktop/Semester 1/Incomplete Data Analysis (IDA)/Assignment2/dataex5.Rdata")

res = em.mixture(y = dataex5, c(0.1, 1, 0.5, 2), 0.00001)
p = res[1]
mu = res[2]
sigma = res[3]
lambda = res[4]

p; mu; sigma; lambda

hist(dataex5, main = "Histogram of dataex5 with estimated density superimposed", xlab = "dataex5", ylab = "Density", cex.main = 1, cex.lab = 0.8, cex.axis = 0.8, freq = F, xlim = c(0, 140), ylim = c(0, 0.1))
curve(p*dlnorm(x, mu, sigma) + (1 - p)*dexp(x, lambda), add = TRUE, lwd = 2, col = "blue2")



