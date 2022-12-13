## ----r------------------------------------------------------------------------
library(mice)       # Just used for the nhanes2 data set
library(INLA)       # INLA modelling
library(dplyr)      # Data wrangling of the results
library(gt)         # Tables
library(tidyverse)  # Data wrangling and plotting
library(showtext)   # Font
library(colorspace) # Color adjustments
library(MCMCpack)   # dinvgamma()


## ----r------------------------------------------------------------------------
inla.setOption(num.threads = "1:1")


## ----r------------------------------------------------------------------------
# Using the nhanes data set found in mice so we can compare to mice:
data(nhanes2)

head(nhanes2)

n <- nrow(nhanes2)

# Manually dummy-code age:
age2 <- ifelse(nhanes2$age == "40-59", 1, 0)
age3 <- ifelse(nhanes2$age == "60-99", 1, 0)

# Center the response and continuous covariates
chl <- scale(nhanes2$chl, scale = FALSE)[,1]
bmi <- scale(nhanes2$bmi, scale = FALSE)[,1]


## ----r------------------------------------------------------------------------
# Priors for model of interest coefficients
prior.beta = c(0, 1e-6) # Gaussian, c(mean, precision)

# Priors for exposure model coefficients
prior.alpha <- c(0, 1e-6) # Gaussian, c(mean, precision)

# Priors for y, measurement error and true x-value precision
# Start by getting a reasonable prior guess for the standard error of the regression and exp. models
summary(lm(chl~bmi+age2+age3))$sigma
summary(lm(bmi~age2+age3))$sigma

# Use those values to create reasonable priors:
prior.prec.y <- c(2.5-1,(2.5)*29.1^2) # Gamma
prior.prec.u_c <- c(0.5, 0.5) # Gamma
prior.prec.x <- c(2.5-1,(2.5)*4.2^2) # Gamma

# We can visualize these priors:
#curve(dinvgamma(x, 2.5-1,(2.5)*29.1^2), 0, 2000)
#abline(v=29.1^2)

#curve(dinvgamma(x,2.5-1,(2.5)*4.1^2), 0, 50)
#abline(v=4.2^2)

# Initial values
prec.y <- 1/29.1^2
prec.u_c <- 1
prec.x <- 1/4.2^2


## ----r------------------------------------------------------------------------
Y <- matrix(NA, 3*n, 3)


Y[1:n, 1] <- chl             # Regression model of interest response
Y[n+(1:n), 2] <- bmi         # Error model response
Y[2*n+(1:n), 3] <- rep(0, n) # Exposure model response

beta.0 <- c(rep(1, n), rep(NA, n), rep(NA, n))
beta.bmi <- c(1:n, rep(NA, n), rep(NA, n))
beta.age2 <- c(age2, rep(NA, n), rep(NA, n))
beta.age3 <- c(age3, rep(NA, n), rep(NA, n))

id.x <- c(rep(NA, n), 1:n, 1:n) 
weight.x <- c(rep(1, n), rep(1, n), rep(-1, n))

alpha.0 <- c(rep(NA, n), rep(NA, n), rep(1, n))
alpha.age2 <- c(rep(NA, n), rep(NA, n), age2)
alpha.age3 <- c(rep(NA, n), rep(NA, n), age3)

dd <- data.frame(Y = Y, 
                 beta.0 = beta.0,
                 beta.bmi = beta.bmi,
                 beta.age2 = beta.age2,
                 beta.age3 = beta.age3,
                 id.x = id.x,
                 weight.x = weight.x,
                 alpha.0 = alpha.0,
                 alpha.age2 = alpha.age2,
                 alpha.age3 = alpha.age3)


## ----r------------------------------------------------------------------------
formula = Y ~ - 1 + beta.0 + beta.age2 + beta.age3 + 
  f(beta.bmi, copy="id.x", 
    hyper = list(beta = list(param = prior.beta, fixed=FALSE))) +
  f(id.x, weight.x, model="iid", values = 1:n, 
    hyper = list(prec = list(initial = -15, fixed=TRUE))) +
  alpha.0 + alpha.age2 + alpha.age3


## ----r------------------------------------------------------------------------
Scale <- c(rep(1, n), rep(10^12, n), rep(1, n))


## ----r------------------------------------------------------------------------
model_missing <- inla(formula, data = dd, scale = Scale,
                     family = c("gaussian", "gaussian", "gaussian"),
                     control.family = list(
                       list(hyper = list(prec = list(initial = log(prec.y), 
                                                     param = prior.prec.y, 
                                                     fixed = FALSE))),
                       list(hyper = list(prec = list(initial = log(prec.u_c), 
                                                     param = prior.prec.u_c, 
                                                     fixed = TRUE))),
                       list(hyper = list(prec = list(initial = log(prec.x), 
                                                     param = prior.prec.x, 
                                                     fixed = FALSE)))
                     ),
                     control.fixed = list(
                       mean = list(beta.0 = prior.beta[1], 
                                   beta.age2 = prior.beta[1], 
                                   beta.age3 = prior.beta[1],  
                                   alpha.0 = prior.alpha[1], 
                                   alpha.age2 = prior.alpha[1],
                                   alpha.age3 = prior.alpha[1]), 
                       prec = list(beta.0 = prior.beta[2], 
                                   beta.age2 = prior.beta[2], 
                                   beta.age3 = prior.beta[2],  
                                   alpha.0 = prior.alpha[2], 
                                   alpha.age2 = prior.alpha[2],
                                   alpha.age3 = prior.alpha[2])),
                     verbose=F)












## ----r------------------------------------------------------------------------
# https://stefvanbuuren.name/RECAPworkshop/Practicals/RECAP_Practical_II.html 






## ----r------------------------------------------------------------------------
summary(model_missing)
