## ----r------------------------------------------------------------------------
library(INLA)
library(inlabru)



## ----r------------------------------------------------------------------------
inla.setOption(num.threads = "1:1")


## ----r------------------------------------------------------------------------
data <- read.csv("data/simulated_data.csv")
n <- nrow(data)


## ----r------------------------------------------------------------------------
# Priors for model of interest coefficients
prior.beta = c(0, 1/1000) # N(0, 10^3)

# Priors for exposure model coefficients
prior.alpha <- c(0, 1/10000) # N(0, 10^4) 

# Priors for y, measurement error and true x-value precision
prior.prec.y <- c(0.5, 0.5) # Gamma(0.5, 0.5)
prior.prec.u_b <- c(0.5, 0.5) # Gamma(0.5, 0.5)
prior.prec.u_c <- c(0.5, 0.5) # Gamma(0.5, 0.5)
prior.prec.r <- c(0.5, 0.5) # Gamma(0.5, 0.5) 
prior.prec.x <- c(0.5, 0.5)


# Initial values
prec.y <- 1
prec.u_b <- 1
prec.u_c <- 1
prec.r <- 1
prec.x <- 1


## ----r------------------------------------------------------------------------
# Regression model of interest 
data_moi <- data.frame(y = data$y, z = data$z, weight.x = 1, weight.r = 1, r = 1:n) 

# Berkson ME model
data_berkson <- data.frame(zero = 0, weight.x = -1, weight.r = 1, r = 1:n) 

# Classical ME model
data_classical <- data.frame(w = data$w, weight.x = 1, weight.r = 1, r = 1:n)

# Imputation model
data_imputation <- data.frame(zero = 0, z = data$z, weight.x = 1, weight.r = -1, r = 1:n)


## ----r------------------------------------------------------------------------
cmp_new <- ~ Intercept(1, model = "linear", prec.linear = prior.beta[2]) +
  beta_z(main = z, model = "linear", prec.linear = prior.beta[2]) +
  x_eff(r, weight.x, model = "iid",  hyper = list(prec = list(initial = -15, fixed=TRUE))) +
  x_eff_copy(r, copy="x_eff", hyper = list(beta = list(param = prior.beta, fixed=FALSE))) +
  r_eff(r, weight.r, model = "iid",  hyper = list(prec = list(initial = -15, fixed=TRUE))) +
  alpha_0(main = 1, model = "linear", prec.linear = prior.alpha[2]) +
  alpha_z(main = z, model = "linear", prec.linear = prior.alpha[2])


lik_moi <- like(formula = y ~ .,
                family = "gaussian",
                include = c("Intercept", "beta_z", "x_eff_copy"),
                control.family = list(
                  hyper = list(prec = list(initial = log(prec.y), 
                                           param = prior.prec.y, 
                                           fixed = FALSE))),
            data = data_moi)

lik_berkson <- like(formula = zero ~ .,
                    family = "gaussian",
                    include = c("x_eff", "r_eff"),
                    control.family = list(
                      hyper = list(prec = list(initial = log(prec.u_b),
                                               param = prior.prec.u_b,
                                               fixed = FALSE))),
                    data = data_berkson)

lik_classical <- like(formula = w ~ .,
                      family = "gaussian",
                      include = c("r_eff"),
                      control.family =  list(
                        hyper = list(prec = list(initial = log(prec.u_c), 
                                                 param = prior.prec.u_c, 
                                                 fixed = FALSE))),
            data  = data_classical)

lik_imputation <- like(formula = zero ~ .,
                       family = "gaussian",
                       include = c("alpha_0", "alpha_z","r_eff"),
                       control.family =   list(
                         hyper = list(prec = list(initial = log(prec.r), 
                                                  param = prior.prec.r, 
                                                  fixed = FALSE))),
             data = data_imputation)


## ----r------------------------------------------------------------------------
bru_options_set(bru_verbose = 1)
fit <- bru(components = cmp_new,
           lik_moi,
           lik_berkson,
           lik_classical,
           lik_imputation,
           options = list(verbose = F,
                          bru_max_iter = 20,
                          inla.mode  = "experimental"))
fit$summary.fixed
fit$summary.hyperpar

