## ----r------------------------------------------------------------------------
library(INLA)
library(inlabru)



## ----r------------------------------------------------------------------------
# Not sure if this is necessary with inlabru (TODO: check)
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
prior.prec.x <- c(0.5, 0.5) # Gamma(0.5, 0.5) 

# Initial values
prec.y <- 1
prec.u_b <- 1
prec.u_c <- 1
prec.x <- 1


## ----r------------------------------------------------------------------------
data1 = data.frame(y = data$y, z = data$z, weight = 1, r = 1:n)
data2 = data.frame(w = data$w, weight = 1, r = 1:n)
data3 = data.frame(zero = 0, z = data$z, weight = -1, r = 1:n)


cmp = ~ Intercept(1, model = "linear", prec.linear = prior.beta[2]) +
  beta_z(main = z, model = "linear", prec.linear = prior.beta[2]) +
  u_b(main = r, model = "iid",  hyper = list(prec = list(initial = log(1), fixed=TRUE))) +
  r_eff(r, weight, model = "iid",  hyper = list(prec = list(initial = -15, fixed=TRUE))) +
  r_eff_copy(r, copy="r_eff", 
             hyper = list(beta = list(param = prior.beta, fixed=FALSE))) +
  alpha_0(main = 1, model = "linear", prec.linear = prior.alpha[2]) +
  alpha_z(main = z, model = "linear", prec.linear = prior.alpha[2])


lik1 = like(formula = y ~ Intercept + beta_z + u_b + r_eff_copy,
            family = "gaussian",
            include = c("Intercept","beta_z","u_b","r_eff_copy"),
            control.family = list(hyper = list(prec = list(initial = log(prec.y), 
                                                           param = prior.prec.y, 
                                                           fixed = FALSE))),
            data = data1)

lik2 = like(formula = w ~ r_eff,
            family = "gaussian",
            include = c("r_eff"),
            control.family =  list(hyper = list(prec = list(initial = log(prec.u_c), 
                                                            param = prior.prec.u_c, 
                                                            fixed = TRUE))),
            data  = data2)

lik3  = like(formula = zero ~ alpha_0 + alpha_z + r_eff,
             family = "gaussian",
             include = c("alpha_0", "alpha_z","r_eff"),
             control.family =   list(hyper = list(prec = list(initial = log(prec.x), 
                                                              param = prior.prec.x, 
                                                              fixed = FALSE))),
             data = data3)
# Note: 
# formula = y ~ .,
# formula = w ~ .,
# formula = zero ~ .,
# does exactly the same in this case.


bru_options_set(bru_verbose = 1)
fit = bru(components = cmp,
          lik1,
          lik2,
          lik3,
          options = list(verbose = F,
                         bru_max_iter = 20,
                         inla.mode  = "experimental"))
fit$summary.fixed


## ----r------------------------------------------------------------------------
# no copy -----------------------------------------------------------------

cmp2 = ~ Intercept(1, model = "linear", prec.linear = prior.beta[2]) +
  beta_z(main = z, model = "linear", prec.linear = prior.beta[2]) +
  u_b(main = r, model = "iid",  hyper = list(prec = list(initial = log(1), fixed=TRUE))) +
  r_eff(r, weight, model = "iid",  hyper = list(prec = list(initial = -15, fixed=TRUE))) +
  #  r_eff_copy(r, copy="r_eff", 
  #             hyper = list(beta = list(param = prior.beta, fixed=FALSE))) +
  beta_u(main = 1, model = "linear", prec.linear = 0.001) +
  alpha_0(main = 1, model = "linear", prec.linear = prior.alpha[2]) +
  alpha_z(main = z, model = "linear", prec.linear = prior.alpha[2])


lik1 = like(formula = y ~ Intercept + beta_z + u_b + beta_u * r_eff ,
            family = "gaussian",
            data = data1,
            include = c("Intercept","beta_z","u_b","beta_u", "r_eff"),
            control.family = list(hyper = list(prec = list(initial = log(prec.y), 
                                                           param = prior.prec.y, 
                                                           fixed = FALSE))))

lik2 = like(formula = w ~ .,
            family = "gaussian",
            include = c("r_eff"),
            control.family =  list(hyper = list(prec = list(initial = log(prec.u_c), 
                                                            param = prior.prec.u_c, 
                                                            fixed = TRUE))),
            data  = data2)

lik3  = like(formula = zero ~ . ,
             data = data3,
             family = "gaussian",
             include = c("alpha_0", "alpha_z","r_eff"),
             control.family =   list(hyper = list(prec = list(initial = log(prec.x), 
                                                              param = prior.prec.x, 
                                                              fixed = FALSE))))



bru_options_set(bru_verbose = 1)
fit2 = bru(components = cmp2,
           lik1,
           lik2,
           lik3,
           options = list(verbose = F,
                          bru_max_iter = 20,
                          inla.mode  = "experimental"))
fit2 = bru_rerun(fit2)

fit2$summary.fixed

