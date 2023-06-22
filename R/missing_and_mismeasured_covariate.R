## ----r------------------------------------------------------------------------
library(survival)   # Survival modelling
library(INLA)       # Modelling
library(tidyverse)  # Data wrangling and visualisation


## ----r------------------------------------------------------------------------
inla.setOption(num.threads = "1:1")


## ----r------------------------------------------------------------------------
view_relevant <- function(INLA_res, model_name){
  fixed <- INLA_res$summary.fixed[c("mean", "0.025quant", "0.975quant")]
  hyper <- INLA_res$summary.hyperpar[c("mean", "0.025quant", "0.975quant")]
  relevant <- c("beta.age", "beta.diabetes", "beta.sex", "beta.smoke")
  
  cat(model_name, "\n")
  beta.x <- hyper[nrow(hyper), ]
  rownames(beta.x) <- "beta.sbp"
  rbind(beta.x, fixed[relevant, ])
}


## ----r------------------------------------------------------------------------
full_data <- read.csv("data/bloodpressure.csv")
head(full_data)


## ----r------------------------------------------------------------------------
data <- full_data %>% drop_na(smoke)
head(data)


## ----r------------------------------------------------------------------------
ccdata <- data %>% drop_na(sbp1)
n <- nrow(ccdata)
head(ccdata)


## ----r------------------------------------------------------------------------
# Note: JAGS and INLA use the same Weibull parameterization
# (provided "variant" = 0 in INLA)

formula.naive <- inla.surv(t/10, d) ~ sbp1 + sex + age + smoke + diabetes

model_naive <- inla(formula.naive,
                    family ="weibullsurv",
                    data = ccdata,
                    control.family = list(list(variant = 0)))

cat("Naive model")
model_naive$summary.fixed[c("mean", "0.025quant", "0.975quant")]


## ----r------------------------------------------------------------------------
# Priors for measurement error variance and true x-value
prior.prec.u <- c(0.5, 0.5) # Gamma(0.5, 0.5) (same as Keogh&Bartlett)
prior.prec.x <- c(0.5, 0.5) # Gamma(0.5, 0.5) (same as K&B)
#prec.u = 1/sigma_uu
prec.u <- 2.8
prec.x = 1

#curve(dgamma(x, shape = 0.5, rate = 0.5))

# Priors for exposure model coefficients
prior.alpha <- c(0, 1/10000) # N(0, 10^4) (same as K&B)

# Priors for model of interest coefficients
prior.beta = c(0, 1/1000) # This has a Gaussian prior
# (K&B specify prior -beta/r ~ N(0,10^6). Since r has prior Exp(0.001),
# the expected value of r is 1000, and so if we fix this,
# we can use prior beta ~ N(0, 1000).)

# Prior for shape parameter of the Weibull survival model
prior.exp <- 0.01 # Gamma(1, 0.001) ~ Exp(0.001) (INLA sets prior on theta, r~Exp(0.1*theta))
exp.init <- 1.4


## ----r------------------------------------------------------------------------
n <- nrow(ccdata)

# Specifying Y object
surv.time <- c(ccdata$t, rep(NA, 3*n))
event <- c(ccdata$d, rep(NA, 3*n))
Y.surv <- inla.surv(surv.time/10, event)
Y.expos.sbp <- c(rep(NA, n), rep(0, n), rep(NA, 2*n))
Y.err.sbp <- c(rep(NA, 2*n), ccdata$sbp1, ccdata$sbp2) # Use all available data
Y <- list(Y.surv, Y.expos.sbp, Y.err.sbp)

beta.0 <- c(rep(1, n), rep(NA, 3*n))
beta.sbp <- c(1:n, rep(NA, 3*n))
beta.sex <- c(ccdata$sex, rep(NA, 3*n))
beta.age <- c(ccdata$age, rep(NA, 3*n))
beta.smoke <- c(ccdata$smoke, rep(NA, 3*n))
beta.diabetes <- c(ccdata$diabetes, rep(NA, 3*n))

# Insert NAs in last model where w is NA
tt <- 1:n
tt[is.na(ccdata$sbp2)] <- NA
beta.sbp.copy <- c(rep(NA, n), 1:n, 1:n, tt)
weight.sbp <- c(rep(NA, n), rep(-1, n), rep(1, 2*n))

alpha.0 <- c(rep(NA, n), rep(1, n), rep(NA, 2*n))
alpha.sex <- c(rep(NA, n), ccdata$sex, rep(NA, 2*n))
alpha.age <- c(rep(NA, n), ccdata$age, rep(NA, 2*n))
alpha.smoke <- c(rep(NA, n), ccdata$smoke, rep(NA, 2*n))
alpha.diabetes <- c(rep(NA, n), ccdata$diabetes, rep(NA, 2*n))

Scale <- c(rep(1, 4*n) )

mat1 <- list(Y = Y,
           beta.0 = beta.0,
           beta.sbp = beta.sbp,
           beta.sex = beta.sex,
           beta.age = beta.age,
           beta.smoke = beta.smoke,
           beta.diabetes = beta.diabetes,
           beta.sbp.copy = beta.sbp.copy,
           weight.sbp = weight.sbp,
           alpha.0 = alpha.0,
           alpha.sex = alpha.sex,
           alpha.age = alpha.age,
           alpha.smoke = alpha.smoke,
           alpha.diabetes = alpha.diabetes,
           Scale = Scale)


## ----r------------------------------------------------------------------------
# INLA formula with copy option
formula1 <- Y ~ beta.0 - 1 +
  f(beta.sbp.copy, weight.sbp, model="iid", values = 1:n,
    hyper = list(prec = list(initial = -15, fixed=TRUE))) +
  f(beta.sbp, copy="beta.sbp.copy",
    hyper = list(beta = list(param = prior.beta, fixed=FALSE))) +
  beta.sex + beta.age + beta.smoke + beta.diabetes +
  alpha.0 + alpha.sex + alpha.age + alpha.smoke + alpha.diabetes


## ----r------------------------------------------------------------------------
model1 <- inla(formula1, data = mat1,
                 family = c("weibull.surv", "gaussian", "gaussian"),
                 control.family = list(
                   list(hyper = list(alpha = list(param = prior.exp,
                                                  initial = log(exp.init),
                                                  fixed = FALSE))),
                   list(hyper = list(prec = list(initial = log(prec.x),
                                                 param = prior.prec.x,
                                                 fixed = FALSE))),
                   list(hyper = list(prec = list(initial = log(prec.u),
                                                 param = prior.prec.u,
                                                 fixed = FALSE)))
                 ),
                 control.predictor=list(link=3),
                 scale = Scale,
                 control.fixed = list(
                   mean = list(beta.0 = prior.beta[1],
                               beta.sex = prior.beta[1],
                               beta.age = prior.beta[1],
                               beta.smoke = prior.beta[1],
                               beta.diabetes = prior.beta[1],
                               alpha.0 = prior.alpha[1],
                               alpha.sex = prior.alpha[1],
                               alpha.age = prior.alpha[1],
                               alpha.smoke = prior.alpha[1],
                               alpha.diabetes = prior.alpha[1]),
                   prec = list(beta.0 = prior.beta[2],
                               beta.sex = prior.beta[2],
                               beta.age = prior.beta[2],
                               beta.smoke = prior.beta[2],
                               beta.diabetes = prior.beta[2],
                               alpha.0 = prior.alpha[2],
                               alpha.sex = prior.alpha[2],
                               alpha.age = prior.alpha[2],
                               alpha.smoke = prior.alpha[2],
                               alpha.diabetes = prior.alpha[2])))

view_relevant(model1, "Repeated measurement")


## ----r------------------------------------------------------------------------
n <- nrow(data)

# Specifying Y object
surv.time <- c(data$t, rep(NA, 3*n))
event <- c(data$d, rep(NA, 3*n))
Y.surv <- inla.surv(surv.time/10, event)
Y.expos.sbp <- c(rep(NA, n), rep(0, n), rep(NA, 2*n))
Y.err.sbp <- c(rep(NA, 2*n), data$sbp1, data$sbp2) # Use all available data
Y <- list(Y.surv, Y.expos.sbp, Y.err.sbp)

beta.0 <- c(rep(1, n), rep(NA, 3*n))
beta.sbp <- c(1:n, rep(NA, 3*n))
beta.sex <- c(data$sex, rep(NA, 3*n))
beta.age <- c(data$age, rep(NA, 3*n))
beta.smoke <- c(data$smoke, rep(NA, 3*n))
beta.diabetes <- c(data$diabetes, rep(NA, 3*n))

# Insert NAs in last model where w is NA
tt <- 1:n
tt[is.na(data$sbp2)] <- NA
beta.sbp.copy <- c(rep(NA, n), 1:n, 1:n, tt)
weight.sbp <- c(rep(NA, n), rep(-1, n), rep(1, 2*n))

alpha.0 <- c(rep(NA, n), rep(1, n), rep(NA, 2*n))
alpha.sex <- c(rep(NA, n), data$sex, rep(NA, 2*n))
alpha.age <- c(rep(NA, n), data$age, rep(NA, 2*n))
alpha.smoke <- c(rep(NA, n), data$smoke, rep(NA, 2*n))
alpha.diabetes <- c(rep(NA, n), data$diabetes, rep(NA, 2*n))

Scale <- c(rep(1, 4*n) )

mat2 <- list(Y = Y,
           beta.0 = beta.0,
           beta.sbp = beta.sbp,
           beta.sex = beta.sex,
           beta.age = beta.age,
           beta.smoke = beta.smoke,
           beta.diabetes = beta.diabetes,
           beta.sbp.copy = beta.sbp.copy,
           weight.sbp = weight.sbp,
           alpha.0 = alpha.0,
           alpha.sex = alpha.sex,
           alpha.age = alpha.age,
           alpha.smoke = alpha.smoke,
           alpha.diabetes = alpha.diabetes,
           Scale = Scale)


## ----r------------------------------------------------------------------------
# INLA formula with copy option
formula2 <- Y ~ beta.0 - 1 +
  f(beta.sbp.copy, weight.sbp, model="iid", values = 1:n,
    hyper = list(prec = list(initial = -15, fixed=TRUE))) +
  f(beta.sbp, copy="beta.sbp.copy",
    hyper = list(beta = list(param = prior.beta, fixed=FALSE))) +
  beta.sex + beta.age + beta.smoke + beta.diabetes +
  alpha.0 + alpha.sex + alpha.age + alpha.smoke + alpha.diabetes


## ----r------------------------------------------------------------------------
model_bloodpressure <- inla(formula2, data = mat2,
                 family = c("weibull.surv", "gaussian", "gaussian"),
                 control.family = list(
                   list(hyper = list(alpha = list(param = prior.exp,
                                                  initial = log(exp.init),
                                                  fixed = FALSE))),
                   list(hyper = list(prec = list(initial = log(prec.x),
                                                 param = prior.prec.x,
                                                 fixed = FALSE))),
                   list(hyper = list(prec = list(initial = log(prec.u),
                                                 param = prior.prec.u,
                                                 fixed = FALSE)))
                 ),
                 control.predictor=list(link=3),
                 scale = Scale,
                 control.fixed = list(
                   mean = list(beta.0 = prior.beta[1],
                               beta.sex = prior.beta[1],
                               beta.age = prior.beta[1],
                               beta.smoke = prior.beta[1],
                               beta.diabetes = prior.beta[1],
                               alpha.0 = prior.alpha[1],
                               alpha.sex = prior.alpha[1],
                               alpha.age = prior.alpha[1],
                               alpha.smoke = prior.alpha[1],
                               alpha.diabetes = prior.alpha[1]),
                   prec = list(beta.0 = prior.beta[2],
                               beta.sex = prior.beta[2],
                               beta.age = prior.beta[2],
                               beta.smoke = prior.beta[2],
                               beta.diabetes = prior.beta[2],
                               alpha.0 = prior.alpha[2],
                               alpha.sex = prior.alpha[2],
                               alpha.age = prior.alpha[2],
                               alpha.smoke = prior.alpha[2],
                               alpha.diabetes = prior.alpha[2])))

view_relevant(model_bloodpressure, "Repeated measurement")


## ----r------------------------------------------------------------------------
# Corresponds to Naive?
cat("Naive model")
naive <- model_naive$summary.fixed[c("mean", "0.025quant", "0.975quant")]
naive

# Corresponds to Bayes^c
me_adjusted <- view_relevant(model1, "Adjusts for ME, but only complete cases for SBP1")
me_adjusted

# Would correspond to Bayes^d if we had exposure model for smoking
missing_adjusted <- view_relevant(model_bloodpressure, "Missingness in both SBP measurements (only complete cases of smoking)")
missing_adjusted


## ----r------------------------------------------------------------------------
naive <- naive %>% 
  mutate(coef_name = rownames(.)) %>% 
  filter(coef_name != "(Intercept)") %>% 
  mutate(coef_name = recode(coef_name, sbp1 = "sbp")) %>% 
  mutate(coef_name = paste0("beta.", coef_name))

me_adjusted <- me_adjusted %>% mutate(coef_name = rownames(.))

missing_adjusted <- missing_adjusted %>% mutate(coef_name = rownames(.))

kb_d <- tibble::tribble(
  ~"mean", ~"0.025quant", ~"0.975quant", ~"coef_name",
  0.121, 0.055, 0.186, "beta.sbp", 
  0.47, 0.36, 0.57, "beta.sex",
  1.02, 0.94, 1.09, "beta.age",
  0.25, 0.08, 0.42, "beta.smoke", 
  0.68, 0.55, 0.82, "beta.diabetes"
)

kb_c <- tibble::tribble(
  ~"mean", ~"0.025quant", ~"0.975quant", ~"coef_name",
  0.114, 0.015, 0.211, "beta.sbp", 
  0.49, 0.30, 0.68, "beta.sex",
  0.87, 0.75, 0.98, "beta.age",
  0.26, 0.07, 0.45, "beta.smoke", 
  0.50, 0.27, 0.72, "beta.diabetes"
)

all_models <- bind_rows("Complete case,\nnaive" = naive, 
                        "Complete case,\nME adjusted" = me_adjusted, 
                        "Imputation,\nME adjusted" = missing_adjusted, 
                        .id = "model") %>% 
  tidyr::separate(coef_name, c("sub_model", "coef_name")) %>% 
  mutate(coef_pretty = paste0("beta[", coef_name, "]")) %>% 
  rename(quant_0.025 = "0.025quant", quant_0.975 = "0.975quant") %>% 
  mutate(model = fct_relevel(model, "Imputation,\nME adjusted",
                                    "Complete case,\nME adjusted",
                                    "Complete case,\nnaive"))




## ----r------------------------------------------------------------------------
ggsave("figures/bloodpressure_figure.pdf", height = 4, width = 10, dpi = 600)
ggsave("figures/bloodpressure_figure.eps", height = 4, width = 10, dpi = 600, 
       device = cairo_ps)

