## --------------------------------------------------------------------------------------------


## --------------------------------------------------------------------------------------------
library(survival)
library(flexsurv)
library(boot)
library(MASS)
library(rjags)
library(INLA)


## --------------------------------------------------------------------------------------------
# Load data   ==================================================================
mydata <- read.csv("data/bloodpressure.csv")
head(mydata)


## --------------------------------------------------------------------------------------------
# Create complete-case data set  ===============================================

# For the naive analysis, only the observations that have both sbp1 measurement
# and smoking status registered are used.
ccaData <- subset(mydata, (is.na(sbp1)==FALSE) & (is.na(smoke)==FALSE))
n <- dim(ccaData)[1]

head(ccaData)
# sbp1 sbp2 sex  age smoke diabetes d        t
# 5   0.5   NA   1 -0.7     1        0 0 12.00000
# 6  -0.1   NA   1 -1.0     1        0 0 21.16666
# 10  1.9   NA   1  2.0     0        0 0  3.25000
# 13  1.1  0.3   1  0.2     0        0 0 13.66666
# 15  0.8   NA   1  1.4     0        0 0  2.75000
# 16 -0.6   NA   1 -0.1     0        0 0  2.50000


# Function to view relevant INLA coefficient estimates =========================
view_relevant <- function(INLA_res, model_name){
  summ <- summary(INLA_res)
  fixed <- summ$fixed[2:5,1:6]
  beta.x <- summ$hyperpar[nrow(summ$hyperpar), ]
  cat(model_name, "\n")
  rbind(beta.x, fixed)
}


## --------------------------------------------------------------------------------------------
# Note: JAGS and INLA use the same Weibull parameterization
# (provided "variant" = 0 in INLA)

formula.naive <- inla.surv(t/10, d) ~ sbp1 + sex + age + smoke + diabetes

model_naive <- inla(formula.naive,
                     family ="weibullsurv",
                     data=ccaData,
                     control.family = list(list(variant = 0)))

summary(model_naive)


## ---- echo = FALSE, eval = FALSE-------------------------------------------------------------
## # Save for package:
## usethis::use_data(model_naive, overwrite = TRUE)


## --------------------------------------------------------------------------------------------
# Estimating sigma_u
W <- ccaData[,1:2]
W_mean <- rowMeans(ccaData[,1:2], na.rm = TRUE)
sigma_uu <- 0
for(i in 1:n){
  k <- sum(!is.na(W[i,])) # count number of repeats (1 or 2)
  for(j in 1:k){
    sigma_uu <- sigma_uu + (W[i,j] - W_mean[i])^2
  }
}
sigma_uu <- sigma_uu/sum(is.na(W$sbp2))

# Alternative calculation
sigma_uu_alt <- (sum((ccaData$sbp1-W_mean)^2) +
                   sum((ccaData$sbp2-W_mean)^2, na.rm = TRUE)) / sum(is.na(ccaData$sbp2))



## --------------------------------------------------------------------------------------------
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
prior.exp <- 0.01 # Exp(0.001) (INLA sets prior on theta, r~Exp(0.1*theta))
exp.init <- 1.4


## --------------------------------------------------------------------------------------------
# Specifying Y object
surv.time <- c(ccaData$t, rep(NA, 2*n))
event <- c(ccaData$d, rep(NA, 2*n))
Y.surv <- inla.surv(surv.time/10, event) # Divide by 10 because of numerical issues otherwise.
Y.expos.sbp <- c(rep(NA, n), rep(0, n), rep(NA, n))
Y.err.sbp <- c(rep(NA, 2*n), ccaData$sbp1) # Use only first measurement from complete case data
Y <- list(Y.surv, Y.expos.sbp, Y.err.sbp)

beta.0 <- c(rep(1, n), rep(NA, 2*n))
beta.sbp <- c(1:n, rep(NA, 2*n))
beta.sex <- c(ccaData$sex, rep(NA, 2*n))
beta.age <- c(ccaData$age, rep(NA, 2*n))
beta.smoke <- c(ccaData$smoke, rep(NA, 2*n))
beta.diabetes <- c(ccaData$diabetes, rep(NA, 2*n))

beta.sbp.copy <- c(rep(NA, n), 1:n, 1:n)
weight.sbp <- c(rep(NA, n), rep(-1, n), rep(1, n))

alpha.0 <- c(rep(NA, n), rep(1, n), rep(NA, n))
alpha.sex <- c(rep(NA, n), ccaData$sex, rep(NA, n))
alpha.age <- c(rep(NA, n), ccaData$age, rep(NA, n))
alpha.smoke <- c(rep(NA, n), ccaData$smoke, rep(NA, n))
alpha.diabetes <- c(rep(NA, n), ccaData$diabetes, rep(NA, n))

dd <- list(Y = Y,
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
           alpha.diabetes = alpha.diabetes)

# INLA formula with copy option
formula = Y ~ beta.0 - 1 +
  f(beta.sbp.copy, weight.sbp, model="iid", values = 1:n,
    hyper = list(prec = list(initial = -15, fixed=TRUE))) +
  f(beta.sbp, copy="beta.sbp.copy",
    hyper = list(beta = list(param = prior.beta, fixed=FALSE))) +
  beta.sex + beta.age + beta.smoke + beta.diabetes +
  alpha.0 + alpha.sex + alpha.age + alpha.smoke + alpha.diabetes

model_bloodpressure0 <- inla(formula, data = dd,
                 family = c("weibullsurv", "gaussian", "gaussian"),
                 control.family = list(
                   list(hyper = list(alpha = list(param = prior.exp,
                                                  initial = log(exp.init),
                                                  fixed = FALSE))),
                   list(hyper = list(prec = list(initial = log(prec.x),
                                                 param = prior.prec.x,
                                                 fixed = FALSE))),
                   list(hyper = list(prec = list(initial = log(prec.u),
                                                 param = prior.prec.u,
                                                 fixed = TRUE)))
                 ),
                 control.predictor = list(compute = TRUE), #compute pred. dist. of the missing obs in the response
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
                               alpha.diabetes = prior.alpha[2])),
                 verbose=F)
summary(model_bloodpressure0)
view_relevant(model_bloodpressure0, model_name = "Single measurement")


## ---- echo = FALSE, eval = FALSE-------------------------------------------------------------
## # Save for package:
## usethis::use_data(model_bloodpressure0, overwrite = TRUE)


## --------------------------------------------------------------------------------------------
# Specifying Y object
surv.time <- c(ccaData$t, rep(NA, 3*n))
event <- c(ccaData$d, rep(NA, 3*n))
Y.surv <- inla.surv(surv.time/10, event)
Y.expos.sbp <- c(rep(NA, n), rep(0, n), rep(NA, 2*n))
Y.err.sbp <- c(rep(NA, 2*n), ccaData$sbp1, ccaData$sbp2) # Use both measurements from complete case data
Y <- list(Y.surv, Y.expos.sbp, Y.err.sbp)

beta.0 <- c(rep(1, n), rep(NA, 3*n))
beta.sbp <- c(1:n, rep(NA, 3*n))
beta.sex <- c(ccaData$sex, rep(NA, 3*n))
beta.age <- c(ccaData$age, rep(NA, 3*n))
beta.smoke <- c(ccaData$smoke, rep(NA, 3*n))
beta.diabetes <- c(ccaData$diabetes, rep(NA, 3*n))

# Insert NAs in last model where w is NA
tt <- 1:n
tt[is.na(ccaData$sbp2)] <- NA
beta.sbp.copy <- c(rep(NA, n), 1:n, 1:n, tt)
weight.sbp <- c(rep(NA, n), rep(-1, n), rep(1, 2*n))

alpha.0 <- c(rep(NA, n), rep(1, n), rep(NA, 2*n))
alpha.sex <- c(rep(NA, n), ccaData$sex, rep(NA, 2*n))
alpha.age <- c(rep(NA, n), ccaData$age, rep(NA, 2*n))
alpha.smoke <- c(rep(NA, n), ccaData$smoke, rep(NA, 2*n))
alpha.diabetes <- c(rep(NA, n), ccaData$diabetes, rep(NA, 2*n))

#Scale <- c(rep(1, 4*n) )

dd <- list(Y = Y,
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
           alpha.diabetes = alpha.diabetes#,
           #Scale = Scale
           )

# INLA formula with copy option
formula = Y ~ beta.0 - 1 +
  f(beta.sbp.copy, weight.sbp, model="iid", values = 1:n,
    hyper = list(prec = list(initial = -15, fixed=TRUE))) +
  f(beta.sbp, copy="beta.sbp.copy",
    hyper = list(beta = list(param = prior.beta, fixed=FALSE))) +
  beta.sex + beta.age + beta.smoke + beta.diabetes +
  alpha.0 + alpha.sex + alpha.age + alpha.smoke + alpha.diabetes

model_bloodpressure1 <- inla(formula, data = dd,
                 family = c("weibullsurv", "gaussian", "gaussian"),
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
                 #scale = Scale,
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
                               alpha.diabetes = prior.alpha[2])),
                 verbose=F)
summary(model_bloodpressure1)
#saveRDS(inla.res, "PaperA_ME_and_missing_data/data/KB_inla_res1b.rds")

view_relevant(model_bloodpressure1, "Repeated measurement")


## ---- echo = FALSE, eval = FALSE-------------------------------------------------------------
## # Save for package:
## usethis::use_data(model_bloodpressure1, overwrite = TRUE)


## --------------------------------------------------------------------------------------------
data2 <- subset(mydata, is.na(smoke) == FALSE)
n <- nrow(data2)

# Specifying Y object
surv.time <- c(data2$t, rep(NA, 3*n))
event <- c(data2$d, rep(NA, 3*n))
Y.surv <- inla.surv(surv.time/10, event)
Y.expos.sbp <- c(rep(NA, n), rep(0, n), rep(NA, 2*n))
Y.err.sbp <- c(rep(NA, 2*n), data2$sbp1, data2$sbp2) # Use all available data
Y <- list(Y.surv, Y.expos.sbp, Y.err.sbp)

beta.0 <- c(rep(1, n), rep(NA, 3*n))
beta.sbp <- c(1:n, rep(NA, 3*n))
beta.sex <- c(data2$sex, rep(NA, 3*n))
beta.age <- c(data2$age, rep(NA, 3*n))
beta.smoke <- c(data2$smoke, rep(NA, 3*n))
beta.diabetes <- c(data2$diabetes, rep(NA, 3*n))

# Insert NAs in last model where w is NA
tt <- 1:n
tt[is.na(data2$sbp2)] <- NA
beta.sbp.copy <- c(rep(NA, n), 1:n, 1:n, tt)
weight.sbp <- c(rep(NA, n), rep(-1, n), rep(1, 2*n))

alpha.0 <- c(rep(NA, n), rep(1, n), rep(NA, 2*n))
alpha.sex <- c(rep(NA, n), data2$sex, rep(NA, 2*n))
alpha.age <- c(rep(NA, n), data2$age, rep(NA, 2*n))
alpha.smoke <- c(rep(NA, n), data2$smoke, rep(NA, 2*n))
alpha.diabetes <- c(rep(NA, n), data2$diabetes, rep(NA, 2*n))

Scale <- c(rep(1, 4*n) )

dd <- list(Y = Y,
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

# INLA formula with copy option
formula = Y ~ beta.0 - 1 +
  f(beta.sbp.copy, weight.sbp, model="iid", values = 1:n,
    hyper = list(prec = list(initial = -15, fixed=TRUE))) +
  f(beta.sbp, copy="beta.sbp.copy",
    hyper = list(beta = list(param = prior.beta, fixed=FALSE))) +
  beta.sex + beta.age + beta.smoke + beta.diabetes +
  alpha.0 + alpha.sex + alpha.age + alpha.smoke + alpha.diabetes

model_bloodpressure <- inla(formula, data = dd,
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
                               alpha.diabetes = prior.alpha[2])),
                 verbose=F)
summary(model_bloodpressure)
#saveRDS(inla.res2, "PaperA_ME_and_missing_data/data/KB_inla_res2.rds")

view_relevant(model_bloodpressure, "Repeated measurement")


## ---- echo = FALSE, eval = FALSE-------------------------------------------------------------
## # Save for package:
## usethis::use_data(model_bloodpressure, overwrite = TRUE)


## ---- eval = FALSE---------------------------------------------------------------------------
## res.naive <- readRDS("PaperA_ME_and_missing_data/data/KB_inla_res0.rds")
## res.singlemes <- readRDS("PaperA_ME_and_missing_data/data/KB_inla_res1a.rds")
## res.repmeas <- readRDS("PaperA_ME_and_missing_data/data/KB_inla_res1b.rds")
## res.complbp <- readRDS("PaperA_ME_and_missing_data/data/KB_inla_res2.rds")


## ---- eval = FALSE---------------------------------------------------------------------------
## # Corresponds to Naive?
## view_relevant(res.naive, "Naive model")
## # No corresponding model in KB
## view_relevant(res.singlemes, "Single SBP measurement (only complete cases of SBP1 and smoking)")
## # Corresponds to Bayes^c?
## view_relevant(res.repmeas, "Repeated SBP measurement (only complete cases of SBP1 and smoking)")
## # Would correspond to Bayes^d if we had exposure model for smoking
## view_relevant(res.complbp, "Missingness in both SBP measurements (only complete cases of smoking)")

