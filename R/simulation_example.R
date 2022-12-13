## ----r------------------------------------------------------------------------
library(INLA)


## ----r------------------------------------------------------------------------
inla.setOption(num.threads = "1:1")


## ----r------------------------------------------------------------------------
set.seed(2022)
n <- 1000

# Covariate without error:
z <- rnorm(n, mean = 0, sd = 1)

# Berkson error:
r <- rnorm(n, mean = 1 + 2*z, sd = 1)
u_b <- rnorm(n, sd = 1)
x <- r + u_b

# Response:
y <- 1 + 2*x + 2*z + rnorm(n)

# Classical error:
u_c <- rnorm(n, sd = 1)
w <- r + u_c 

# Missingness:
m_pred <- -1.5 - 0.5*z # This gives a mean probability of missing of ca 0.2.
m_prob <- exp(m_pred)/(1 + exp(m_pred))

m_index <- rbinom(n, 1, prob = m_prob) # MAR
# m_index <- sample(1:n, 0.2*n, replace = FALSE) # MCAR
w[m_index] <- NA

simulated_data <- data.frame(y = y, w = w, z = z)


## ----r------------------------------------------------------------------------
attach(simulated_data)
n <- nrow(simulated_data)


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
  
# Initial values
prec.y <- 1
prec.u_b <- 1
prec.u_c <- 1
prec.r <- 1


## ----r------------------------------------------------------------------------
Y <- matrix(NA, 4*n, 4)

Y[1:n, 1] <- y                   # Regression model of interest response
Y[n+(1:n), 2] <- rep(0, n)       # Berkson error model response
Y[2*n+(1:n), 3] <- w             # Classical error model response
Y[3*n+(1:n), 4] <- rep(0, n)     # Imputation model response

beta.0 <- c(rep(1, n), rep(NA, 3*n))
beta.x <- c(1:n, rep(NA, 3*n))
beta.z <- c(z, rep(NA, 3*n))

id.x <- c(rep(NA, n), 1:n, rep(NA, n), rep(NA, n))
weight.x <- c(rep(NA, n), rep(-1, n), rep(NA, n), rep(NA, n))

id.r <- c(rep(NA, n), 1:n, 1:n, 1:n)
weight.r <- c(rep(NA, n), rep(1, n), rep(1, n), rep(-1, n))

alpha.0 = c(rep(NA, 3*n), rep(1, n))
alpha.z = c(rep(NA, 3*n), z)


## ----r------------------------------------------------------------------------
dd <- list(Y = Y,
           beta.0 = beta.0,
           beta.x = beta.x,
           beta.z = beta.z,
           id.x = id.x, 
           weight.x = weight.x,
           id.r = id.r,
           weight.r = weight.r,
           alpha.0 = alpha.0,
           alpha.z = alpha.z)


## ----r------------------------------------------------------------------------
formula = Y ~ - 1 + beta.0 + beta.z +
  f(beta.x, copy = "id.x",  
    hyper = list(beta = list(param = prior.beta, fixed = FALSE))) +
  f(id.x, weight.x, model = "iid", values = 1:n, 
    hyper = list(prec = list(initial = -15, fixed = TRUE))) +
  f(id.r, weight.r, model="iid", values = 1:n, 
    hyper = list(prec = list(initial = -15, fixed = TRUE))) + 
  alpha.0 + alpha.z


## ----r------------------------------------------------------------------------
model_sim <- inla(formula, data = dd, scale = scale.vec,
                  family = c("gaussian", "gaussian", "gaussian", "gaussian"),
                  control.family = list(
                    list(hyper = list(prec = list(initial = log(prec.y),
                                                  param = prior.prec.y,
                                                  fixed = FALSE))),
                    list(hyper = list(prec = list(initial = log(prec.u_b),
                                                  param = prior.prec.u_b,
                                                  fixed = FALSE))),
                    list(hyper = list(prec = list(initial = log(prec.u_c),
                                                  param = prior.prec.u_c,
                                                  fixed = FALSE))),
                    list(hyper = list(prec = list(initial = log(prec.r),
                                                  param = prior.prec.r,
                                                  fixed = FALSE)))
                  ),
                  control.predictor = list(compute = TRUE), 
                  control.fixed = list(
                    mean = list(beta.0 = prior.beta[1],
                                beta.z = prior.beta[1],
                                alpha.0 = prior.alpha[1],
                                alpha.z = prior.alpha[1]),
                    prec = list(beta.0 = prior.beta[2],
                                beta.z = prior.beta[2],
                                alpha.0 = prior.alpha[2],
                                alpha.z = prior.alpha[2]))
               )


## ----r------------------------------------------------------------------------
# Summary of fixed effects:
fixed <- model_sim$summary.fixed[1:5]
fixed[c("mean", "0.025quant", "0.975quant")]

# Summary of random effects:
hyper <- model_sim$summary.hyperpar[1:5]
hyper[c("mean", "0.025quant", "0.975quant")]


## ----r------------------------------------------------------------------------
library(tidyverse)
library(showtext)

true_values <- tibble::tribble(
  ~coef_name, ~mean,
  "Beta for beta.x", 2, 
  "beta.z", 2, 
  "beta.0", 1, 
  "alpha.0", 1, 
  "alpha.z", 2
)

# Prepare data for plotting
post_estimates <- bind_rows(fixed, hyper) %>% 
  janitor::clean_names() %>% 
  mutate(coef_name = rownames(.)) %>% 
  bind_rows(me_model = ., "truth" = true_values, .id = "model") %>% 
  filter(coef_name %in% c("beta.0", "beta.z", "Beta for beta.x", "alpha.0", "alpha.z")) %>% 
  mutate(coef_name = recode_factor(coef_name, "Beta for beta.x" = "beta.x")) %>%
  separate(coef_name, c("sub_model", "coef_name")) %>% 
  mutate(coef_pretty = paste0(sub_model, "[", coef_name, "]")) %>% 
  mutate(model = recode(model, "me_model" = "ME model", "truth" = "True value")) %>% 
  mutate(coef_pretty = fct_relevel(coef_pretty, levels = c("beta[0]", "beta[z]", "beta[x]", "alpha[0]", "alpha[z]"))) %>% 
  mutate(limit_lower = ifelse(coef_name==0, 0.5, 1.5),
         limit_upper = ifelse(coef_name==0, 1.2, 2.2))


# Colors
col_bgr <- "white" #"#fbf9f4"
col_text <- "#191919"
color_pal <- c("#BB5566", "#DDAA33", "#004488") # red, yellow, blue

# Loading fonts
f1 <- "Open Sans"
f2 <- "Open Sans"
font_add_google(name = f1, family = f1)
font_add_google(name = f2, family = f2)

# Plot theme
theme_model_summary <- theme_minimal(base_size = 18, base_family = "Open Sans") + 
  theme(
  axis.title.y = element_blank(),
  axis.text = element_text(size = 10, color = col_text),
  #axis.text.x = element_blank(),
  axis.ticks = element_blank(),
  legend.title = element_blank(),
  legend.text = element_text(size = 10),
  panel.background = element_rect(fill = col_bgr, color = col_bgr),
  #plot.background = element_rect(fill = col_bgr, color = "grey75", size = 1),
  legend.position = "none",
  strip.placement = "outside",
  strip.text = element_text(color = col_text),
  panel.grid.major.y = element_blank(),
  panel.grid.minor = element_blank(),
  plot.title.position = "plot",
  axis.line.x = element_line(size = 1, color = "grey65"),
  plot.margin = margin(rep(15, 4))
)

showtext.auto()

ggplot(post_estimates, aes(y = model)) +
  geom_linerange(aes(xmin = limit_lower, xmax = limit_upper), alpha = 0) +
  geom_linerange(aes(xmin = x0_025quant, xmax = x0_975quant, color = model), size = 1) +
  geom_point(aes(x = mean, color = model), size = 3) +
  geom_text(aes(x = mean, y = model, label = round(mean, 2)), 
            vjust = -1, family = "Open Sans") +
  scale_color_manual(values = color_pal[2:3], 
                     #labels = c("ME model", "True value"),
    guide = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  facet_wrap(vars(coef_pretty), nrow = 2, #switch = "x",
             labeller = label_parsed, scales = "free_x"
             ) +
  labs(x = "Posterior mean") +
  coord_cartesian(clip = "off") +
  theme_model_summary


## ----r------------------------------------------------------------------------
ggsave("figures/simulation_ex_figure.png", width = 7, height = 5)


## ----r------------------------------------------------------------------------
# Save the INLA-model so it can be summarized in the paper.
saveRDS(model_sim, file = "results/model_simulation.rds")
