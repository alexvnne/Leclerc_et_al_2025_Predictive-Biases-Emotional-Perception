#######################################################################
# Fitting the Bayesian multilevel DDM                                 #
# ------------------------------------------------------------------- #
# Written by Ladislas Nalborczyk                                      #
# E-mail: ladislas.nalborczyk@gmail.com                               #
# Last updated on May 23, 2025                                        #
#######################################################################

######################
# healthy participants
######################

# load packages

library(marginaleffects)
library(easystats)
library(tidyverse)
library(MetBrewer)
library(emmeans)
library(brms)

######################
# importing the data #
######################

# set working directory
setwd('/home/al274392/Desktop/DDM')

# fit a model to estimate the 4 DDM parameters for block * priming * mood on morphs
# mood comes from a common PCA with patients

df <- read.csv(file = "ddm_morph_HC_patients-Alexane.csv") %>%
  rename(rt = RT, response = Response, mood = PC1, 
         participant = ID, block = Block, prime_emotion = Prime) %>%
  filter(rt > 0.1 & rt < 3) %>%
  filter(Group == 'HC') %>%
  mutate(block = factor(block)) 
df$X <- NULL

# computing RT mean and summary per condition
df %>% summarise(mean = mean(rt), sd = sd(rt), .by = c(block, prime_emotion) )

# plotting distributions of RTs per condition
df %>%
  ggplot(aes(x = rt, colour = prime_emotion, fill = prime_emotion) ) +
  geom_density(alpha = 0.4) +
  facet_wrap(~block) +
  theme_bw(base_size = 12, base_family = "Open Sans") +
  labs(x = "Reaction time (s)", y = "Density")

# defining the contrasts
contrasts(df$block) <- c(-0.5, +0.5)


###########################
# defining the brms model #
###########################

# defining the model formula (one "linear model" per parameter)
formula <- brmsformula(
  # drift rate (delta)
  rt | dec(response) ~ 1 + prime_emotion * block * mood + (1 + prime_emotion * block | participant),
  # boundary separation parameter (alpha)
  bs ~ 1 + prime_emotion * block * mood + (1 + prime_emotion * block | participant),
  # non-decision time (tau)
  ndt ~ 1 + prime_emotion * block * mood + (1 + prime_emotion * block | participant),
  # starting point or bias (beta)
  bias ~ 1 + prime_emotion * block * mood + (1 + prime_emotion * block | participant)
)

# defining the priors
priors <- c(
  # priors for the intercepts
  prior(normal(0, 1), class = "Intercept"),
  prior(normal(0, 1), class = "Intercept", dpar = "bs"),
  prior(normal(0, 1), class = "Intercept", dpar = "ndt"),
  prior(normal(0, 1), class = "Intercept", dpar = "bias"),
  
  # priors for the slopes
  prior(normal(0, 1), class = "b"),
  prior(normal(0, 1), class = "b", dpar = "bs"),
  prior(normal(0, 1), class = "b", dpar = "ndt"),
  prior(normal(0, 1), class = "b", dpar = 'bias'),
  
  # priors on the SD of the varying effects
  prior(exponential(1), class = "sd"),
  prior(exponential(1), class = "sd", dpar = "bs"),
  prior(exponential(1), class = "sd", dpar = "ndt"),
  prior(exponential(1), class = "sd", dpar = "bias")
)

# specify initial values to help the model start sampling
# (with small variation between chains)
chains <- 32 # number of chains
#cores <- ?
epsilon <- 0.1 # variability in starting value for the NDT intercept
get_init_value <- function (x) list(Intercept_ndt = rnorm(n = 1, mean = x, sd = epsilon) )
inits_drift <- replicate(chains, get_init_value(-3), simplify = FALSE)

# fitting the model
fit_wiener <- brm(
  formula = formula,
  data = df,
  # specifying the family and link functions for each parameter
  family = wiener(
    link = "identity", link_bs = "log",
    link_ndt = "log", link_bias = "logit"
  ),
  # comment this line to use default priors
  prior = priors,
  # list of initialisation values
  init = inits_drift,
  init_r = 0.01,
  warmup = 1000, iter = 3000, 
  chains = chains, cores = cores,
  # control = list(adapt_delta = 0.99, max_treedepth = 15),
  control = list(adapt_delta = 0.99),
  # saves the model (as .rds) or loads it if it already exists
  file = "models/ddm_HC_23-05-25.rds",
  # needed for hypothesis testing
  sample_prior = TRUE,
  # speed optimisation
  # https://discourse.mc-stan.org/t/30-40-drop-in-sampling-time-using-stanc-o1-optimizations/28347
  # backend = "cmdstanr", stan_model_args = list(stanc_options = list("O1") )
)

# retrieving model estimates
summary(fit_wiener)

# plotting MCMC diagnostics
# combo can be hist, dens, dens_overlay, trace, trace_highlight...
# cf. https://mc-stan.org/bayesplot/reference/MCMC-overview.html
plot(
  x = fit_wiener, combo = c("dens_overlay", "trace"),
  variable = variables(fit_wiener)[1:4],
  ask = FALSE
)

# posterior predictive checking
pp_check(object = fit_wiener, ndraws = 10) +
  labs(
    title = "Posterior predictive check",
    subtitle = "Observed and simulated distributions of RTs",
    x = "Reaction time", y = "Density"
  )
