#######################################################################
# Calculate BF                                                        #
# ------------------------------------------------------------------- #
# Written by Ladislas Nalborczyk & Alexane Leclerc                    #
# E-mail: ladislas.nalborczyk@gmail.com                               #
# Last updated on May 23, 2025                                        #
#######################################################################

#######################
# healthy participants
#######################

# load packages
library(marginaleffects)
library(easystats)
library(tidyverse)
library(MetBrewer)
library(emmeans)
library(brms)

# loading model
fit_wiener <- readRDS(file = "models/ddm_HC_23-05-25.rds")

### 1: analyses for block x prime_emotion

# retrieving predictions at the individual level for drift rate
cond_effects_bias <- conditional_effects(
  # defining the model
  x = fit_wiener,
  # defining the interaction of interest
  effects = "prime_emotion:block",
  # defining the "conditions" (participant) on which computing predictions
  conditions = make_conditions(x = fit_wiener$data, vars = "participant"),
  # getting predictions about p (i.e., between 0 and 1)
  method = "posterior_epred",
  # including all random/varying effects
  re_formula = NULL,
  # needed to distinguish between positive and negative RTs
  negative_rt = TRUE
)[[1]]

# retrieving predictions at the individual level for bias
cond_effects_bias <- conditional_effects(
  # defining the model
  x = fit_wiener,
  # defining the interaction of interest
  effects = "prime_emotion:block",
  # defining the "conditions" (participant) on which computing predictions
  conditions = make_conditions(x = fit_wiener$data, vars = "participant"),
  # getting predictions about p (i.e., between 0 and 1)
  method = "posterior_epred",
  # including all random/varying effects
  re_formula = NULL,
  # needed to distinguish between positive and negative RTs
  negative_rt = TRUE,
  # defining the parameter we are interested in
  dpar = "bias"
)[[1]]


cond_effects_bias <- cond_effects[!cond_effects$prime_emotion == 'Morph',]
cond_effects_drift <- cond_effects_drift[!cond_effects_drift$prime_emotion == 'Morph',]


##### parameter estimates ####
# retrieving per-condition estimates for the drift rate
drift_samples_per_condition <- fit_wiener %>%
  # retrieving drift rate values per condition
  emmeans(~prime_emotion * block * participant, type = "response") %>%
  # retrieving posterior sample for each cell
  # gather_emmeans_draws()
  summary() %>%
  mutate(parameter = "drift") %>%
  rename(response = emmean)

# retrieving per-condition estimates for the boundary separation
bs_samples_per_condition <- fit_wiener %>%
  # retrieving drift rate values per condition
  emmeans(~prime_emotion * block, dpar = "bs", type = "response") %>%
  # retrieving posterior sample for each cell
  # gather_emmeans_draws()
  summary() %>%
  mutate(parameter = "bs")

# retrieving per-condition estimates for the non-decision time
ndt_samples_per_condition <- fit_wiener %>%
  # retrieving drift rate values per condition
  emmeans(~prime_emotion * block, dpar = "ndt", type = "response") %>%
  # retrieving posterior sample for each cell
  # gather_emmeans_draws()
  summary() %>%
  mutate(parameter = "ndt")

# retrieving per-condition estimates for the starting point
bias_samples_per_condition <- fit_wiener %>%
  # retrieving drift rate values per condition
  emmeans(~prime_emotion * block, dpar = "bias", type = "response") %>%
  # retrieving posterior sample for each cell
  # gather_emmeans_draws()
  summary() %>%
  mutate(parameter = "bias")

# binding everything together
parameter_estimates <- bind_rows(
  drift_samples_per_condition, bs_samples_per_condition,
  ndt_samples_per_condition, bias_samples_per_condition
)

# plotting it
parameter_estimates %>%
  ggplot(aes(x = prime_emotion, y = response, colour = block) ) +
  geom_errorbar(
    aes(ymin = lower.HPD, ymax = upper.HPD),
    width = 0,
    position = position_dodge(0.1)
  ) +
  geom_line(aes(group = block), position = position_dodge(0.1) ) +
  geom_point(position = position_dodge(0.1) ) +
  facet_wrap(~parameter, scales = "free") +
  labs(x = "Prime", y = "Parameter value") +
  theme_nice()



#####################
#### calculate BF ###
#####################

# specify priors for model
prior_model <- unupdate(fit_wiener)

### model for bias
prior_emmgrid <- emmeans(
  object = prior_model,
  specs = pairwise ~ prime_emotion | block,
  dpar = "bias"
)
posterior_emmgrid <- emmeans(
  object = fit_wiener,
  specs = pairwise ~ prime_emotion | block,
  dpar = "bias"
)
model_bfs_bias <- bayesfactor_parameters(posterior = posterior_emmgrid, prior = prior_emmgrid)

model_bfs_bias$BF <- exp(model_bfs_bias$log_BF)
model_bfs_bias$BF <- sprintf("%.4f", model_bfs_bias$BF)
model_bfs_bias <- model_bfs_bias %>% rename(BF_bias = BF)
model_bfs_bias <- model_bfs_bias %>% select(-log_BF)


posterior_summary <- summary(posterior_emmgrid, point.est = "mean")

# posterior_samples <- as.data.frame(posterior_emmgrid)  # extract posterior samples
# posterior_means <- colMeans(posterior_samples)  # compute the mean for each contrast


### model for drift rate
prior_emmgrid <- emmeans(
  object = prior_model,
  specs = pairwise ~ prime_emotion | block
)
posterior_emmgrid <- emmeans(
  object = fit_wiener,
  specs = pairwise ~ prime_emotion | block
)
model_bfs_drift <- bayesfactor_parameters(posterior = posterior_emmgrid, prior = prior_emmgrid)

model_bfs_drift$BF <- exp(model_bfs_drift$log_BF)
model_bfs_drift$BF <- sprintf("%.4f", model_bfs_drift$BF)
model_bfs_drift <- model_bfs_drift %>% rename(BF_drift = BF)
model_bfs_drift <- model_bfs_drift %>% select(-log_BF)

bf_table <- merge(model_bfs_drift, model_bfs_bias, by='Parameter')



### 2: analyses for prime_emotion x mood

### control only - both blocks
prior_trends <- emtrends(
  object = prior_model,
  specs = ~ prime_emotion,
  var = "mood",
)
posterior_trends <- emtrends(
  object = fit_wiener,
  specs = ~ prime_emotion,
  var = "mood"
)
model_mood_drift <- bayesfactor_parameters(posterior = posterior_trends, prior = prior_trends)
print(model_mood_drift)

model_mood_drift$BF <- exp(model_mood_drift$log_BF)
model_mood_drift$BF <- sprintf("%.4f", model_mood_drift$BF)
model_mood_drift <- model_mood_drift %>% rename(BF_drift = BF)
model_mood_drift <- model_mood_drift %>% select(-log_BF)

posterior_summary <- summary(posterior_trends, point.est = "mean")

# model mood bias
prior_trends <- emtrends(
  object = prior_model,
  specs = ~ prime_emotion,
  var = "mood",
  dpar='bias'
)
posterior_trends <- emtrends(
  object = fit_wiener,
  specs = ~ prime_emotion,
  var = "mood",
  dpar='bias'
)
model_mood_bias <- bayesfactor_parameters(posterior = posterior_trends, prior = prior_trends)
print(model_mood_bias)


model_mood_bias$BF <- exp(model_mood_bias$log_BF)
model_mood_bias$BF <- sprintf("%.4f", model_mood_bias$BF)
model_mood_bias <- model_mood_bias %>% rename(BF_bias = BF)
model_mood_bias <- model_mood_bias %>% select(-log_BF)

bf_table_mood_control <- merge(model_mood_drift, model_mood_bias, by='Parameter')

posterior_summary <- summary(posterior_trends, point.est = "mean")








# 
# 
# 
# #### retrieving predictions ####
# # retrieving predictions (about RTs) for prime_emotion
# conditional_effects(fit_wiener)$prime_emotion
# 
# # retrieving predictions for the drift rate
# conditional_effects(
#   x = fit_wiener,
#   # specifying effects of interest
#   effects = "mood",
#   negative_rt = TRUE,
#   method = "posterior_linpred"
# )
# 
# # retrieving predictions for the boundary separation
# conditional_effects(
#   x = fit_wiener,
#   # specifying effects of interest
#   # effects = "movie",
#   negative_rt = TRUE,
#   dpar = "bs",
#   method = "posterior_linpred"
# )
# 
# # retrieving predictions for the non-decision time
# conditional_effects(
#   x = fit_wiener,
#   # specifying effects of interest
#   # effects = "movie",
#   negative_rt = TRUE,
#   dpar = "ndt",
#   method = "posterior_linpred"
# )
# 
# # retrieving predictions for the bias (starting position)
# conditional_effects(
#   x = fit_wiener,
#   # specifying effects of interest
#   #effects = "mood",
#   negative_rt = TRUE,
#   dpar = "bias",
#   method = "posterior_linpred"
# )
# 
