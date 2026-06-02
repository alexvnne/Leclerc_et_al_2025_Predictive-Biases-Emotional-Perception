#######################################################################
# Calculate BF                                                        #
# ------------------------------------------------------------------- #
# Written by Ladislas Nalborczyk & Alexane Leclerc                    #
# E-mail: ladislas.nalborczyk@gmail.com                               #
# Last updated on May 23, 2025                                        #
#######################################################################

###########################################
# patients (bipolar disorder & depression)
###########################################

# load packages
library(marginaleffects)
library(easystats)
library(tidyverse)
library(MetBrewer)
library(emmeans)
library(brms)

### 1: analyses for group x block x prime_emotion

# load model
fit_wiener <- readRDS(file = "models/ddm_morph_patients.rds") 

# retrieving predictions at the individual level for drift rate
conditions_df <- unique(fit_wiener$data[, c("participant", "group")])
cond_effects_drift <- conditional_effects(
  # defining the model
  x = fit_wiener,
  # defining the interaction of interest
  effects = "prime_emotion:block",  
  # defining the "conditions" (participant) on which computing predictions
  conditions = conditions_df, 
  # getting predictions about p (i.e., between 0 and 1)
  method = "posterior_epred",
  # including all random/varying effects
  re_formula = NULL,
  # needed to distinguish between positive and negative RTs
  negative_rt = TRUE,
)[[1]] 

# retrieving predictions at the individual level for bias
cond_effects_bias <- conditional_effects(
  # defining the model
  x = fit_wiener,
  # defining the interaction of interest
  effects = "prime_emotion:block",  
  # defining the "conditions" (participant) on which computing predictions
  conditions = conditions_df, 
  # getting predictions about p (i.e., between 0 and 1)
  method = "posterior_epred",
  # including all random/varying effects
  re_formula = NULL,
  # needed to distinguish between positive and negative RTs
  negative_rt = TRUE,
  # defining the parameter we are interested in
  dpar = 'bias'
)[[1]] 

# parameter estimates + plot
# retrieving per-condition estimates for the drift rate
drift_samples_per_condition <- fit_wiener %>%
  # retrieving drift rate values per condition
  emmeans(~prime_emotion * block , type = "response") %>%
  #emmeans(~prime_emotion * block * group, type = "response") %>%
  # retrieving posterior sample for each cell
  #gather_emmeans_draws()
  summary() %>%
  mutate(parameter = "drift") %>%
  rename(response = emmean)

# retrieving per-condition estimates for the boundary separation
bs_samples_per_condition <- fit_wiener %>%
  # retrieving drift rate values per condition
  emmeans(~prime_emotion * block , dpar = "bs", type = "response") %>%
  #emmeans(~prime_emotion * block, type = "response") %>%
  # retrieving posterior sample for each cell
  # gather_emmeans_draws()
  summary() %>%
  mutate(parameter = "bs")

# retrieving per-condition estimates for the non-decision time
ndt_samples_per_condition <- fit_wiener %>%
  # retrieving drift rate values per condition
  emmeans(~prime_emotion * block , dpar = "ndt", type = "response") %>%
  #emmeans(~prime_emotion * block, type = "response") %>%
  # retrieving posterior sample for each cell
  # gather_emmeans_draws()
  summary() %>%
  mutate(parameter = "ndt")

# retrieving per-condition estimates for the starting point
bias_samples_per_condition <- fit_wiener %>%
  # retrieving drift rate values per condition
  emmeans(~prime_emotion * block , dpar = "bias", type = "response") %>%
  #emmeans(~prime_emotion * block, type = "response") %>%
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
  #facet_wrap(~ group * parameter, scales = "free") +
  facet_wrap(~ parameter, scales = "free")+
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
  #specs = pairwise ~ prime_emotion | block | group,
  specs = pairwise ~ prime_emotion | block,
  dpar = "bias"
)
posterior_emmgrid <- emmeans(
  object = fit_wiener,
  #specs = pairwise ~ prime_emotion | block | group,
  specs = pairwise ~ prime_emotion | block,
  var = 'group',
  dpar = "bias"
)
model_bfs_bias <- bayesfactor_parameters(posterior = posterior_emmgrid, prior = prior_emmgrid)
print(model_bfs_bias)


model_bfs_bias$BF <- exp(model_bfs_bias$log_BF)
model_bfs_bias$BF <- sprintf("%.4f", model_bfs_bias$BF)
model_bfs_bias <- model_bfs_bias %>% rename(BF_bias = BF)
model_bfs_bias <- model_bfs_bias %>% select(-log_BF)

emm_bias <- emmeans(
  object = fit_wiener,
  specs  = pairwise ~ prime_emotion | block + group,
  dpar='bias'
)

bias_by_group <- summary(emm_bias$contrasts, level = 0.95)


### model for drift rate
prior_emmgrid <- emmeans(
  object = prior_model,
  specs = pairwise ~ prime_emotion | block | group,
)
posterior_emmgrid <- emmeans(
  object = fit_wiener,
  specs = pairwise ~ prime_emotion | block | group,
)
model_bfs_drift <- bayesfactor_parameters(posterior = posterior_emmgrid, prior = prior_emmgrid)
print(model_bfs_drift)

model_bfs_drift$BF <- exp(model_bfs_drift$log_BF)
model_bfs_drift$BF <- sprintf("%.4f", model_bfs_drift$BF)
model_bfs_drift <- model_bfs_drift %>% rename(BF_drift = BF)
model_bfs_drift <- model_bfs_drift %>% select(-log_BF)

bf_table <- merge(model_bfs_drift, model_bfs_bias, by='Parameter')



############################
# 2: analyses for mood*group
############################

# # load model 
fit_wiener <- readRDS(file = "models/ddm_patients_mood_16-05-25.rds")

# retrieving predictions for the drift rate
conditional_effects(
  x = fit_wiener,
  negative_rt = TRUE,
  method = "posterior_linpred"
)

# retrieving predictions for the bias (starting position)
conditional_effects(
  x = fit_wiener,
  negative_rt = TRUE,
  dpar = "bias",
  method = "posterior_linpred"
)

#####################
#### calculate BF ###
#####################

# # specify priors for model
prior_model <- unupdate(fit_wiener)

prior_trends <- emtrends(
  object = prior_model,
  specs = ~ group,
  var = "mood",
)
posterior_trends <- emtrends(
  object = fit_wiener,
  specs = ~ group,
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
  specs = ~ group,
  var = "mood",
  dpar='bias'
)
posterior_trends <- emtrends(
  object = fit_wiener,
  specs = ~ group,
  var = "mood",
  dpar='bias'
)
model_mood_bias <- bayesfactor_parameters(posterior = posterior_trends, prior = prior_trends)
print(model_mood_bias)

model_mood_bias$BF <- exp(model_mood_bias$log_BF)
model_mood_bias$BF <- sprintf("%.4f", model_mood_bias$BF)
model_mood_bias <- model_mood_bias %>% rename(BF_bias = BF)
model_mood_bias <- model_mood_bias %>% select(-log_BF)

posterior_summary <- summary(posterior_trends, point.est = "mean")


bf_table_mood_control <- merge(model_mood_drift, model_mood_bias, by='Parameter')


