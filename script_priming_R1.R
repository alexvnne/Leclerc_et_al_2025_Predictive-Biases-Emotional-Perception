# Analyses for "Predictive Biases in Emotional Perception: Differential Influence
# of Mood and Affective Primes in Individuals with and without Mood Disorders
# Written by Alexane Leclerc - 06.02.2026

# clear work space
rm(list = ls())

# set working directory
setwd('path...')

# loading / installing Packages
requiredPackages <- c('plyr','dplyr','tidyr','sjPlot' ,'ggplot2','r2glmm'
                      ,'lme4', 'lmerTest' ,'multcomp', 'broom', 'RLRsim', 'esquisse',
                      'grid', 'gridExtra', 'scales','car','data.table','emmeans', 'tidyverse','ggrain',
                      'scales','car','sm','forestplot','cowplot','psych','ggsignif') 
idx <- which(!requiredPackages %in% installed.packages())
if (length(idx)>0) {
  install.packages(requiredPackages[idx])
}
for (pkg in requiredPackages) {
  require(pkg, character.only = TRUE)
}



# Function to sort and extract data ####
read_trials_folder <- function(folder, add_group_from_filename = FALSE, default_group = "HC") {
  files <- list.files(folder, pattern="\\.csv$", full.names=TRUE)
  if (length(files) == 0) stop(paste("Aucun csv trouvé dans:", folder))
  
  columns_to_extract <- c(
    'participant','frameRate','responseMain.keys','responseMain.rt',
    'mask.tStart','block','prime_exp.tStart','supraliminal2.tStart',
    'main_block.thisN','OS','stim_file','prime_emotion','target_emotion',
    'slider_panas.response','slider_panas.rt','words','valence',
    'slider_mathys1.response','slider_mathys1.rt','item_type','direction',
    'slider_mathys2.response','slider_mathys2.rt','words_mathys',
    'mathys3.text','order','response.keys','key_resp_2.keys',
    'gender_target','target','left_side','right_side',
    'mask_8.tStart','prime_exp_7.tStart'
  )
  
  out <- data.table()
  
  for (f in files) {
    df <- fread(f, na.strings=c('','NA'))
    
    # --- Group
    if (add_group_from_filename) {
      fn <- basename(f)
      if (startsWith(fn, "mdd")) {
        df$Group <- "MDD"
      } else if (startsWith(fn, "bp")) {
        df$Group <- "BD"
      } else {
        df$Group <- "Control"
      }
    } else {
      df$Group <- default_group
    }
    
    # --- Recode keys
    if (!"vb_positive" %in% names(df)) stop(paste("vb_positive manquant dans", f))
    
    if (df$vb_positive[1] == 1) {
      df$responseMain.keys <- ifelse(df$responseMain.keys == 'v', 0, 1)
      df$response.keys     <- ifelse(df$response.keys == 'v', 0, 1)
      df$key_resp_2.keys   <- ifelse(df$key_resp_2.keys == 'v', 1, 0)
    } else {
      df$responseMain.keys <- ifelse(df$responseMain.keys == 'b', 0, 1)
      df$response.keys     <- ifelse(df$response.keys == 'b', 0, 1)
      df$key_resp_2.keys   <- ifelse(df$key_resp_2.keys == 'b', 1, 0)
    }
    
    df$responseMain.keys <- as.numeric(df$responseMain.keys)
    df$response.keys     <- as.numeric(df$response.keys)
    df$key_resp_2.keys   <- as.numeric(df$key_resp_2.keys)
    
    # --- Select columns
    keep <- intersect(columns_to_extract, names(df))
    selected <- df[, ..keep]
    
    # Check whether group is present in selected
    if (!"Group" %in% names(selected)) selected[, Group := df$Group]
    
    out <- rbind(out, selected, fill=TRUE)
  }
  
  out
}

# Apply to HC and patients and combine them
df.patients <- read_trials_folder("data_patients/", add_group_from_filename = TRUE)
df.hc       <- read_trials_folder("data/", add_group_from_filename = FALSE, default_group = "HC")

df.all <- rbind(df.patients, df.hc, fill = TRUE)
colnames(df.all)[colnames(df.all) == 'participant'] <- 'ID'



# Functions to read demographics ####
read_demo_patients_folder <- function(folder) {
  files <- list.files(folder, pattern="\\.csv$", full.names=TRUE)
  out <- data.table()
  
  for (f in files) {
    df <- read.csv(f, na.strings=c('',' ','NA'))
    
    cols <- c('Participant.id','Time.taken','Age','Sex','Ethnicity.simplified','Nationality')
    keep <- intersect(cols, names(df))
    tmp <- as.data.table(df[, keep, drop=FALSE])
    
    out <- rbind(out, tmp, fill=TRUE)
  }
  
  setnames(out, old="Participant.id", new="ID", skip_absent=TRUE)
  out
  
  out <- out[Age != "CONSENT_REVOKED"]
}

read_demo_hc_file <- function(file) {
  df <- read.csv(file, na.strings=c('',' ','NA'))
  df <- df[df$Completion.code=='CMCK1AHR' & df$Status=='APPROVED',]
  
  cols <- c('Participant.id','Time.taken','Age','Sex','Ethnicity.simplified',
            'Nationality')
  keep <- intersect(cols, names(df))
  out <- as.data.table(df[, keep, drop=FALSE])
  
  setnames(out, old="Participant.id", new="ID", skip_absent=TRUE)
  out
}

demo.patients <- read_demo_patients_folder("demographics_patients/")
demo.hc <- read_demo_hc_file("prolific_export_65df4ddfffb481c878ff618a.csv")

demo.all <- rbind(demo.patients, demo.hc, fill=TRUE)



# Prepare final dataframe ####

df.all <- df.all %>%
  rename(Response = responseMain.keys,
         RT = responseMain.rt,
         Target = target_emotion,
         Prime = prime_emotion,
         trial = main_block.thisN,
         Block = block,
         Valence = response.keys,
         Seen = key_resp_2.keys)

df.all <- df.all %>%
  group_by(ID) %>%
  mutate(stim_file = first(stim_file)) %>%
  ungroup()


# rename values in columns and define factors
df.all <- df.all %>%
  mutate(
    Target = fct_recode(Target,
                        Sad   = "sad",
                        Happy = "happy",
                        Morph = "morph"),
    Prime = fct_recode(Prime,
                       Sad     = "sad",
                       Happy   = "happy",
                       Morph   = "morph",
                       Control = "control"),
    Block = fct_recode(Block,
                       Masked = 'subliminal',
                       Unmasked = 'supraliminal'),
    
    ID        = factor(ID),
    stim_file = factor(stim_file),
    Block     = factor(Block)
  )

# convert RT in ms
df.all$RT <- df.all$RT*1000

# get target data
df.target <- df.all[!is.na(df.all$Response),]

df.target <- df.target %>%
  mutate(
    Correct = case_when(
      Target == "Morph" ~ NA,                       # always NA
      Target == "Sad"   & Response == 0 ~ 1,        # correct
      Target == "Sad"   & Response == 1 ~ 0,        # incorrect
      Target == "Happy" & Response == 1 ~ 1,        # correct
      Target == "Happy" & Response == 0 ~ 0,        # incorrect
      TRUE ~ NA                                      # fallback (shouldn't be needed)
    )
  )

# number of participants per group
df.target %>%
  distinct(ID, Group) %>%
  count(Group)


# Exclusions ####

# get unique id per line to get rid of trials later
df.target <- df.target %>%
  mutate(row_extract = row_number())

# 1) Issues with presentation times

df.target <- df.target[!df.target$RT > 3000,]

df.target$primeSubli <- df.target$mask.tStart - df.target$prime_exp.tStart
df.target$primeSupra <- df.target$supraliminal2.tStart - df.target$prime_exp.tStart

# check for trials that were below 20ms and above 45ms
out_subli <- subset(df.target, Block == 'Masked' & (primeSubli <= 0.020 | primeSubli >= 0.045))
out_supra <- subset(df.target, Block == 'Unmasked' & (primeSupra <= 0.020 | primeSupra >= 0.045))

count_subli <- as.data.frame(table(out_subli$ID))
count_supra <- as.data.frame(table(out_supra$ID))

out_timing <- unique(c(
  as.character(count_subli$Var1[count_subli$Freq > 0.15 * 432]),
  as.character(count_supra$Var1[count_supra$Freq > 0.15 * 432])
))

df.target <- df.target %>% filter(!ID %in% out_timing)

df.target <- df.target %>%
  filter(((Block == "Masked" & primeSubli > 0.025 & primeSubli < 0.040) | 
            (Block == "Unmasked" & primeSupra > 0.025 & primeSupra < 0.040)))

# 2) Aberrant ratings of happy and sad targets

out_rating <- df.target %>%
  filter(Target %in% c("Happy","Sad"), Group %in% c("MDD","BD","HC")) %>%
  group_by(Group, Target, ID) %>%
  summarise(Response = mean(Response, na.rm = TRUE), .groups = "drop") %>%
  group_by(Group, Target) %>%
  mutate(Z = as.numeric(scale(Response))) %>%
  ungroup() %>%
  filter(
    (Target == "Happy" & Z < -2) |
      (Target == "Sad"   & Z >  2)
  ) %>%
  count(Group, ID) #%>%       # per group


# 3) Aberrant RTs to targets

out_RT <- df.target %>%
  filter(Target %in% c("Happy","Sad"),
         Group %in% c("MDD","BD","HC")) %>%
  group_by(Group, ID) %>%
  summarise(meanRT = mean(RT, na.rm = TRUE), .groups = "drop") %>%
  group_by(Group) %>%
  mutate(Z_RT = as.numeric(scale(meanRT))) %>%
  ungroup() %>%
  filter(abs(Z_RT) > 2) %>%
  #dplyr::select(Group, ID) %>%
  distinct()

all_outliers <- merge(out_RT,out_rating, all=T)
df.target <- df.target %>% filter(!ID %in% all_outliers$ID)

# 4) Exclude trials above or below 3 SD from the mean for each target and group

trials_out <- df.target %>%
  filter(Group %in% c("HC","BD","MDD"),
         Target %in% c("Happy","Sad","Morph")) %>%
  group_by(Group, Target) %>%
  mutate(
    Z_RT = as.numeric(scale(RT)),
    outlier_trial = abs(Z_RT) > 3
  ) %>%
  ungroup() %>%
  filter(outlier_trial)

df.target <- df.target[!df.target$row_extract %in% trials_out$row_extract,]


# Load PCA ####
pca <- read.table('path to pca table')
colnames(pca)[colnames(pca) == 'participant_id'] <- 'ID'
columns_to_extract <- c('ID','Group','PC1','PC2','PC3','PC4')
pca <- pca[, columns_to_extract]

df.target <- merge(df.target, pca, by=c('ID','Group'), all=T)
df <- merge(df.target, demo.all, by='ID', all=T)


ggplot(df, aes(x = PC1)) +
  geom_density(alpha = 0.4)

ggplot(df, aes(sample = PC1)) +
  stat_qq() +
  stat_qq_line() +
  facet_wrap(~ Group)


# Dataframe for Bayesian Modeling ####
df <- df.target %>%
  left_join(demo.all, by = "ID")
df <- df[df$Prime!='Morph',]
df.ddm <- df[df$Target=='Morph',]

################################################################################
################################################################################

# sort questionnaire data

# PANAS
df.panas <- df.all[!is.na(df.all$slider_panas.response),]
extract_panas <- c('slider_panas.response','ID','words','Group','valence')
df.panas <- df.panas[,extract_panas]
df.panas <- df.panas[(df.panas$ID %in% df.target$ID),]

summary_panas <- df.panas %>%
  group_by(ID, valence, Group) %>%
  summarise(Total_Score = sum(slider_panas.response))
summary_panas <- pivot_wider(summary_panas, names_from = valence, values_from = Total_Score)

panas_pca <- pivot_wider(df.panas, names_from = words, values_from = slider_panas.response)


# MATHYS
df.mathys1 <- df.all[!is.na(df.all$slider_mathys1.response),]
extract_mathys1 <- c('ID','slider_mathys1.response','item_type','direction','left_side','right_side','Group')
df.mathys1 <- df.mathys1[,extract_mathys1]
# reverse scores
df.mathys1$slider_mathys1.response[df.mathys1$direction == "inverse"] <- 10 - df.mathys1$slider_mathys1.response[df.mathys1$direction == "inverse"]
df.mathys1 <- df.mathys1 %>% mutate(item = ifelse(direction == "normal", right_side, left_side))
extract_mathys1 <- c('ID','slider_mathys1.response','item','Group','item_type')
df.mathys1 <- df.mathys1[,extract_mathys1]
df.mathys1 <- df.mathys1[(df.mathys1$ID %in% df.target$ID),]

summary_mathys1 <- df.mathys1 %>%
  group_by(ID,item_type, Group) %>%
  summarise(Total_Score = sum(slider_mathys1.response))

# pivot the data to have a separate columns for all item types
summary_mathys1 <- pivot_wider(summary_mathys1, names_from = item_type, values_from = Total_Score)
summary_mathys1$Total_Mathys <- summary_mathys1$CO + summary_mathys1$VO + summary_mathys1$MO + 
  summary_mathys1$EM + summary_mathys1$SE 

mathys1_pca <- pivot_wider(df.mathys1, names_from = item, values_from = slider_mathys1.response)


df.mathys2 <- df.all[!is.na(df.all$slider_mathys2.response),]
extract_mathys2 <- c('ID','slider_mathys2.response','words_mathys','Group')
df.mathys2 <- df.mathys2[,extract_mathys2]
df.mathys2 <- df.mathys2[(df.mathys2$ID %in% df.target$ID),]
summary_mathys2 <- pivot_wider(df.mathys2, names_from = words_mathys, values_from = slider_mathys2.response)

mathys2_pca <- summary_mathys2

summary_questionnaires <- merge(summary_mathys1, summary_mathys2, by = c('ID','Group'))
summary_questionnaires <- merge(summary_questionnaires, summary_panas, by = c('ID','Group'))


summary_questionnaires$Group <- factor(summary_questionnaires$Group, levels = c("HC", "BD", "MDD"))

quest <- summary_questionnaires %>% 
  group_by(Group) %>% 
  summarise(across(where(is.numeric),
                   ~sprintf("%.2f ± %.2f", mean(.x, na.rm = TRUE), sd(.x, na.rm = TRUE)),
                   .names = "{.col}")) %>% 
  pivot_longer(-Group, names_to = "Variable", values_to = "value") %>% 
  pivot_wider(names_from = Group, values_from = value)

write_xlsx(results_table, file.path(fig_outdir, "table_quest.xlsx"))

# t-tests and multiple comparisons ? 

model <- aov(Negative ~ Group, data = summary_questionnaires)
summary(model)
emmeans(model, pairwise ~ Group, adjust = 'tukey')


# variables à tester
vars <- summary_questionnaires %>%
  dplyr::select(-Group) %>%
  dplyr::select(where(is.numeric)) %>%
  names()

results_table <- map_dfr(vars, function(var) {
  
  form <- as.formula(paste(var, "~ Group"))
  model <- aov(form, data = summary_questionnaires)
  
  # test global
  anova_tidy <- tidy(model) %>% filter(term == "Group")
  
  anova_res <- tibble(
    Variable = var,
    df1 = anova_tidy$df[1],
    df2 = df.residual(model),
    F = anova_tidy$statistic[1],
    p_global = anova_tidy$p.value[1]
  )
  
  # descriptifs
  descriptives <- summary_questionnaires %>%
    group_by(Group) %>%
    summarise(
      mean = mean(.data[[var]], na.rm = TRUE),
      sd = sd(.data[[var]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(value = sprintf("%.2f ± %.2f", mean, sd)) %>%
    dplyr::select(Group, value) %>%
    pivot_wider(names_from = Group, values_from = value)
  
  # post-hoc Tukey
  posthoc <- emmeans(model, pairwise ~ Group, adjust = "tukey")$contrasts %>%
    summary() %>%
    as_tibble() %>%
    dplyr::select(contrast, p.value) %>%
    mutate(contrast = gsub(" - ", "_vs_", contrast)) %>%
    pivot_wider(names_from = contrast, values_from = p.value, names_prefix = "p_")
  
  bind_cols(anova_res, descriptives, posthoc)
})

results_table


# Add the proportion of correct trials across conditions and groups
df.hs <- df[df$Target!='Morph',]
df.hs <- df.hs[df.hs$Prime %in% c('Sad', 'Happy'), ]
df.hs$Congruence <- ifelse((df.hs$Target == 'Happy' & df.hs$Prime == 'Happy') |
                             (df.hs$Target == 'Sad' & df.hs$Prime == 'Sad'),
                           'Congruent', 'Incongruent')

df_prop <- df.hs %>%
  group_by(Block, Group, Congruence, Target) %>%
  summarise(
    n_trials = n(),
    prop_correct = mean(Correct, na.rm = TRUE),
    .groups = "drop"
  )

library(writexl)
fig_outdir <- "fig outdir"
write_xlsx(df_prop, file.path(fig_outdir, "table_prop_correct.xlsx"))



# Detection task ####

df.detection <- df.all[!is.na(df.all$Valence),]
extract_detection <- c('ID','Valence','Seen','Prime','Target',
                       'prime_exp_7.tStart','mask_8.tStart','Group')
df.detection <- df.detection[,extract_detection]
df.detection <- df.detection[df.detection$ID %in% df$ID,]

df.detection <- df.detection %>%
  mutate(row_extract = row_number())

# check presentation time
df.detection$pres_prime <- df.detection$mask_8.tStart - df.detection$prime_exp_7.tStart
out_det <- subset(df.detection, pres_prime <= 0.020 | pres_prime >= 0.045)
df.detection <- df.detection[!df.detection$row_extract %in% out_det$row_extract,]

### dprime

df.detection$Stim <- ifelse(df.detection$Prime %in% c("Happy", "Sad", "Morph"), "Target", "Control")
df.detection$Participant <- df.detection$ID

std <- function(x) sd(x)/sqrt(length(x))

dprime <- function(hit,fa) {
  qnorm(hit) - qnorm(fa)
}

beta <- function(hit,fa) {
  zhr <- qnorm(hit)
  zfar <- qnorm(fa)
  exp(-zhr*zhr/2+zfar*zfar/2)
}

#create a function to compute dprime at the group level

compdprime <- function(subdata,anscol, ans, stimcol, stimsame, stimdiff) {
  subdata <- as.data.frame(subdata)
  subdata$hit <- ifelse(subdata[stimcol] == stimsame & subdata[anscol] ==  ans, 1,0)
  subdata$fa <- ifelse(subdata[stimcol] == stimdiff & subdata[anscol] == ans, 1,0)
  if (length(unique(subdata$Participant))==1) {
    hit.rate <- mean(subdata[subdata[stimcol]==stimsame,]$hit)
    fa.rate <- mean(subdata[subdata[stimcol]==stimdiff,]$fa)
  }
  else{
    hit.rate <- with(subdata[subdata[stimcol]==stimsame,], tapply(hit, Participant, mean))
    fa.rate <- with(subdata[subdata[stimcol]==stimdiff,], tapply(fa, Participant, mean))
  }
  
  # apply corrections
  hit.rate[(hit.rate)==1]<- (nrow(subdata[subdata[stimcol]==stimsame,])-0.5)/nrow(subdata[subdata[stimcol]==stimsame,])
  hit.rate[(hit.rate)==0]<- (0.5)/nrow(subdata[subdata[stimcol]==stimdiff,])
  fa.rate[(fa.rate)==0]<- (0.5)/nrow(subdata[subdata[stimcol]==stimdiff,])
  fa.rate[(fa.rate)==1]<- (nrow(subdata[subdata[stimcol]==stimsame,])-0.5)/nrow(subdata[subdata[stimcol]==stimsame,])
  
  
  return(list(dprime(hit.rate, fa.rate),hit.rate,fa.rate, beta(hit.rate, fa.rate)))
  
}


# before calculating dprime, check who has aberrant responses to detection task
# check proportion of seen trials per ID then per prime
seenID <- aggregate(df.detection, Seen ~ ID + Group + Stim, mean)

ggplot(seenID, aes(x = (ID), y = Seen, fill = Stim)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "ID", y = "Proportion Seen", title = "Proportion of Seen per ID by Stim Type") +
  scale_fill_manual(values = c("Target" = "#69b3a2", "Control" = "#f6ae2d")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# exclude from visibility analyses those who answered always seen or have control > target
seenID.wide <- seenID %>%
  pivot_wider(names_from = Stim, values_from = Seen)

seenID.wide$diff <- seenID.wide$Control - seenID.wide$Target
seen.out <- seenID.wide[seenID.wide$diff >= 0 | seenID.wide$Control>0.75,]
seen.out <- seen.out[seen.out$Target!=0 | seen.out$Control!=0,]

df.detection <- df.detection[!df.detection$ID %in% seen.out$ID,]


# now, calculate visibility dprime 

subjdprime_visibility <- data.frame(compdprime(df.detection,'Seen', 1, 'Stim' , 'Target', 'Control')[1])
names(subjdprime_visibility) <- "dprime"
subjdprime_visibility['hit'] <- data.frame(compdprime(df.detection,'Seen', 1, 'Stim' , 'Target', 'Control')[2])
subjdprime_visibility['fa'] <- data.frame(compdprime(df.detection,'Seen', 1, 'Stim' , 'Target', 'Control')[3])
subjdprime_visibility$Participant <- as.factor(rownames(subjdprime_visibility))
subjdprime_visibility <- subjdprime_visibility[subjdprime_visibility$Participant %in% df.detection$ID,]
subjdprime_visibility$ID <- subjdprime_visibility$Participant

ggplot(subjdprime_visibility, aes(x = dprime)) +
  geom_histogram(aes(y = stat(density)), bins = 30, fill = "#69b3a2", alpha = 0.7) +
  geom_density(color = "#f6ae2d", alpha = 0.8) +
  ggtitle("Distribution of d-prime visibility per participant") +
  theme_minimal()

# mean dp
mean(subjdprime_visibility$dprime, na.rm=T) # 3.19

# compute accuracy
accuracy <- subjdprime_visibility %>%
  summarise(
    mean_accuracy = mean(hit * 100, na.rm = TRUE) # 75.6.
  )

# t-test
t_test_dprime <- t.test(subjdprime_visibility$dprime, mu = 0)
print(t_test_dprime) # significant


### 2 : happy versus sad prime
################################

df.dprime.HS <- df.detection[df.detection$Prime== 'Happy' | df.detection$Prime== 'Sad' ,]
df.dprime.HS <- df.dprime.HS[df.dprime.HS$Target== 'Happy' | df.dprime.HS$Target== 'Sad' ,]

valenceID <- aggregate(df.dprime.HS, Valence ~ ID + Prime + Target , mean)
ggplot(valenceID)+aes(x=Target,y=Valence, fill=Prime)+geom_violin()
model <- glmer(Valence ~ Prime * Target + (1 | ID), data= df.dprime.HS, family = binomial(link='logit'))
summary(model)

# now calculate d' for happy vs. sad prime perception
subjdprime_happySad <- data.frame(compdprime(df.detection,'Valence', 1, 'Prime' , 'Happy', 'Sad')[1])
names(subjdprime_happySad) <- "dprime"
subjdprime_happySad['hit'] <- data.frame(compdprime(df.detection,'Valence', 1, 'Prime' , 'Happy', 'Sad')[2])
subjdprime_happySad['fa'] <- data.frame(compdprime(df.detection,'Valence', 1, 'Prime' , 'Happy', 'Sad')[3])
subjdprime_happySad$Participant <- as.factor(rownames(subjdprime_happySad))
subjdprime_happySad <- subjdprime_happySad[subjdprime_happySad$Participant %in% df.detection$ID,]
subjdprime_happySad$ID <- subjdprime_happySad$Participant
ggplot(subjdprime_happySad, aes(x = dprime)) +
  geom_histogram(aes(y = stat(density)), bins = 30, fill = "#69b3a2", alpha = 0.7) +
  geom_density(color = "#f6ae2d", alpha = 0.8) +
  ggtitle("Distribution of d-prime happy/sad per participant") +
  theme_minimal()

# mean happy-sad dprime
mean(subjdprime_happySad$dprime, na.rm=T) # 1.63

# compute accuracy
accuracy <- subjdprime_happySad %>% 
  summarise(mean_accuracy = mean(hit * 100, na.rm = TRUE)) # 72.68

print(t.test(subjdprime_happySad$dprime, mu = 0)) # significant

# d' happy vs. sad par groupe
groups <- unique(df.detection$Group)

dprime_by_group <- lapply(groups, function(g) {
  subdata <- df.detection[df.detection$Group == g, ]
  
  dprime_vals <- data.frame(
    compdprime(subdata, 'Valence', 1, 'Prime', 'Happy', 'Sad')[1]
  )
  names(dprime_vals) <- "dprime"
  dprime_vals$ID <- rownames(dprime_vals)
  dprime_vals$Group <- g
  return(dprime_vals)
})

dprime_by_group <- do.call(rbind, dprime_by_group)

# Moyenne et t-test par groupe
dprime_by_group %>%
  group_by(Group) %>%
  summarise(
    mean_dprime = mean(dprime, na.rm = TRUE),
    sd_dprime = sd(dprime, na.rm = TRUE),
    t_stat = t.test(dprime, mu = 0)$statistic,
    p_value = t.test(dprime, mu = 0)$p.value
  )

# d' happy (happy = signal, control = noise)
dprime_happy <- data.frame(compdprime(df.detection, 'Valence', 1, 'Prime', 'Happy', 'Control')[1])
names(dprime_happy) <- "dprime_happy"
dprime_happy$ID <- rownames(dprime_happy)

# d' sad (sad = signal, control = noise)
dprime_sad <- data.frame(compdprime(df.detection, 'Valence', 0, 'Prime', 'Sad', 'Control')[1])
names(dprime_sad) <- "dprime_sad"
dprime_sad$ID <- rownames(dprime_sad)

# Merge et ajouter le groupe
dprime_both <- merge(dprime_happy, dprime_sad, by = "ID")
dprime_both <- merge(dprime_both, unique(df.detection[, c("ID", "Group")]), by = "ID")

# Comparer par groupe
dprime_both %>%
  group_by(Group) %>%
  summarise(
    mean_happy = mean(dprime_happy, na.rm = TRUE),
    mean_sad = mean(dprime_sad, na.rm = TRUE),
    t_stat = t.test(dprime_happy, dprime_sad, paired = TRUE)$statistic,
    p_value = t.test(dprime_happy, dprime_sad, paired = TRUE)$p.value
  )

table_dprime <- dprime_both %>%
  group_by(Group) %>%
  summarise(
    Mean_dprime_happy = round(mean(dprime_happy, na.rm = TRUE), 2),
    SD_dprime_happy = round(sd(dprime_happy, na.rm = TRUE), 2),
    Mean_dprime_sad = round(mean(dprime_sad, na.rm = TRUE), 2),
    SD_dprime_sad = round(sd(dprime_sad, na.rm = TRUE), 2),
    t_stat = round(t.test(dprime_happy, dprime_sad, paired = TRUE)$statistic, 2),
    df = t.test(dprime_happy, dprime_sad, paired = TRUE)$parameter,
    p_value = round(t.test(dprime_happy, dprime_sad, paired = TRUE)$p.value, 3)
  )
write.csv(table_dprime, "path to write table", row.names = FALSE)
