# Script to Analyze Data Healthy Participants - AL
# Healthy participants
# v0.2 : AL 28.05.2025

# clear work space
rm(list = ls())

# set working directory
setwd('path/to/directory')

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


##########################
#### Load data frames ####
##########################

# load priming task data
df.full <- read.csv('/data_priming_task_HC.csv')

# load detection task data
df.detection <- read.csv('/data_detection_task_HC.csv')

# load PCA data
pca <- read.table('/PCA_data.csv')
# only keep data from healthy participants
pca <- pca[pca$Group=='HC',]
# rename the column to ID
colnames(pca)[colnames(pca) == 'participant_id'] <- 'ID'

# add PCA scores to the dataframe
df.full <- merge(df.full, pca, by='ID')


#########################################
#### Define data frames per interest ####
#########################################

# re-define each df per target 
df.happy <- df.full[df.full$Target == 'Happy',]
df.sad <- df.full[df.full$Target == 'Sad',]
df.morph <- df.full[df.full$Target == 'Morph',]

# now df for RT analyses with only correct trials included
df.happy.corr <- df.happy[df.happy$correct.happy==1,]
df.sad.corr <- df.sad[df.sad$correct.sad==1,]



# set palette for graphical representations
lPalette <- c("#2f4858",'#f6ae2d',"#f26419")
rtPalette <- c('#00416A','#86bbd8')



########### Analysis 1: Primes Effect on Morphed targets categorization 
# Graphs:

# revalue subli supra into masked unmasked
df.morph$Block <- plyr::revalue(df.morph$Block, c("subliminal"="Masked","supraliminal"="Unmasked"))

# aggregate data for graphs
agg.morph <- aggregate(df.morph, Response ~ Prime + Block + ID, mean)
# only interested in control, happy, sad
agg.morph <- agg.morph[agg.morph$Prime != 'Morph',]

# plot effect of primes on morph categorization by priming_type 
ggplot(agg.morph) + 
  aes(y = Response, x = Prime, fill = Prime, color = Prime) +
  geom_boxplot(alpha = 0.6, width = 0.4, outlier.shape = NA, 
               position = position_dodge(width = 0.25)) +
  geom_jitter(aes(color = Prime), alpha = 0.6, size = 1.4, width = 0.15) +
  facet_wrap(~ Block, strip.position = 'bottom') +  
  stat_summary(fun = mean, shape = 23, size = 0.8, stroke = 0.8, color = "black", alpha = 0.8) +  
  theme_blank() +  
  scale_color_manual(values = lPalette) +
  scale_fill_manual(values = lPalette) +
  ggtitle("Healthy participants") +
  theme(plot.title = element_text(size = 21),
        axis.title.x = element_text(size = 20),   # x‐axis label
        axis.title.y = element_text(size = 20),   # y‐axis label
        axis.text.x  = element_text(size = 18),   # x‐axis tick labels
        axis.text.y  = element_text(size = 18),
        strip.text.x = element_text(size = 18)) +
  
  theme(legend.position = "none")+
  geom_hline(yintercept = 0.5,
             linetype   = "dashed",
             color      = "grey40",
             size       = 0.45) +
  
  labs(tag = "A") +
  theme(
    plot.tag.position = c(0.02, 0.98),  # x, y in [0,1], (0,1) = top-left
    plot.tag = element_text(size = 21, face = "bold"))+
  
  # Masked
  geom_signif(data = subset(agg.morph, Block == "Masked"), 
              comparisons = list(c("Happy", "Control")), 
              annotations = c("**"), 
              y_position = 1.05, 
              tip_length = 0.02, 
              textsize = 6, 
              vjust = 0.4, 
              color = "black") +
  
  # Unmasked
  geom_signif(data = subset(agg.morph, Block == "Unmasked"), 
              comparisons = list(c("Happy", "Control")), 
              annotations = c("***"), 
              y_position = 1.05, 
              tip_length = 0.02, 
              textsize = 6, 
              vjust = 0.4, 
              color = "black") +
  
  # Unmasked
  geom_signif(data = subset(agg.morph, Block == "Unmasked"), 
              comparisons = list(c("Sad", "Control")), 
              annotations = c("**"), 
              y_position = 1.15, 
              tip_length = 0.02, 
              textsize = 6, 
              vjust = 0.4, 
              color = "black")


### Stats: 

# same here, delete morph primes from the target morph df
df.morph.model <- df.morph[!df.morph$Prime == 'Morph',]

# 1. Model on Masked trials
# only keep masked trials
morph.subli <- df.morph.model[df.morph.model$Block == 'Masked',]
# fit glmm
model.prime.subli <- glmer(Response ~ Prime + (1 | ID), data= morph.subli, family = binomial(link='logit'))
summary(model.prime.subli)
# OR
exp(fixef(model.prime.subli))
# 95 CI for OR values
exp(confint(model.prime.subli, parm = c("PrimeHappy", "PrimeSad"), level = 0.95))

# 2. Model on Unmasked trials
# only keep unmasked trials
morph.supra <- df.morph.model[df.morph.model$Block == 'Unmasked',]
#fit glmm
model.prime.supra <- glmer(Response ~ Prime + (1 | ID), data= morph.supra, family = binomial(link='logit'))
summary(model.prime.supra)
# OR
exp(fixef(model.prime.supra))
# 95 CI for OR values
exp(confint(model.prime.supra, parm = c("PrimeHappy", "PrimeSad"), level = 0.95))


# fit models with covariates: age and sex
model.prime.subli.conf <- glmer(Response ~ Prime + Age + Sex + (1 | ID), data= morph.subli, family = binomial(link='logit'))
summary(model.prime.subli.conf)
model.prime.supra.conf <- glmer(Response ~ Prime + Age + Sex + (1 | ID), data= morph.supra, family = binomial(link='logit'))
summary(model.prime.supra.conf)




########### Analysis 2: Mood Effect on Morphed targets categorization 

# only keep trials with control primes
df.morph.control <- df.morph[df.morph$Prime=='Control',]
# aggregate data
agg.morph.control <- aggregate(df.morph.control, Response ~ ID + PC1, mean)

### Graph
ggplot(df.morph.control, aes(x = PC1, y = Response)) + 
  geom_point() + stat_smooth(method = "glm",method.args = list(family = binomial(link='logit')), se = TRUE, color='#2f4858') + 
  xlab("Mood (PC1)") + ylab("Probability of Positive Answer") + 
  ggtitle("Healthy participants") +
  theme_blank()+
  theme(plot.title = element_text(size = 21),
        axis.title.x = element_text(size = 20),   # x‐axis label
        axis.title.y = element_text(size = 20),   # y‐axis label
        axis.text.x  = element_text(size = 18),   # x‐axis tick labels
        axis.text.y  = element_text(size = 18))+
  
  labs(tag = "A") +
  theme(
    plot.tag.position = c(0.02, 0.98),  # x, y in [0,1], (0,1) = top-left
    plot.tag = element_text(size = 21, face = "bold"))+
  
  annotate("text", x = 4, y = 0.9, label = "p = .025*", 
           size = 6, fontface = "italic")

### Stats:
# fit glmm
model.morph.mood <- glmer(Response ~ PC1 +  (1 | ID), data=df.morph.control,family = binomial(link='logit'))
summary(model.morph.mood)

# same model with age and sex as covariates
model.morph.mood <- glmer(Response ~ PC1 + Age + Sex + (1 | ID), data=df.morph.control,family = binomial(link='logit'))
summary(model.morph.mood)



########### Analysis 3: RTs to targets preceded by congruent primes 

# only keep df of correct responses to happy and sad targets
df.hs <- rbind(df.happy.corr, df.sad.corr)
# keep happy and sad primes only
df.hs <- df.hs[df.hs$Prime %in% c('Sad', 'Happy'), ]
# create a congruence column -> e.g., if prime = happy and target = happy -> congruent
df.hs$Congruence <- ifelse((df.hs$Target == 'Happy' & df.hs$Prime == 'Happy') |
                             (df.hs$Target == 'Sad' & df.hs$Prime == 'Sad'),
                           'Congruent', 'Incongruent')
# only keep columns that we are interested in
extract <- c('Target','Prime','Block','ID','Congruence','RT','Age','Sex')
df.hs <- df.hs[,extract]
# rename subli to masked and supra to unmasked
df.hs$Block <- plyr::revalue(df.hs$Block, c("subliminal"="Masked","supraliminal"="Unmasked"))
# aggregate for graphical representation
agg.df.hs <- aggregate(RT ~ ID + Block + Target + Prime + Congruence, data=df.hs, mean)

### Graph:
ggplot(agg.df.hs) + 
  aes(y = RT, x = Congruence, fill = Congruence, color = Congruence) +
  geom_violin(alpha = 0.5) +
  labs(x = "Prime Congruence", y = "RT to Targets (ms)", 
       title = 'Healthy participants') +
  facet_wrap(~ Block, strip.position = 'bottom') +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, color = "black", fatten = 2) +  
  geom_jitter(width = 0.2, alpha = 1) +
  theme_blank() +  
  theme(
    #plot.title = element_text(hjust = 0.5, size = 12), 
    legend.position = "none",
    plot.title = element_text(size = 21, margin = margin(b = 35)),
    axis.title.x = element_text(size = 20),   # x‐axis label
    axis.title.y = element_text(size = 20),   # y‐axis label
    axis.text.x  = element_text(size = 18),   # x‐axis tick labels
    axis.text.y  = element_text(size = 18),
    strip.text.x = element_text(size = 18)
  ) +
  scale_color_manual(values = rtPalette) +
  scale_fill_manual(values = rtPalette) +
  
  labs(tag = "A") +
  theme(
    plot.tag.position = c(0.02, 0.98),  # x, y in [0,1], (0,1) = top-left
    plot.tag = element_text(size = 21, face = "bold"))+
  
  # Add significance bars for each condition separately
  # Masked
  geom_signif(data = subset(agg.df.hs, Block == "Masked"), 
              comparisons = list(c("Congruent", "Incongruent")), 
              annotations = c("***"),  # Adjust p-value annotation
              y_position = 820,        # Adjust height of significance bar
              tip_length = 0.02, 
              textsize = 6, 
              vjust = 0.4, 
              size=0.8,
              color = "black")+
  
  # Unmasked
  geom_signif(data = subset(agg.df.hs, Block == "Unmasked"), 
              comparisons = list(c("Congruent", "Incongruent")), 
              annotations = c("***"),  
              y_position = 820,        
              tip_length = 0.02, 
              textsize = 6, 
              vjust = 0.4, 
              size=0.8,
              color = "black") 

### Stats:
# fit a lmm
model.rt <- lmer(RT ~ Congruence * Block * Target + (1 | ID), data=df.hs)
summary(model.rt)
# planned comparisons
print(emmeans(model.rt, pairwise ~ Congruence | Block)) # priming effect in masked then unmasked trials
contrast(emmeans(model.rt, ~ Congruence * Block), interaction = "pairwise") # stronger priming effect in masked or unmasked? 

# get summary stats: means and sd per priming type and congruence
summary_stats <- df.hs %>%
  group_by(Block,Congruence) %>%
  summarise(
    mean = mean(RT),
    se = sd(RT) / sqrt(n())
  )


# same model with covariates: age + sex
model.rt.conf <- lmer(RT ~ Congruence * Block * Target + Age + Sex + (1 | ID), data=df.hs)
summary(model.rt.conf)
print(emmeans(model.rt.conf, pairwise ~ Congruence | Block)) # priming effect in masked then unmasked trial
contrast(emmeans(model.rt.conf, ~ Congruence * Block), interaction = "pairwise") # stronger priming effect in masked or unmasked? 






########### Analysis 4: RTs to morphed targets preceded by primes congruent with response

# only keep happy and sad primes
df.morph.hs <- df.morph[df.morph$Prime %in% c('Happy', 'Sad'), ]
# define congruence -> e.g., if response to morph = 1 and prime = happy then it is congruent
df.morph.hs$Congruence <- ifelse((df.morph.hs$Response == 1 & df.morph.hs$Prime == 'Happy') |
                                   (df.morph.hs$Response == 0 & df.morph.hs$Prime == 'Sad'),
                                 'Congruent', 'Incongruent')
# aggregate for graphical representation
agg.morph <- aggregate(RT ~ ID + Prime + Block + Response + Congruence, data=df.morph.hs, mean)

### Graphs: 
ggplot(agg.morph) + aes(y=RT, x= Congruence, fill=Congruence, color=Congruence)+
  geom_violin(alpha=0.5)+
  labs(x = "Prime", y = "RT (ms)", title = 'Healthy participants')+
  facet_wrap( ~ Block, strip.position = 'bottom')+
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, color = "black", fatten = 2) +  
  geom_jitter(width = 0.2, alpha =1)+
  theme_blank()+
  theme(
    legend.position = "none",
    plot.title = element_text(size = 21, margin = margin(b = 35)),
    axis.title.x = element_text(size = 20),   # x‐axis label
    axis.title.y = element_text(size = 20),   # y‐axis label
    axis.text.x  = element_text(size = 20),   # x‐axis tick labels
    axis.text.y  = element_text(size = 18),
    strip.text.x = element_text(size = 18)
  ) +
  scale_color_manual(values=rtPalette)+
  scale_fill_manual(values = rtPalette)+
  
  # Add significance bars for each condition separately
  # Masked
  geom_signif(data = subset(agg.morph, Block == "Masked"), 
              comparisons = list(c("Congruent", "Incongruent")), 
              annotations = c("***"),  
              y_position = 1300,        
              tip_length = 0.02, 
              textsize = 6, 
              vjust = 0.4, 
              size = 0.8,
              color = "black")


# summary statistics
summary_stats_morph <- agg.morph %>%
  group_by(Block,Congruence) %>%
  summarise(
    mean = mean(RT),
    se = sd(RT) / sqrt(n())
  )


# check for outliers (as seen on the graph)
agg.morph$Z <- scale(agg.morph$RT)
agg.morph.out <- agg.morph[agg.morph$Z < (-3) |agg.morph$Z > 3 ,] # one outlier, 3 SD above the mean

# create a new df for morphs without outliers
agg.morph.NoOutliers <- agg.morph[!agg.morph$ID %in% agg.morph.out$ID,]

# new graph without outlier
ggplot(agg.morph.NoOutliers) + aes(y=RT, x= Congruence, fill=Congruence, color=Congruence)+
  geom_violin(alpha=0.5)+
  labs(x = "Prime Congruence", y = "RT (ms)", title = 'Healthy participants')+
  facet_wrap( ~Block, strip.position = 'bottom')+
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, color = "black", fatten = 2) +  
  geom_jitter(width = 0.2, alpha =1)+
  theme_blank()+
  theme(
    legend.position = "none",
    plot.title = element_text(size = 21, margin=margin(b=35)),
    axis.title.x = element_text(size = 20),   # x‐axis label
    axis.title.y = element_text(size = 20),   # y‐axis label
    axis.text.x  = element_text(size = 18),   # x‐axis tick labels
    axis.text.y  = element_text(size = 18),
    strip.text.x = element_text(size = 18)
  ) +
  #plot.title = element_text(hjust = 0.5, size = 12), legend.position = "none") +
  scale_color_manual(values=rtPalette)+
  scale_fill_manual(values = rtPalette)+
  coord_cartesian(ylim = c(430, 1350))+
  
  labs(tag = "A") +
  theme(
    plot.tag.position = c(0.02, 0.98),  # x, y in [0,1], (0,1) = top-left
    plot.tag = element_text(size = 21, face = "bold"))+
  
  # Add significance bars for each condition separately
  # Masked
  geom_signif(data = subset(agg.morph, Block == "Masked"), 
              comparisons = list(c("Congruent", "Incongruent")), 
              annotations = c("***"),  
              y_position = 1050,        
              tip_length = 0.02, 
              textsize = 6, 
              vjust = 0.4, 
              size=0.8,
              color = "black")

# new summary statistics without outliers
summary_stats_morph <- agg.morph.NoOutliers %>%
  group_by(Block,Congruence) %>%
  summarise(
    mean = mean(RT),
    se = sd(RT) / sqrt(n())
  )


### Stats:
# with outliers
# fit lmm
model.morph <- lmer(RT ~ Congruence * Block * Response + (1 | ID), data=df.morph.hs)
summary(model.morph)
print(emmeans(model.morph, pairwise ~ Congruence | Block)) # priming effect in masked then unmasked trial
contrast(emmeans(model.morph, ~ Congruence * Block), interaction = "pairwise") # stronger priming effect in masked or unmasked trials? 

# without outliers
df.morph.hs.noOut <- df.morph.hs[!df.morph.hs$ID %in% agg.morph.out$ID,]
model.morph.out <- lmer(RT ~ Congruence * Block * Response + (1 | ID), data=df.morph.hs.noOut)
summary(model.morph.out)
print(emmeans(model.morph.out, pairwise ~ Congruence | Block)) # priming effect in masked then unmasked trial
contrast(emmeans(model.morph.out, ~ Block * Congruence), interaction = "pairwise") # stronger priming effect in masked or unmasked trials? 

# same model with covariates: Age + Sex
model.morph.out.conf <- lmer(RT ~ Congruence * Block * Response + Age + Sex + (1 | ID), data=df.morph.hs.noOut)
summary(model.morph.out.conf)
print(emmeans(model.morph.out.conf, pairwise ~ Congruence | Block))
contrast(emmeans(model.morph.out.conf, ~ Congruence * Block), interaction = "pairwise")





####################################################
#### Signal detection theory for detection task ####
####################################################


### dprime

# add a column Target for happy sad and morph primes and control for control primes
df.detection$Stim <- ifelse(df.detection$Prime %in% c("Happy", "Sad", "Morph"), "Target", "Control")
df.detection$Participant <- df.detection$ID

# define functions
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
seenID <- aggregate(df.detection, Seen ~ ID + Stim, mean)
# plot it
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


### 1. calculate visibility dprime 
###################################

subjdprime_visibility <- data.frame(compdprime(df.detection,'Seen', 1, 'Stim' , 'Target', 'Control')[1])
names(subjdprime_visibility) <- "dprime"
subjdprime_visibility['hit'] <- data.frame(compdprime(df.detection,'Seen', 1, 'Stim' , 'Target', 'Control')[2])
subjdprime_visibility['fa'] <- data.frame(compdprime(df.detection,'Seen', 1, 'Stim' , 'Target', 'Control')[3])
subjdprime_visibility$Participant <- as.factor(rownames(subjdprime_visibility))
subjdprime_visibility <- subjdprime_visibility[subjdprime_visibility$Participant %in% df.detection$ID,]
subjdprime_visibility$ID <- subjdprime_visibility$Participant

# mean dp
mean(subjdprime_visibility$dprime, na.rm=T)

# compute accuracy
accuracy <- subjdprime_visibility %>%
  summarise(
    mean_accuracy = mean(hit * 100, na.rm = TRUE)
  )

# t-test
t_test_dprime <- t.test(subjdprime_visibility$dprime, mu = 0)
print(t_test_dprime) 


### 2 : visibility happy versus sad prime
##########################################

# only keep happy and sad primes
df.dprime.HS <- df.detection[df.detection$Prime== 'Happy' | df.detection$Prime== 'Sad' ,]
df.dprime.HS <- df.dprime.HS[df.dprime.HS$Target== 'Happy' | df.dprime.HS$Target== 'Sad' ,]

# calculate d' for happy vs. sad prime perception
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
  summarise(mean_accuracy = mean(hit * 100, na.rm = TRUE)) 

# t-test
print(t.test(subjdprime_happySad$dprime, mu = 0)) 


# correlate subjdprime happy sad with priming effect
# check correlation between the priming effect and the prime visibility on masked trials
## happy
happy.corr <- df.happy.corr[df.happy.corr$Block=='subliminal',]
happy.corr <- happy.corr[happy.corr$Prime=='Happy' | happy.corr$Prime=='Sad',]
happy.cong <- happy.corr[, c("ID", "Prime", "RT", "Target","Block")]
happy.cong <- aggregate(happy.cong, RT ~ Prime + ID, mean)

happy.cong <- happy.cong %>%
  pivot_wider(
    names_from = Prime,    
    values_from = RT)

happy.cong$deltaH <- happy.cong$Sad - happy.cong$Happy

## sad targets
sad.corr <- df.sad.corr[df.sad.corr$Block=='subliminal',]
sad.cong <- sad.corr[, c("ID", "Prime", "RT", "Target","Block")]
sad.cong <- aggregate(sad.cong, RT ~ Prime + ID, mean)

sad.cong <- sad.cong %>%
  pivot_wider(
    names_from = Prime,    
    values_from = RT)

sad.cong$deltaS <- sad.cong$Happy - sad.cong$Sad


# combine happy and sad targets for priming effects
overall.cong <- merge(happy.cong,sad.cong, by='ID')
overall.cong$overall.priming <- rowMeans(overall.cong[, c("deltaH", "deltaS")])

# combine with dprime
overall.cong <- merge(overall.cong,subjdprime_happySad, by="ID",)

### Graph: 
ggplot(overall.cong, aes(x = dprime, y = overall.priming)) + 
  geom_point(color='#2f4858', alpha=0.8) + 
  geom_smooth(method = 'lm', alpha = 0.3, color ='#2f4858') +
  xlab("dprime Happy-Sad") + ylab("RT priming effects (ms)") + 
  ggtitle("Healthy participants") +
  theme_blank()+
  theme(plot.title = element_text(size = 16),
        axis.title.x = element_text(size = 16),   # x‐axis label
        axis.title.y = element_text(size = 16),   # y‐axis label
        axis.text.x  = element_text(size = 14),   # x‐axis tick labels
        axis.text.y  = element_text(size = 14))+
  
  labs(tag = "A") +
  theme(
    plot.tag.position = c(0.02, 0.98),  # x, y in [0,1], (0,1) = top-left
    plot.tag = element_text(size = 18, face = "bold"))+
  
  annotate("text", x = 3.8, y = 100, label = "p = .780", 
           size = 5, fontface = "italic")


### Stats:

# check normality before correlation
shapiro.test(overall.cong$overall.priming) 
shapiro.test(overall.cong$dprime) 

cor.test(overall.cong$overall.priming, overall.cong$dprime)

# check the significance of the intercept of this regression
model.cong <- lm(overall.priming ~ dprime, data = overall.cong)
residuals <- residuals(model.cong)
shapiro.test(residuals)
summary(model.cong)


