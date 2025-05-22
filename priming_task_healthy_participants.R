# Script for Priming Task - AL
# Healthy participants
# v0.2 : AL 27.05.2024


##########################
#### Load data frames ####
##########################

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

# load files
files <- list.files('data/')
csv_files <- grep ('.csv',files)
files <- files[csv_files]

# initialize an empty data.table to store all df
df.all <- data.table()

# loop through each file name in my file list
for (filename in paste0('data/',files)) {
  
  df <- fread(filename, na.strings=c('','NA'))
  
  columns_to_extract <- c('participant','frameRate','responseMain.keys','responseMain.rt',
                          'mask.tStart','block','prime_exp.tStart','supraliminal2.tStart',
                          'main_block.thisN','OS','stim_file','prime_emotion','target_emotion',
                          'slider_panas.response','slider_panas.rt','words','valence',
                          'slider_mathys1.response','slider_mathys1.rt','item_type','direction',
                          'slider_mathys2.response','slider_mathys2.rt','words_mathys',
                          'mathys3.text','order','response.keys','key_resp_2.keys',
                          'gender_target','target','left_side','right_side','mask_8.tStart',
                          'prime_exp_7.tStart')
  
  # convert keys if vb_positive is 1
  if (df$vb_positive[1] == 1) {
    df$responseMain.keys <- ifelse(df$responseMain.keys == 'v', 0, 1)
    df$response.keys <- ifelse(df$response.keys == 'v', 0, 1)
    df$key_resp_2.keys <- ifelse(df$key_resp_2.keys == 'v', 1, 0)
    
  } else {
    df$responseMain.keys <- ifelse(df$responseMain.keys == 'b', 0, 1)
    df$response.keys <- ifelse(df$response.keys == 'b', 0, 1)
    df$key_resp_2.keys <- ifelse(df$key_resp_2.keys == 'b', 1, 0)
  }
  
  # convert to numeric
  df$responseMain.keys <- as.numeric(df$responseMain.keys)
  df$response.keys <- as.numeric(df$response.keys)
  df$key_resp_2.keys <- as.numeric(df$key_resp_2.keys)
  
  # select desired columns
  selected_data <- df[, ..columns_to_extract]
  
  # combine selected data with df
  df.all <- rbind(df.all, selected_data)
}

# check how many participants we have
subject_counts <- table(df.all$participant)

# also load demographic data
demographics <- read.csv('prolific_export_65df4ddfffb481c878ff618a.csv',na.strings=c('',' ','NA')) # put demographics here
# only keep data from participants who completed the study
demographics <- demographics[demographics$Completion.code=='CMCK1AHR' & demographics$Status=='APPROVED',]
# extract desired columns
extract <- c('Participant.id','Time.taken','Age','Sex','Ethnicity.simplified','Nationality','Student.status','Employment.status','Language')
demographics <- demographics[, extract]
# rename the column to ID
colnames(demographics)[colnames(demographics) == 'Participant.id'] <- 'ID'

# load PCA data
pca <- read.table('path/to/PCA/data/PCA.csv')
# only keep data from healthy participants
pca <- pca[pca$Group=='HC',]
# rename the column to ID
colnames(pca)[colnames(pca) == 'participant_id'] <- 'ID'



########################
### Data definition ####
########################

# rename columns
colnames(df.all)[colnames(df.all) == 'participant'] <- 'ID'
colnames(df.all)[colnames(df.all) == 'responseMain.keys'] <- 'Response'
colnames(df.all)[colnames(df.all) == 'responseMain.rt'] <- 'RT'
colnames(df.all)[colnames(df.all) == 'target_emotion'] <- 'Target'
colnames(df.all)[colnames(df.all) == 'prime_emotion'] <- 'Prime'
colnames(df.all)[colnames(df.all) == 'main_block.thisN'] <- 'trial'
colnames(df.all)[colnames(df.all) == 'block'] <- 'Block'
colnames(df.all)[colnames(df.all) == 'response.keys'] <- 'Valence'
colnames(df.all)[colnames(df.all) == 'key_resp_2.keys'] <- 'Seen'

# rename values in columns
df.all$Target <- plyr::revalue(df.all$Target, c("sad"="Sad","happy"="Happy","morph"="Morph"))
df.all$Prime <- plyr::revalue(df.all$Prime, c("sad"="Sad","happy"="Happy","morph"="Morph","control"="Control"))

# define factors
df.all$ID <- factor(df.all$ID) 
df.all$Target <- factor(df.all$Target)
df.all$Prime <- factor(df.all$Prime)
df.all$Block <- factor(df.all$Block)
df.all$stim_file <- factor(df.all$stim_file)
demographics$ID <- factor(demographics$ID)
demographics$Sex <- factor(demographics$Sex)

# convert RT in ms
df.all$RT <- df.all$RT*1000

# add a correct column for happy and sad targets
df.all$correct.happy <- ifelse(df.all$Target=='Happy' & df.all$Response == 1, 1, 0)
df.all$correct.sad <- ifelse(df.all$Target=='Sad' & df.all$Response == 0, 1, 0)



##################################################
### Exclude participants with technical issues ###
##################################################

# get target data
df.target <- df.all[!is.na(df.all$Response),]
# get unique id per line to get rid of trials later
df.target <- df.target %>%
  mutate(row_extract = row_number())

# calculate prime presentation for all trials
df.target$primeSubli <- df.target$mask.tStart - df.target$prime_exp.tStart
df.target$primeSupra <- df.target$supraliminal2.tStart - df.target$prime_exp.tStart
# create two different df for subli and supra
df.subli <- df.target[df.target$Block == 'subliminal',]
df.supra <- df.target[df.target$Block == 'supraliminal',]

# check for trials that were below 20ms and above 45ms
prime.out.subli <- df.subli[df.subli$primeSubli <= 0.020 | df.subli$primeSubli >= 0.045,]
prime.out.supra <- df.supra[df.supra$primeSupra <= 0.020 | df.supra$primeSupra >= 0.045,]
# count per participant
prime.out.subli$ID <- as.character(prime.out.subli$ID) # in subli
count <- table(prime.out.subli$ID)
count.fr.subli <- as.data.frame(count)
prime.out.supra$ID <- as.character(prime.out.supra$ID) # in supra
count <- table(prime.out.supra$ID)
count.fr.supra <- as.data.frame(count)

# exlude participants who have more than 5% of their trial below 20ms and above 45ms
x <- 0.15 * 432 
out.ppt.subli <- count.fr.subli[count.fr.subli$Freq > x,]
out.ppt.supra <- count.fr.supra[count.fr.supra$Freq > x,]
colnames(out.ppt.subli)[colnames(out.ppt.subli) == 'Var1'] <- 'ID' # rename columns
colnames(out.ppt.supra)[colnames(out.ppt.supra) == 'Var1'] <- 'ID' 

out.ppt <- merge(out.ppt.supra, out.ppt.subli, by='ID', all=T)

#### exclude participants who aberrantly rated happy and sad targets

# first define each df per target
df.happy <- df.target[df.target$Target == 'Happy',]
df.sad <- df.target[df.target$Target == 'Sad',]
df.morph <- df.target[df.target$Target == 'Morph',]

# get accuracy scores for happy targets per participant
df.happy.acc <- aggregate(Response ~ID, data=df.happy, mean) # happy
df.happy.acc$Z.Response <- scale(df.happy.acc$Response)
HR.out <- df.happy.acc[df.happy.acc$Z.Response > 2 | df.happy.acc$Z.Response < (-2),]

# get accuracy scores for sad targets per participant
df.sad.acc <- aggregate(Response ~ID, data=df.sad, mean) # sad
df.sad.acc$Z.Response <- scale(df.sad.acc$Response)
SR.out <- df.sad.acc[df.sad.acc$Z.Response > 2 | df.sad.acc$Z.Response < (-2),]

# merge them (by ID and only show the ones that were both outliers in happy and sad responses)
outliers.acc <- merge(HR.out, SR.out, by='ID',all=T)

out.ppt <- merge(out.ppt,outliers.acc, by='ID', all=T)

df.target <- df.target[!df.target$ID %in% out.ppt$ID,]


# check aberrant RT per targets per participant
df.hs <- df.target[!df.target$Target=='Morph',]
df.hs.out <- aggregate(RT ~ID, data=df.hs, mean)
df.hs.out$Z.RT <- scale(df.hs.out$RT)
hs.out <- df.hs.out[df.hs.out$Z.RT > 2 | df.hs.out$Z.RT < (-2),]

out.ppt <- merge(hs.out, out.ppt, by='ID', all=T)

# delete outliers from the df
df.target <- df.target[!df.target$ID %in% out.ppt$ID,]


### outliers rt per trial per target
####################################

# exclude trials with RT > 3s (meaning technical problem)
df.target <- df.target[!df.target$RT > 3000,]

# define 3 different df for rt trimming
df.happy <- df.target[df.target$Target == 'Happy',]
df.sad <- df.target[df.target$Target == 'Sad',]
df.morph <- df.target[df.target$Target == 'Morph',]

# plot RT distribution
ggplot(df.target) + aes(x=RT) + geom_histogram() + facet_grid(df.target$Target)+
  labs(x='Reaction time',
       y='Count')

df.sad$Z.RT <- scale(df.sad$RT)
trials.sad.out <- df.sad[df.sad$Z.RT > 3 | df.sad$Z.RT < (-3),]

df.happy$Z.RT <- scale(df.happy$RT)
trials.happy.out <- df.happy[df.happy$Z.RT > 3 | df.happy$Z.RT < (-3),]

df.morph$Z.RT <- scale(df.morph$RT)
trials.morph.out <- df.morph[df.morph$Z.RT > 3 | df.morph$Z.RT < (-3),]


# group the 3 dataframes together
trials.out <- merge(trials.sad.out,trials.happy.out, all=T)
trials.out <- merge(trials.out,trials.morph.out, all=T)

# exclude aberrant trials from the analysis
df.target <- df.target[!df.target$row_extract %in% trials.out$row_extract,]

# replot histogram
ggplot(df.target) + aes(x=RT) + geom_histogram() + facet_grid(df.target$Target)+
  labs(x='Reaction time',
       y='Count')




#################################
########## Demographics #########
#################################

# convert age to numerical data
demographics$Age <- as.numeric(demographics$Age)
# convert ethnicity to factor
demographics$Ethnicity.simplified <- factor(demographics$Ethnicity.simplified)

# get a df for demographics (age and sex)
desc_stats <- demographics %>%
  summarise(
    mean_age = mean(Age, na.rm = T),
    sd_age = sd(Age, na.rm = T),
    mean_time = mean(Time.taken, na.rm=T),
    sd_time = sd(Time.taken, na.rm=T)
  )

# exclude outlier participants from demographics df
demographics <- demographics[!demographics$ID %in% out.ppt$ID,]

# check that demographics and df.targets have the same number of ppts
length(unique(df.target$ID))
length(unique(demographics$ID))

# merge demographic data to clean df
df.full <- merge(df.target, demographics, by='ID')
# merge it all with PCA scores
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


################################################
##### Extract a data frame for DDM analysis ####
################################################

# get morph df because DDM is only fitted on morphed targets data
ddm <- df.morph
# only keep columns that are relevant for DDM
extract_ddm <- c('ID','Response','RT','Block','Prime','PC1')
ddm <- ddm[,extract_ddm]

# reconvert ddm RT into seconds
ddm$RT <- ddm$RT/1000
write.csv(ddm, 'path/to/write/csv/for/DDM.csv')


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
  geom_boxplot(alpha = 0.55, width = 0.4, outlier.shape = NA, 
               position = position_dodge(width = 0.25)) +
  geom_jitter(aes(color = Prime), alpha = 0.6, size = 1.4, width = 0.15) +
  facet_wrap(~ Block, strip.position = 'bottom') +  
  stat_summary(fun = mean, shape = 23, size = 0.8, stroke = 0.8, color = "black", alpha = 0.8) +  
  theme_blank() +  
  scale_color_manual(values = lPalette) +
  scale_fill_manual(values = lPalette) +
  ggtitle("Healthy participants") +
  theme(plot.title = element_text(size = 16),
        axis.title.x = element_text(size = 16),   # x‐axis label
        axis.title.y = element_text(size = 16),   # y‐axis label
        axis.text.x  = element_text(size = 14),   # x‐axis tick labels
        axis.text.y  = element_text(size = 14),
        strip.text.x = element_text(size = 14)) +
  
  theme(legend.position = "none")+
  geom_hline(yintercept = 0.5,
             linetype   = "dashed",
             color      = "grey40",
             size       = 0.45) +
  
  labs(tag = "A") +
  theme(
    plot.tag.position = c(0.02, 0.98),  # x, y in [0,1], (0,1) = top-left
    plot.tag = element_text(size = 18, face = "bold"))+
  
  # Masked
  geom_signif(data = subset(agg.morph, Block == "Masked"), 
              comparisons = list(c("Happy", "Control")), 
              annotations = c("**"), 
              y_position = 1.05, 
              tip_length = 0.02, 
              textsize = 4, 
              vjust = 0.4, 
              color = "black") +
  
  # Unmasked
  geom_signif(data = subset(agg.morph, Block == "Unmasked"), 
              comparisons = list(c("Happy", "Control")), 
              annotations = c("***"), 
              y_position = 1.05, 
              tip_length = 0.02, 
              textsize = 4, 
              vjust = 0.4, 
              color = "black") +
  
  # Unmasked
  geom_signif(data = subset(agg.morph, Block == "Unmasked"), 
              comparisons = list(c("Sad", "Control")), 
              annotations = c("**"), 
              y_position = 1.15, 
              tip_length = 0.02, 
              textsize = 4, 
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
  theme(plot.title = element_text(size = 16),
        axis.title.x = element_text(size = 16),   # x‐axis label
        axis.title.y = element_text(size = 16),   # y‐axis label
        axis.text.x  = element_text(size = 14),   # x‐axis tick labels
        axis.text.y  = element_text(size = 14))+
  
  labs(tag = "A") +
  theme(
    plot.tag.position = c(0.02, 0.98),  # x, y in [0,1], (0,1) = top-left
    plot.tag = element_text(size = 18, face = "bold"))+
  
  annotate("text", x = 4, y = 0.9, label = "p = .025*", 
           size = 5, fontface = "italic")

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
    plot.title = element_text(size = 16),
    axis.title.x = element_text(size = 16),   # x‐axis label
    axis.title.y = element_text(size = 16),   # y‐axis label
    axis.text.x  = element_text(size = 14),   # x‐axis tick labels
    axis.text.y  = element_text(size = 14),
    strip.text.x = element_text(size = 14)
  ) +
  scale_color_manual(values = rtPalette) +
  scale_fill_manual(values = rtPalette) +
  
  labs(tag = "A") +
  theme(
    plot.tag.position = c(0.02, 0.98),  # x, y in [0,1], (0,1) = top-left
    plot.tag = element_text(size = 18, face = "bold"))+
  
  # Add significance bars for each condition separately
  # Masked
  geom_signif(data = subset(agg.df.hs, Block == "Masked"), 
              comparisons = list(c("Congruent", "Incongruent")), 
              annotations = c("***"),  # Adjust p-value annotation
              y_position = 820,        # Adjust height of significance bar
              tip_length = 0.02, 
              textsize = 4, 
              vjust = 0.4, 
              color = "black")+
  
  # Unmasked
  geom_signif(data = subset(agg.df.hs, Block == "Unmasked"), 
              comparisons = list(c("Congruent", "Incongruent")), 
              annotations = c("***"),  
              y_position = 820,        
              tip_length = 0.02, 
              textsize = 4, 
              vjust = 0.4, 
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
  facet_wrap( ~ Block)+
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, color = "black", fatten = 2) +  
  geom_jitter(width = 0.2, alpha =1)+
  theme_blank()+
  theme(
    legend.position = "none",
    plot.title = element_text(size = 16),
    axis.title.x = element_text(size = 16),   # x‐axis label
    axis.title.y = element_text(size = 16),   # y‐axis label
    axis.text.x  = element_text(size = 14),   # x‐axis tick labels
    axis.text.y  = element_text(size = 14),
    strip.text.x = element_text(size = 14)
  ) +
  scale_color_manual(values=rtPalette)+
  scale_fill_manual(values = rtPalette)+
  
  labs(tag = "A") +
  theme(
    plot.tag.position = c(0.02, 0.98),  # x, y in [0,1], (0,1) = top-left
    plot.tag = element_text(size = 18, face = "bold"))+
  
  # Add significance bars for each condition separately
  # Masked
  geom_signif(data = subset(agg.morph, Block == "Masked"), 
              comparisons = list(c("Congruent", "Incongruent")), 
              annotations = c("***"),  
              y_position = 1300,        
              tip_length = 0.02, 
              textsize = 4, 
              vjust = 0.4, 
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
  labs(x = "Prime", y = "RT (ms)", title = 'Healthy participants')+
  facet_wrap( ~Block, strip.position = 'bottom')+
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, color = "black", fatten = 2) +  
  geom_jitter(width = 0.2, alpha =1)+
  theme_blank()+
  theme(
    legend.position = "none",
    plot.title = element_text(size = 16),
    axis.title.x = element_text(size = 16),   # x‐axis label
    axis.title.y = element_text(size = 16),   # y‐axis label
    axis.text.x  = element_text(size = 14),   # x‐axis tick labels
    axis.text.y  = element_text(size = 14),
    strip.text.x = element_text(size = 14)
  ) +
  #plot.title = element_text(hjust = 0.5, size = 12), legend.position = "none") +
  scale_color_manual(values=rtPalette)+
  scale_fill_manual(values = rtPalette)+
  coord_cartesian(ylim = c(430, 1350))+
  
  labs(tag = "A") +
  theme(
    plot.tag.position = c(0.02, 0.98),  # x, y in [0,1], (0,1) = top-left
    plot.tag = element_text(size = 18, face = "bold"))+
  
  # Add significance bars for each condition separately
  # Masked
  geom_signif(data = subset(agg.morph, Block == "Masked"), 
              comparisons = list(c("Congruent", "Incongruent")), 
              annotations = c("***"),  
              y_position = 1050,        
              tip_length = 0.02, 
              textsize = 4, 
              vjust = 0.4, 
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

