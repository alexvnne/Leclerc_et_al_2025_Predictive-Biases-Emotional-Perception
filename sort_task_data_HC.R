# Script to Sort Data Healthy Participants - AL
# Healthy participants
# v0.2 : AL 28.05.2025


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
demographics <- read.csv('demographics_file.csv',na.strings=c('',' ','NA')) # put demographics file here
# only keep data from participants who completed the study
demographics <- demographics[demographics$Completion.code=='CMCK1AHR' & demographics$Status=='APPROVED',]
# extract desired columns
extract <- c('Participant.id','Time.taken','Age','Sex','Ethnicity.simplified','Nationality','Student.status','Employment.status','Language')
demographics <- demographics[, extract]
# rename the column to ID
colnames(demographics)[colnames(demographics) == 'Participant.id'] <- 'ID'


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
# write priming task data csv
write.csv(df.full, '/data_priming_task_HC.csv')


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
# write the DDM csv
write.csv(ddm, '/data_DDM_HC.csv')


##################################
#### Sort detection task data ####
##################################

# extract df for detection task
df.detection <- df.all[!is.na(df.all$Valence),]
# only keep columns that we need
extract_detection <- c('ID','Valence','Seen','Prime','Target','prime_exp_7.tStart','mask_8.tStart','stim_file')
df.detection <- df.detection[,extract_detection]
# exclude from detection df outlier participants 
df.detection <- df.detection[!df.detection$ID %in% out.ppt$ID,]
length(unique(df.detection$ID)) # check how many participants there are
#write the detection csv
write.csv(df.detection, '/data_detection_task_HC.csv')


#########################
#### extract data for PCA
#########################

# make the df for the PCA
# #panas
df.panas <- df.all[!is.na(df.all$slider_panas.response),]
# only keep columns that we need
extract_panas <- c('slider_panas.response','ID','words')
df.panas <- df.panas[,extract_panas]
panas_pca <- pivot_wider(df.panas, names_from = words, values_from = slider_panas.response)
# delete excluded participants
panas_pca <- panas_pca[!(panas_pca$ID %in% out.ppt$ID),]

# #mathys
df.mathys1 <- df.all[!is.na(df.all$slider_mathys1.response),]
# only keep columns that we need
extract_mathys1 <- c('ID','slider_mathys1.response','item_type','direction','left_side','right_side')
df.mathys1 <- df.mathys1[,extract_mathys1]
# some mathys items are reversed, inverse scores
df.mathys1$slider_mathys1.response[df.mathys1$direction == "inverse"] <- 10 - df.mathys1$slider_mathys1.response[df.mathys1$direction == "inverse"]
df.mathys1 <- df.mathys1 %>% mutate(item = ifelse(direction == "normal", right_side, left_side))
extract_mathys1 <- c('ID','slider_mathys1.response','item')
df.mathys1 <- df.mathys1[,extract_mathys1]
mathys1_pca <- pivot_wider(df.mathys1, names_from = item, values_from = slider_mathys1.response)
# delete excluded participants
mathys1_pca <- mathys1_pca[!(mathys1_pca$ID %in% out.ppt$ID),]

df.mathys2 <- df.all[!is.na(df.all$slider_mathys2.response),]
# only keep columns that we need
extract_mathys2 <- c('ID','slider_mathys2.response','words_mathys')
df.mathys2 <- df.mathys2[,extract_mathys2]
mathys2_pca <- pivot_wider(df.mathys2, names_from = words_mathys, values_from = slider_mathys2.response)
# delete excluded participants
mathys2_pca <- mathys2_pca[!(mathys2_pca$ID %in% out.ppt$ID),]

# # merge them all
pca_all <- merge(panas_pca, mathys1_pca, by='ID')
pca_all <- merge(pca_all, mathys2_pca, by='ID')

colnames(pca_all)[colnames(pca_all) == 'ID'] <- 'participant_id'
pca_all$Group <- 'Healthy_Controls'
# write csv with this df to later run the pca
write.table(pca_all,"/pca_healthy_controls.csv")

