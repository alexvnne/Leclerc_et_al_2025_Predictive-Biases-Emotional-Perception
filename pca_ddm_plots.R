# Script to plot PCA and DDM results - AL
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




#############
#### PCA ####
#############

# load data frame
pca.grouped <- read.table("/pca_score_HC_patients.csv")
# rename MDD to DEP
pca.grouped$Group <- plyr::revalue(pca.grouped$Group, c("MDD"="DEP"))
# relevel
pca.grouped <-  pca.grouped %>%
  filter(Group %in% c("HC","DEP","BD")) %>%
  mutate(Group = factor(Group, levels = c("HC","DEP","BD")))


# define palette - one color per group
disPalette <- c(
  HC  = "#613F75",  
  MDD = "#2f4858",  
  BD  = "#86bbd8"   
)

## Group differences on PCA scores
# PC1
# Graph: 
ggplot(pca.grouped, aes(y = PC1, x = Group, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.6) +      
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.3,
    color = "black",
    fatten = 2
  ) +                                           
  geom_jitter(
    aes(color = Group),
    width = 0.1,
    height = 0.1,
    size = 1.5,
    alpha = 0.8
  ) +                                           
  labs(
    y = "PC1 score",
    x = "Group",
    title = "Group differences on PC1 scores"
  ) +
  theme_blank() +
  theme(legend.position = "none")+
  scale_fill_manual(values=disPalette)+
  scale_color_manual(values=disPalette)+
  theme(plot.title = element_text(size = 28, face = 'bold',margin = margin(b=30)),
        axis.title.x = element_text(size = 26, color='black'),   # x‐axis label
        axis.title.y = element_text(size = 26, color='black'),   # y‐axis label
        axis.text.x  = element_text(size = 26, color='black'),   # x‐axis tick labels
        axis.text.y  = element_text(size = 26, color='black'),
        legend.text = element_text(size = 26, color='black'),
        legend.title = element_text(size = 26, color='black'))+
  coord_cartesian(ylim=c(-15,18))+
  
  geom_signif(
    comparisons = list(c("DEP", "HC")),
    annotations = "p = .012*",
    y_position = 15,    
    tip_length = 0.05,
    textsize = 7,
    size = 0.8)


# Stats:
summary(lm(PC1 ~ Group, data = pca.grouped))


# PC2
# Graph:
ggplot(pca.grouped, aes(y = PC2, x = Group, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.6) +      
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.3,
    color = "black",
    fatten = 2
  ) +                                          
  geom_jitter(
    aes(color = Group),
    width = 0.1,
    height = 0.1,
    size = 1.5,
    alpha = 0.8
  ) +                                           
  labs(
    y = "PC2 score",
    x = "Group",
    title = "Group differences on PC2 scores"
  ) +
  theme_blank() +
  theme(legend.position = "none")+
  scale_fill_manual(values=disPalette)+
  scale_color_manual(values=disPalette)+
  theme(plot.title = element_text(size = 26, face='bold'),
        axis.title.x = element_text(size = 24),   # x‐axis label
        axis.title.y = element_text(size = 24),   # y‐axis label
        axis.text.x  = element_text(size = 18),   # x‐axis tick labels
        axis.text.y  = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 22))+
  
  geom_signif(
    comparisons = list(c("DEP", "HC")),
    annotations = "p = .013*",
    y_position = 12,    # adjust this based on data range
    tip_length = 0.05,
    textsize = 6
  ) +
  
  geom_signif(
    comparisons = list(c("BD", "HC")),
    annotations = "p < .001***",
    y_position = 14.5,    # adjust this based on data range
    tip_length = 0.05,
    textsize = 6
  )

# Stats:
summary(lm(PC2 ~ Group, data = pca.grouped))



# PC3
# inverse PC3 scores
pca.grouped$PC3_inv <- - pca.grouped$PC3

# Graph:
ggplot(pca.grouped, aes(y = PC3_inv, x = Group, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.6) +      
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.3,
    color = "black",
    fatten = 2
  ) +                                          
  geom_jitter(
    aes(color = Group),
    width = 0.1,
    height = 0.1,
    size = 1.5,
    alpha = 0.8
  ) +                                           
  labs(
    y = "PC3 score",
    x = "Group",
    title = "Group differences on PC3 scores"
  ) +
  theme_blank() +
  theme(legend.position = "none")+
  scale_fill_manual(values=disPalette)+
  scale_color_manual(values=disPalette)+
  theme(plot.title = element_text(size = 26, face='bold'),
        axis.title.x = element_text(size = 24),   # x‐axis label
        axis.title.y = element_text(size = 24),   # y‐axis label
        axis.text.x  = element_text(size = 18),   # x‐axis tick labels
        axis.text.y  = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 22))


# Stats:
summary(lm(PC3 ~ Group, data = pca.grouped))



# PC4
# Graph:
gplot(pca.grouped, aes(y = PC4, x = Group, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.6) +      
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.3,
    color = "black",
    fatten = 2
  ) +                                          
  geom_jitter(
    aes(color = Group),
    width = 0.1,
    height = 0.1,
    size = 1.5,
    alpha = 0.8
  ) +                                           
  labs(
    y = "PC4 score",
    x = "Group",
    title = "Group differences on PC4 scores"
  ) +
  theme_blank() +
  theme(legend.position = "none")+
  scale_fill_manual(values=disPalette)+
  scale_color_manual(values=disPalette)+
  theme(plot.title = element_text(size = 26, face='bold'),
        axis.title.x = element_text(size = 24),   # x‐axis label
        axis.title.y = element_text(size = 24),   # y‐axis label
        axis.text.x  = element_text(size = 18),   # x‐axis tick labels
        axis.text.y  = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 22))

# Stats:
summary(lm(PC4 ~ Group, data = pca.grouped))




# Item contribution to PC loadings
# load data frame
diffrep <- read.csv("/PCA.csv")

# define palette
pal <- c("#80B192","#f26419",'#f6ae2d')

# reorder: PANAS first
diffrep$Scale <- relevel(as.factor(diffrep$Scale), ref = "PANAS")

# PC1

# sort from higher to lower
diffrep.pc1 <- diffrep %>%
  mutate(Feature = reorder(Feature, PC1)) %>%
  arrange(desc(PC1))

# threshold - only show items with feature contribution below -0.15 and above 0.15
diffrep.pc1 <- diffrep.pc1[diffrep.pc1$PC1 > 0.15 | diffrep.pc1$PC1 < -0.15,]

# plot
ggplot(diffrep.pc1) +
  aes(x = PC1, y = Feature, fill = Scale) +
  geom_bar(stat = "summary", alpha=0.7, width = 0.7) +
  labs(x='PC1 scores', y="Items")+
  scale_fill_manual(values=pal)+
  theme_blank()+
  ggtitle("Items loadings on PC1") +
  theme(plot.title = element_text(size = 41, face='bold'),
        axis.title.x = element_text(size = 39,color='black'),   # x‐axis label
        axis.title.y = element_text(size = 39,color='black'),   # y‐axis label
        axis.text.x  = element_text(size = 36,color='black'),   # x‐axis tick labels
        axis.text.y  = element_text(size = 28, color='#393E46'),
        legend.text = element_text(size = 38,color='black'),
        legend.title = element_text(size = 39,color='black')) 





# PC2

# sort from higher to lower
diffrep.pc2 <- diffrep %>%
  mutate(Feature = reorder(Feature, PC2)) %>%
  arrange(desc(PC2))

# threshold - only show items with feature contribution below -0.15 and above 0.15
diffrep.pc2 <- diffrep.pc2[diffrep.pc2$PC2 > 0.15 | diffrep.pc2$PC2 < -0.15,]

# plot
ggplot(diffrep.pc2) +
  aes(x = PC2, y = Feature, fill = Scale) +
  geom_bar(stat = "summary", alpha=0.7, width = 0.7) +
  labs(x='PC2 scores', y="Items")+
  scale_fill_manual(values=pal)+
  theme_blank()+
  ggtitle("Items loadings on PC2") +
  theme(plot.title = element_text(size = 26, face='bold'),
        axis.title.x = element_text(size = 24),   # x‐axis label
        axis.title.y = element_text(size = 24),   # y‐axis label
        axis.text.x  = element_text(size = 18),   # x‐axis tick labels
        axis.text.y  = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 22))



# PC3
# inverse PC3
diffrep.pc3$PC3_inv <- - diffrep.pc3$PC3

# threshold - only show items with feature contribution below -0.15 and above 0.15
diffrep.pc3 <- diffrep.pc3[diffrep.pc3$PC3_inv > 0.15 | diffrep.pc3$PC3_inv < -0.15,]

# sort from higher to lower
diffrep.pc3 <- diffrep.pc3 %>%
  mutate(Feature = reorder(Feature, PC3_inv)) %>%
  arrange(desc(PC3_inv))

# plot
ggplot(diffrep.pc3) +
  aes(x = PC3_inv, y = Feature, fill = Scale) +
  geom_bar(stat = "summary", alpha=0.7, width = 0.5) +
  labs(x='PC3 scores', y="Items")+
  scale_fill_manual(values=pal)+
  theme_blank()+
  ggtitle("Items loadings on PC3") +
  theme(plot.title = element_text(size = 26, face='bold'),
        axis.title.x = element_text(size = 24),   # x‐axis label
        axis.title.y = element_text(size = 24),   # y‐axis label
        axis.text.x  = element_text(size = 18),   # x‐axis tick labels
        axis.text.y  = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 22))


# PC4

# sort from higher to lower
diffrep.pc4 <- diffrep %>%
  mutate(Feature = reorder(Feature, PC4)) %>%
  arrange(desc(PC4))

# threshold - only show items with feature contribution below -0.15 and above 0.15
diffrep.pc4 <- diffrep.pc4[diffrep.pc4$PC4 > 0.15 | diffrep.pc4$PC4 < -0.15,]

# plot
ggplot(diffrep.pc4) +
  aes(x = PC4, y = Feature, fill = Scale) +
  geom_bar(stat = "summary", alpha=0.7, width = 0.7) +
  labs(x='PC4 scores', y="Items")+
  scale_fill_manual(values=pal)+
  theme_blank()+
  ggtitle("Items loadings on PC4") +
  theme(plot.title = element_text(size = 26, face='bold'),
        axis.title.x = element_text(size = 24),   # x‐axis label
        axis.title.y = element_text(size = 24),   # y‐axis label
        axis.text.x  = element_text(size = 18),   # x‐axis tick labels
        axis.text.y  = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 22))




#############
#### DDM ####
#############

# load data frame healthy controls
cond_effects_bias <- read.csv("/bias_HC.csv")
cond_effects_bias$block <- plyr::revalue(cond_effects_bias$block, c("subliminal"="Masked","supraliminal"="Unmasked"))

# plot effect of primes on starting point
ggplot(cond_effects_bias) +
  aes(y=estimate__, x=prime_emotion, fill = prime_emotion)+
  geom_violin(alpha=0.6)+
  stat_summary(fun = mean, geom = "crossbar", width = 0.3, color = "black") +
  geom_jitter(aes(color = prime_emotion), width = 0.2, alpha = 0.7,shape = 21, stroke = 0.05, color='black') +
  facet_wrap(~block, strip.position='bottom')+
  theme_blank()+
  scale_fill_manual(values = c("Control" = "#7F7F7F", "Happy" = "#F4A300", "Sad" = "#E87D3E")) +
  scale_color_manual(values = c("Control" = "#7F7F7F", "Happy" = "#F4A300", "Sad" = "#E87D3E")) +
  guides(fill = "none", color = "none")+
  coord_cartesian(ylim = c(0.35, 0.80)) +
  labs(
    x = "Prime",  
    y = "Estimated Bias",
    title = "Healthy participants")+
  theme(plot.title = element_text(size = 21),
        axis.title.x = element_text(size = 20),   # x‐axis label
        axis.title.y = element_text(size = 20),   # y‐axis label
        axis.text.x  = element_text(size = 18),   # x‐axis tick labels
        axis.text.y  = element_text(size = 18),
        strip.text.x = element_text(size = 18)) +
  
  labs(tag = "A") +
  theme(
    plot.tag.position = c(0.02, 0.98),  # x, y in [0,1], (0,1) = top-left
    plot.tag = element_text(size = 21, face = "bold"))+
  
  # Add significance bars for each condition separately
  # Masked HS
  geom_signif(data = subset(cond_effects_bias, block == "Masked"), 
              comparisons = list(c("Control", "Happy")), 
              annotations = c("BF = 10.19"),  
              y_position = 0.75,       
              tip_length = 0.02, 
              textsize = 5, 
              #vjust = 0.05, 
              color = "black")+
  
  # Masked HC
  geom_signif(data = subset(cond_effects_bias, block == "Masked"), 
              comparisons = list(c("Happy", "Sad")), 
              annotations = c("BF = 72204.20"),  
              y_position = 0.70,        
              tip_length = 0.02, 
              textsize = 5, 
              #vjust = 0.05, 
              color = "black")



# load data frame patients
cond_effects_bias <- read.csv("bias_patients.csv")


# plot BD
cond_effects_bias_bd <- cond_effects_bias[cond_effects_bias$group=='BD',]

ggplot(cond_effects_bias_bd) +
  aes(y=estimate__, x=prime_emotion, fill = prime_emotion)+
  #geom_boxplot(alpha=0.6)+
  geom_violin(alpha=0.6)+
  stat_summary(fun = mean, geom = "crossbar", width = 0.3, color = "black") +
  geom_jitter(aes(color = prime_emotion), width = 0.2, alpha = 0.7,shape = 21, stroke = 0.1, color='black') +
  facet_wrap(~block, strip.position = 'bottom')+
  theme_blank()+
  scale_fill_manual(values = c("Control" = "#7F7F7F", "Happy" = "#F4A300", "Sad" = "#E87D3E")) +
  scale_color_manual(values = c("Control" = "#7F7F7F", "Happy" = "#F4A300", "Sad" = "#E87D3E")) +
  coord_cartesian(ylim = c(0.3, 0.85)) +
  guides(fill = "none", color = "none")+
  labs(
    x = "Prime",  
    y = "Estimated Bias",
    title = "BD patients")+
  theme(plot.title = element_text(size = 21),
        axis.title.x = element_text(size = 20),   # x‐axis label
        axis.title.y = element_text(size = 20),   # y‐axis label
        axis.text.x  = element_text(size = 18),   # x‐axis tick labels
        axis.text.y  = element_text(size = 18),
        strip.text.x = element_text(size = 18)) +
  
  labs(tag = "B") +
  theme(
    plot.tag.position = c(0.02, 0.98),  # x, y in [0,1], (0,1) = top-left
    plot.tag = element_text(size = 21, face = "bold"))+
  
  # Add significance bars for each condition separately
  # BD Masked HS
  geom_signif(data = subset(cond_effects_bias, group == "BD" & block == "Masked"), 
              comparisons = list(c("Sad", "Happy")), 
              annotations = c("BF = 3535.40"),  
              y_position = 0.68,       
              tip_length = 0.02, 
              textsize = 5, 
              #vjust = 0.02, 
              color = "black")+
  
  # BD Masked HC
  geom_signif(data = subset(cond_effects_bias, group == "BD" & block == "Masked"), 
              comparisons = list(c("Happy", "Control")), 
              annotations = c("BF = 8.02"),  
              y_position = 0.82,        
              tip_length = 0.02, 
              textsize = 5, 
              #vjust = 0.01, 
              color = "black")+
  
  # BD Masked HC
  geom_signif(data = subset(cond_effects_bias, group == "BD" & block == "Masked"), 
              comparisons = list(c("Sad", "Control")), 
              annotations = c("BF = 3.96"),  
              y_position = 0.75,        
              tip_length = 0.02, 
              textsize = 5, 
              #vjust = 0.01, 
              color = "black")+
  
  # BD Masked HC
  geom_signif(data = subset(cond_effects_bias, group == "BD" & block == "Unmasked"), 
              comparisons = list(c("Sad", "Happy")), 
              annotations = c("BF = 4.27"),  
              y_position = 0.68,        
              tip_length = 0.02, 
              textsize = 5, 
              #vjust = 0.01, 
              color = "black")


# plot DEP
cond_effects_bias_mdd <- cond_effects_bias[cond_effects_bias$group=='MDD',]

ggplot(cond_effects_bias_mdd) +
  aes(y=estimate__, x=prime_emotion, fill = prime_emotion)+
  #geom_boxplot(alpha=0.6)+
  geom_violin(alpha=0.6)+
  stat_summary(fun = mean, geom = "crossbar", width = 0.3, color = "black") +
  geom_jitter(aes(color = prime_emotion), width = 0.2, alpha = 0.7,shape = 21, stroke = 0.1, color='black') +
  facet_wrap(~block, strip.position = 'bottom')+
  theme_blank()+
  scale_fill_manual(values = c("Control" = "#7F7F7F", "Happy" = "#F4A300", "Sad" = "#E87D3E")) +
  scale_color_manual(values = c("Control" = "#7F7F7F", "Happy" = "#F4A300", "Sad" = "#E87D3E")) +
  guides(fill = "none", color = "none")+
  labs(
    x = "Prime",  
    y = "Estimated Bias",
    title = "DEP patients")+
  theme(plot.title = element_text(size = 21),
        axis.title.x = element_text(size = 20),   # x‐axis label
        axis.title.y = element_text(size = 20),   # y‐axis label
        axis.text.x  = element_text(size = 18),   # x‐axis tick labels
        axis.text.y  = element_text(size = 18),
        strip.text.x = element_text(size = 18)) +
  
  labs(tag = "C") +
  theme(
    plot.tag.position = c(0.02, 0.98),  # x, y in [0,1], (0,1) = top-left
    plot.tag = element_text(size = 21, face = "bold"))


