# Step 3: Interpolation PrC
library(tidyverse)
library(vegan)
library(ggplot2)
library(gratia)
library(tidypaleo) #remotes::install_github("paleolimbot/tidypaleo")
library(patchwork)
library(mvgam)

## Read in full PrC dataset
prc_pollen_diatoms_df <- read.csv("outputs/prc_pollen_diatoms_F7F9_df_February2025.csv")[-1]
head(prc_pollen_diatoms_df)

## PrC interpolation
### Interpolation to the coarser dataset: diatoms
data_to_int <- prc_pollen_diatoms_df %>%
  filter(proxy=="aquatic_pollen")

diatomsPrC <- prc_pollen_diatoms_df %>%
  filter(proxy=="diatoms")

PrC_i <- as.data.frame(approx(data_to_int$upper_age, data_to_int$PrC, diatomsPrC$upper_age)$y)
PrC_i$age <- diatomsPrC$upper_age
colnames(PrC_i) <- c("aquatic_pollen", "age")

prc_F7F9_df_i <- PrC_i
prc_F7F9_df_i$aquatic <- PrC_i$aquatic_pollen
prc_F7F9_df_i$diatoms <- diatomsPrC$PrC
prc_F7F9_df_i$lower_age <- diatomsPrC$lower_age

data_to_int <- prc_pollen_diatoms_df %>%
  filter(proxy=="terrestrial")
PrC_i <- as.data.frame(approx(data_to_int$upper_age, data_to_int$PrC, diatomsPrC$upper_age)$y)
PrC_i$age <- diatomsPrC$upper_age
colnames(PrC_i) <- c("terrestrial", "age")

prc_F7F9_df_i$terrestrial <- PrC_i$terrestrial

## Climate reconstructions interpolation to diatoms
### Interpolation to the coarser dataset: diatoms
climate_rec <- read.csv("data/Reconstructions_raw_reconstructions.csv", sep=",") %>%
  rename(age_bp=age..kyr.BP. ) %>%
  mutate(age_bp=age_bp*1000) %>%
  filter(age_bp<280000)

clim_i <- as.data.frame(approx(climate_rec$age_bp, climate_rec$mat, diatomsPrC$upper_age)$y)
clim_i$age <- diatomsPrC$upper_age
colnames(clim_i) <- c("mat_i", "age")
#clim_i$diat_prc <- diatomsPrC$PrC

prc_F7F9_df_i$mat_i <- clim_i$mat_i

clim_i <- as.data.frame(approx(climate_rec$age_bp, climate_rec$pann, diatomsPrC$upper_age)$y)
clim_i$age <- diatomsPrC$upper_age
colnames(clim_i) <- c("map_i", "age")

prc_F7F9_df_i$map_i <- clim_i$map_i

clim_i <- as.data.frame(approx(climate_rec$age_bp, climate_rec$ai, diatomsPrC$upper_age)$y)
clim_i$age <- diatomsPrC$upper_age
colnames(clim_i) <- c("ai_i", "age")

prc_F7F9_df_i$ai_i <- clim_i$ai_i
#write.csv(prc_F7F9_df_i, "outputs/F7_F9/prc_F7F9_df_i_April2025.csv")

prc_F7F9_df_i <- read.csv("outputs/F7_F9/prc_F7F9_df_i_April2025.csv")[-1]

# ## Read in the geochemistry dataset
# F9_geochem_grainsize <- read.csv("outputs/F9_geochem_grainsize_i.csv", sep = ",")[-1] 

geochem_i <- as.data.frame(approx(F9_geochem_grainsize$upper_age, F9_geochem_grainsize$TOC_percent, prc_F7F9_df_i$upper_age)$y)
#colnames(geochem_i) <- c("TOC_percent", "age")

geochem_i$TC_TN <- as.vector(approx(F9_geochem_grainsize$upper_age, F9_geochem_grainsize$TC_TN, prc_F7F9_df_i$upper_age)$y)
geochem_i$d13C <- as.vector(approx(F9_geochem_grainsize$upper_age, F9_geochem_grainsize$d13C, prc_F7F9_df_i$upper_age)$y)
geochem_i$LOI375 <- as.vector(approx(F9_geochem_grainsize$upper_age, F9_geochem_grainsize$LOI375, prc_F7F9_df_i$upper_age)$y)
geochem_i$sand <- as.vector(approx(F9_geochem_grainsize$upper_age, F9_geochem_grainsize$sand, prc_F7F9_df_i$upper_age)$y)
geochem_i$coarsesilt <- as.vector(approx(F9_geochem_grainsize$upper_age, F9_geochem_grainsize$coarsesilt, prc_F7F9_df_i$upper_age)$y)
geochem_i$finesilt <- as.vector(approx(F9_geochem_grainsize$upper_age, F9_geochem_grainsize$finesilt, prc_F7F9_df_i$upper_age)$y)
geochem_i$clay <- as.vector(approx(F9_geochem_grainsize$upper_age, F9_geochem_grainsize$clay, prc_F7F9_df_i$upper_age)$y)
geochem_i$C_percent <- as.vector(approx(F9_geochem_grainsize$upper_age, F9_geochem_grainsize$C_percent, prc_F7F9_df_i$upper_age)$y)
# geochem_i$upper_age <- prc_F7F9_df_i$upper_age
# geochem_i$lower_age <- prc_F7F9_df_i$lower_age

names(geochem_i) <- c("TOC_percent", "TC_TN","d13C","LOI375","sand","coarsesilt", "finesilt","clay","C_percent")
prc_geochem_F7F9_i <- cbind(prc_F7F9_df_i,geochem_i) %>%
  arrange(upper_age)
# 
# write.csv(prc_geochem_F7F9_i, "outputs/F7_F9/geochem_F9_df_i_April2025.csv")


## Relationships among PrCs
# Plot interpolated Pcurves with depth and ages
# Read in previous saved dataset
prc_F7F9_df_i <- read.csv("outputs/F7_F9/prc_geochem_F7F9_df_i_April2025.csv")[-1]

prc_df_long <- gather(prc_F7F9_df_i, proxy, PrC, -upper_age, -lower_age)

PrC_plot <- ggplot(prc_df_long, aes(x = upper_age, y = PrC)) +
  geom_line() + geom_point() +
  facet_wrap(~proxy, ncol=3,scales = "free_y")+
  #geom_smooth() +
  labs(y = "Principal curve interpolated", x = "Age (cal yr BP)") +
  scale_x_continuous(breaks=seq(0, 284000, by=40000)) +
  ggtitle("Fuquene F7&F9 aquatic and terrestrial trends") +
  geom_vline(xintercept = 24000, linetype=2) #F7 onset
PrC_plot

terr <- prc_F7F9_df_i %>% 
  ggplot(aes(x = terrestrial, y = diatoms, colour=upper_age)) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 12),
              col = 'darkred', fill = "#A25050") +
  labs(title = '',
       y = "Diatoms PrC", 
       x = 'Terrestrial Pollen PrC') 

aquatic <- prc_F7F9_df_i %>% 
  ggplot(aes(x = aquatic, y = diatoms,colour=upper_age)) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 12),
              col = 'darkred', fill = "#A25050") +
  labs(title = '',
       y = NULL, 
       x = 'Aquatic pollen PrC') + theme(legend.position = "none")

mat <- prc_F7F9_df_i %>% 
  ggplot(aes(x = mat_i, y = diatoms,colour=upper_age)) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 12),
              col = 'darkred', fill = "#A25050") +
  labs(title = '',
       y = NULL, 
       x = 'MAT') + theme(legend.position = "none") 

map <- prc_F7F9_df_i %>% 
  ggplot(aes(x = map_i, y = diatoms,colour=upper_age)) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 12),
              col = 'darkred', fill = "#A25050") +
  labs(title = '',
       y = NULL, 
       x = 'MAP') + theme(legend.position = "none")

ai <- prc_F7F9_df_i %>% 
  ggplot(aes(x = ai_i, y = diatoms,colour=upper_age)) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 12),
              col = 'darkred', fill = "#A25050") +
  labs(title = '',
       y = NULL, 
       x = 'Aridity index') + theme(legend.position = "none")

d13C <- prc_F7F9_df_i %>% 
  ggplot(aes(x = d13C, y = diatoms, colour=upper_age)) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 12),
              col = 'darkred', fill = "#A25050") +
  labs(title = '',
       y = "Diatoms PrC", 
       x = 'd13C') 

egg::ggarrange(terr,aquatic,mat,map,ai,nrow = 2)
