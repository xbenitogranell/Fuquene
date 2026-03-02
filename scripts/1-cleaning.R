
## Step 1: Clean data and estimate Principal Curves (PrC) of the diatom composite core (M1, F7 and F9)
# PrC is a nonlinear ordination technique that extracts a single gradient of variation from multivariate data.

library(tidyverse)
library(vegan)
library(ggplot2)
library(gratia)
library(tidypaleo) #remotes::install_github("paleolimbot/tidypaleo")
library(patchwork)
library(mvgam)

F_diatoms <- read.csv("data/Fuquene_diatoms_counts.csv") %>%
  filter(!core=="M1") %>%
  drop_na(upper_age)

nms_diat <- read.csv("data/nms_changes_diatoms_M1_F7_F9.csv", sep=";")
#rename(core=ï..core) #rename extraneous column name 

agedepth <- F_diatoms[, names(F_diatoms) %in% c("depth", "upper_age", "lower_age", "core")] 
diat <- F_diatoms[, !names(F_diatoms) %in% c("depth", "upper_age", "lower_age", "core")]
diat[is.na(diat)] <- 0

diatoms_save <- cbind(agedepth, diat)
colnames(diatoms_save[,5:ncol(diatoms_save)]) -> nms_diat[,1]

#this is to transform to tidy format, calculate % and subset more common species
new <- diatoms_save %>% 
  gather(key = taxa, value = count, -depth, -upper_age, -lower_age, -core) %>%
  mutate(taxa = plyr::mapvalues(taxa, from = nms_diat[,1], to = nms_diat[,3])) %>%
  mutate(taxa = plyr::mapvalues(taxa, from = nms_diat[,3], to=nms_diat[,4])) %>%
  group_by(depth, taxa, upper_age, lower_age) %>%
  summarise(count = sum(count)) %>%
  filter(!count == "0") %>% #this is to remove empty samples (rows)
  filter(!upper_age==0.0) %>% #this is to drop extraneous ages
  filter(!is.na(lower_age)) %>%
  ungroup() %>%
  group_by(depth, upper_age, lower_age) %>%
  mutate(relative_abundance_percent = count / sum(count) * 100) %>%
  mutate(total_sample=sum(count)) %>%
  filter(total_sample>100) %>%
  mutate(elapsedTime = abs(upper_age - lower_age)) %>%
  mutate(negAge=-upper_age) %>%
  ungroup() 

unique(new$taxa)

# filter more abundant taxa; format need to be on long-wide format-->no spreaded 
core_common_taxa <- new %>%
  group_by(taxa) %>%
  summarise(max_rel_abund = max(relative_abundance_percent)) %>%
  filter(max_rel_abund >= 10) %>%
  arrange(max_rel_abund) %>%
  pull(taxa)

# select from initial table
core_counts_common <- new %>%
  filter(taxa %in% core_common_taxa) %>%
  mutate(taxa = factor(taxa, levels = core_common_taxa)) %>%
  arrange(taxa)

#mutate(n_valves=ifelse(total_sample < 100, "low", "high"))

#make it wide with RA
core_ra_wide_diatoms <- core_counts_common %>%
  select(depth, upper_age, lower_age, taxa, relative_abundance_percent) %>%
  spread(key = taxa, value = relative_abundance_percent) %>%
  arrange(upper_age) #sort by increasing time

levels(core_counts_common$taxa)

#make it wide with counts
core_counts_wide_diatoms <- core_counts_common %>%
  select(depth, upper_age, lower_age, taxa, count) %>%
  spread(key = taxa, value = count) %>%
  arrange(upper_age) #sort by increasing time


### Check temporal resolution of diatom record
median(diff(core_counts_wide_diatoms$upper_age))
plot.ts(diff(core_counts_wide_diatoms$upper_age))
range(core_counts_wide_diatoms$upper_age)

### Estime PrC and plot trends
## extract agedepth variables
agedepth <- core_ra_wide_diatoms[,names(core_ra_wide_diatoms) %in% c("depth", "upper_age", "lower_age")]
diatoms <- core_ra_wide_diatoms[,!(names(core_ra_wide_diatoms) %in% c("upper_age", "lower_age", "depth"))]

# Transform data to Hellinger form
diatoms[is.na(diatoms)] <- 0 #Replace NA (if any) by 0
diat_hell <- decostand(diatoms, method="hellinger")

# Run Principal Curves
diat_prc <- prcurve(diat_hell, method = "ca", trace = TRUE, vary = TRUE, penalty = 1.4)

## Extract position on the curve
scrs_prc <- scores(diat_prc, display = "curve")

# Combine dataframe with ages and depths
diatomsPrC <- cbind(agedepth, scrs_prc)
diatomsPrC$proxy <- "diatoms"

# Plot Pcurves with depth and ages
diat_plt_prc <- ggplot(diatomsPrC, aes(x = upper_age, y = PrC)) +
  geom_line() + geom_point() +
  labs(y = "PrC", x = "Age (cal yr BP)", title = "") +
  ggtitle("Diatoms PrC") +
  theme_bw()
diat_plt_prc


## Clean data and Estimate Principal Curves (PrC) of the F9 pollen core
#Read in F9 pollen record
F9_pollen <- read.csv("data/Fuquene_F9_pollen.csv") %>%
  #rename(depth=ï..depth) %>% #rename extraneous column name  
  select(-c(contains("Lycopo") | contains("Tot.LocTot.100"))) %>%
  #select(-c(66,67,68,70)) %>% #exclude spores 
  rename(depth=ï..depth)

names(F9_pollen)

#this is to transform to tidy format, calculate % and subset more common species
pollen_long <- F9_pollen %>% 
  gather(key = taxa, value = count, -depth, -upper_age, -lower_age) %>%
  #mutate(taxa = plyr::mapvalues(taxa, from = changes[,1], to = changes$new_2)) %>%
  group_by(depth, taxa, upper_age, lower_age) %>%
  summarise(count = sum(count)) %>%
  filter(!count == "0" ) %>% #this is to remove empty samples (rows)
  filter(!upper_age==0.0) %>% #this is to drop extraneous ages
  ungroup() %>%
  group_by(depth) %>%
  mutate(relative_abundance_percent = count / sum(count) * 100) %>%
  ungroup()

range(pollen_long$upper_age)

# filter more abundant taxa; format need to be on long-wide format-->no spreaded 
core_common_taxa <- pollen_long %>%
  group_by(taxa) %>%
  summarise(max_rel_abund = max(relative_abundance_percent)) %>%
  filter(max_rel_abund >= 2) %>%
  arrange(max_rel_abund) %>%
  pull(taxa)

# select from initial table
core_counts_common_F9 <- pollen_long %>%
  filter(taxa %in% core_common_taxa) %>%
  mutate(taxa = factor(taxa, levels = core_common_taxa)) %>%
  arrange(taxa)

#make it wide--relative abundance
core_counts_wide_F9 <- core_counts_common_F9 %>%
  select(depth, upper_age, lower_age, taxa, relative_abundance_percent) %>%
  spread(key = taxa, value = relative_abundance_percent) %>%
  arrange(upper_age)

#make it wide--counts
core_counts_wide_F9 <- core_counts_common_F9 %>%
  select(depth, upper_age, lower_age, taxa, count) %>%
  spread(key = taxa, value = count) %>%
  arrange(upper_age)

# check age core ranges
range(core_counts_wide_F9$upper_age)


### Check temporal resolution of pollen F9 record
mean(diff(core_counts_wide_F9$upper_age))
plot.ts(diff(core_counts_wide_F9$upper_age))

### Run and Plot PrC trends
# remove agedepth columns for Hellinger transformation and bind back
agedepth <- core_counts_wide_F9[,names(core_counts_wide_F9) %in% c("depth", "upper_age", "lower_age")]
pollen <- core_counts_wide_F9[,!(names(core_counts_wide_F9) %in% c("upper_age", "lower_age", "depth"))]
pollen[is.na(pollen)] <- 0 #Replace NA (if any) by 0

# Transform data to Hellinger form
pollen_hell <- decostand(pollen, method="hellinger")

# Run Principal Curves and extract scores for then combining the resulting dataframe with agedepth
pollen.prc <- prcurve(pollen_hell, method = "ca", trace = TRUE, vary = TRUE, penalty = 1.4)
scrs_prc <- scores(pollen.prc, display = "curve")

PollenPrC <- cbind(agedepth, scrs_prc)

# Plot Pcurves with depth and ages
PollenPlot <- ggplot(PollenPrC, aes(x = upper_age, y = PrC)) +
  geom_line() + geom_point() +
  labs(y = "Pollen PrC", x = "Age (cal yr BP)", title = "") +
  ggtitle("Pollen PrC")
PollenPlot

##Clean data and estimate PrC of the Fuquene 7 pollen record
#Read in F7 pollen record
F7_pollen <- read.csv("data/Fuquene_F7_pollen.csv", sep = ",") %>%
  #rename(depth=ï..depth) %>% #rename extraneous column name  
  select(-c(contains("Lycopo"))) %>%
  select(-c(115,116,117,118,119)) %>% #exclude secondary variables (tablets, influx factor, etc)
  filter(row_number() <= n()-1) %>% #Last observation F7 matrix lower age=0. Remove it
  rename(depth=ï..depth)

names(F7_pollen)
str(F7_pollen)

#this is to transform to tidy format, calculate % and subset more common species
pollen_long <- F7_pollen %>% 
  gather(key = taxa, value = count, -depth, -upper_age, -lower_age) %>%
  #mutate(taxa = plyr::mapvalues(taxa, from = changes[,1], to = changes$new_2)) %>%
  group_by(depth, taxa, upper_age, lower_age) %>%
  summarise(count = sum(count)) %>%
  filter(!count == "0" ) %>% #this is to remove empty samples (rows)
  filter(!upper_age==0.0) %>% #this is to drop extraneous ages
  ungroup() %>%
  group_by(depth) %>%
  mutate(relative_abundance_percent = count / sum(count) * 100) %>%
  ungroup()

range(pollen_long$upper_age) 

# filter more abundant taxa; format need to be on long-wide format-->no spreaded 
core_common_taxa <- pollen_long %>%
  group_by(taxa) %>%
  summarise(max_rel_abund = max(relative_abundance_percent)) %>%
  filter(max_rel_abund >= 2) %>%
  arrange(max_rel_abund) %>%
  pull(taxa)

# select from initial table
core_counts_common_F7 <- pollen_long %>%
  filter(taxa %in% core_common_taxa) %>%
  mutate(taxa = factor(taxa, levels = core_common_taxa)) %>%
  arrange(taxa)

#make it wide--relative abundance
core_counts_wide_F7 <- core_counts_common_F7 %>%
  select(depth, upper_age, lower_age, taxa, relative_abundance_percent) %>%
  spread(key = taxa, value = relative_abundance_percent) %>%
  arrange(upper_age)

#make it wide--counts
core_counts_wide_F7 <- core_counts_common_F7 %>%
  select(depth, upper_age, lower_age, taxa, count) %>%
  spread(key = taxa, value = count) %>%
  arrange(upper_age)

# check age core ranges
range(core_counts_wide_F7$upper_age)
median(diff(core_counts_wide_F7$upper_age))
plot.ts(diff(core_counts_wide_F7$upper_age))


## Join F7 and F9 records by common taxa names (columns)
#combine cores (join's analogue package function)
df <- analogue::join(core_counts_wide_F7, core_counts_wide_F9, verbose = TRUE)

#check NA in the list and name the list
listnans <- lapply(df, function(x) sum(is.na(x)))
names(df) <- c("F7", "F9")

# extract merged cores and sort by age
merged <- plyr::ldply(df, data.frame)
# merged <- merged[,-1] #remove .id variable

pollen_F7_F9 <- merged[order(merged$upper_age),] %>% #arrange in increasing depth
  rename(core=.id)

# remove extraneous samples with negative ages
pollen_F7_F9 <- pollen_F7_F9 %>%
  filter(upper_age > -100) %>%
  filter(!lower_age==0.0) #this is to drop extraneous ages

# write results
#write.csv(pollen_F7_F9, "outputs/pollen_F7_F9_counts_2RA.csv")