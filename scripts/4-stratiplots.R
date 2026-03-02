library(tidyverse)
library(tidypaleo)
library(palinsol)
library(egg)


rm(list = ls())

### Insolation curve
# 1. Define parameters
lat <- 5 * pi/180          # Convert 5N to radians
# July 15 is approximately day 196 of the year. 
# We convert this to 'True Solar Longitude' (lambda) for accuracy in the past.
# At the summer solstice, lambda = pi/2 (90°). mid-July is ~115°.
target_lambda <- day2l(ber78(0), 196) 

# 3. Create a time vector from 240,000 years ago to present (0)
# Time is in years; negative values are in the past.
times <- seq(-240000, 0, by = 1000)

# 4. Calculate daily mean insolation for each time step
# We use the 'la04' (Laskar 2004) solution.
july_insolation <- sapply(times, function(t) {
  orbit <- la04(t) # Get orbital parameters for time t
  Insol(orbit, long = target_lambda, lat = lat)
})

# 5. Plot the results
plot(times / 1000, july_insolation, type = "l", col = "blue",
     xlab = "Thousands of years before present (kyr BP)",
     ylab = expression(paste("July Insolation (W/m"^2, ")")),
     main = expression(paste("July Insolation at 5"^{o}, "N (Last 240ka)")),
     panel.first = grid())

# create a dataframe
insolation <- data.frame(age_bp=abs(times),insolation=july_insolation)

insol_plt <- ggplot(insolation, aes(x = age_bp, y = insolation)) +
  # geom_areah() +
  geom_line() +
  # geom_point() +
  scale_x_reverse(breaks=seq(0, 240000, by=20000)) +
  scale_y_reverse() +
  coord_flip() +
  #theme_article() +
  labs(x = NULL, y = "Annual Mean \n Insolation (W/m2)") +
  theme(axis.text=element_text(size = 12),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y.left = element_blank()) 
insol_plt  

#### Global d18O stack
global_stack_d18O <- read.csv("data/Global_stack_d18O.csv", sep=";")

global_stack_d18O <- global_stack_d18O %>%
  rename(d18O=X18O) %>%
  mutate(age_bp=age_kabp*1000) %>%
  filter(age_bp<250000)

d18O_plt <- ggplot(global_stack_d18O, aes(x = age_bp, y = d18O)) +
  geom_line() +
  scale_y_reverse() +
  scale_x_reverse(breaks=seq(0, 240000, by=20000)) +
  #theme_article() + 
  coord_flip() +
  theme(axis.text=element_text(size = 12),
        axis.text.y = element_blank(),
        axis.ticks.y.left = element_blank(),
        #axis.title.y.left = element_text(size=14),
        axis.title.y.left = element_blank()) +
  labs(x = "Age (cal yr BP)", y = expression(italic(delta)^18*O)) 
d18O_plt  

#### Chinese d18O stack (Cheng et al. 2016)
chinese_stack_d18O <- read.csv("data/chinese_stack_d18O.csv", sep=",")

chinese_stack_d18O <- chinese_stack_d18O %>%
  mutate(age_bp=age_kabp*1000) %>%
  filter(age_bp<250000)

chinese_d18O_plt <- ggplot(chinese_stack_d18O, aes(x = age_bp, y = d18O)) +
  geom_line() +
  scale_y_reverse() +
  scale_x_reverse(breaks=seq(0, 240000, by=20000)) +
  #theme_article() + 
  coord_flip() +
  theme(axis.text=element_text(size = 12),
        axis.text.y = element_blank(),
        axis.ticks.y.left = element_blank(),
        #axis.title.y.left = element_text(size=14),
        axis.title.y.left = element_blank()) +
  labs(x = "Age (cal yr BP)", y = expression(italic(delta)^18*O)) 
chinese_d18O_plt  

Penon_d18O <- read.csv("data/Penon_et_al_2023_d18O.csv", sep=";") %>%
  mutate(Age.ka. = as.numeric(gsub(",", ".", Age.ka., fixed = TRUE)),
         Depth.mm. = as.numeric(gsub(",", ".", Depth.mm., fixed = TRUE)),
         d18O.permil... = as.numeric(gsub(",", ".", d18O.permil..., fixed = TRUE))) %>%
  mutate(age_bp=Age.ka.*1000) 

Penon_d18O_plt <- ggplot(Penon_d18O, aes(x = d18O.permil..., y = age_bp)) +
  geom_lineh() +
  scale_y_reverse(breaks=seq(0, 240000, by=20000)) +
  #geom_lineh(data=chinese_stack_d18O, aes(x = d18O, y = age_bp), colour="blue") +
  #theme_article() +
  theme(strip.text.x=element_text(angle=75, size=12, hjust=0),
        strip.background = element_rect(fill = "white", colour = "white"),
        legend.position = "bottom",
        axis.text.y = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.text.x = element_text(size=12, vjust=0)) +
  labs(x = expression(italic(delta)^18*O), y = NULL) +
  annotate(geom="rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=11800, alpha=0.2, fill="white") + #MIS1
  annotate(geom="rect", xmin=-Inf, xmax=Inf, ymin=70000, ymax=135000, alpha=0.2, fill="white") + #MIS5
  annotate(geom="rect", xmin=-Inf, xmax=Inf, ymin=191000, ymax=243000, alpha=0.2, fill="white") #MIS7
Penon_d18O_plt

## Climate reconstructions (Chevalier et al. 2025)
climate_rec <- read.csv("data/Reconstructions_raw_reconstructions.csv", sep=",") %>%
  rename(age_bp=age..kyr.BP. ) %>%
  mutate(age_bp=age_bp*1000) %>%
  filter(age_bp<240000) %>%
  gather(var, value, -age_bp) %>%
  mutate(var = case_when(
    var == "mat" ~ "Mean Annual Temperature (deg)",
    var == "pann" ~ "Mean Annual Precipitation (mm)",
    var == "ai" ~ "Aridity index",
    TRUE ~ var  # Keep other values unchanged
  ))

clim_rec_plt <- ggplot(climate_rec, aes(x = value, y = age_bp)) +
  geom_lineh() +
  scale_y_reverse(breaks=seq(0, 240000, by=20000)) +
  #scale_y_reverse(breaks=seq(20000, 280000, by=10000)) +
  facet_geochem_gridh(vars(var), rotate_axis_labels=0) +
  #theme_article() +
  theme(strip.text.x=element_text(size=12, vjust=0, hjust=0),
        strip.background = element_rect(fill = "white", colour = "white"),
        legend.position = "bottom",
        axis.text.y = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.text.x = element_text(size=12)) +
  annotate(geom="rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=11800, alpha=0.2, fill="#F8766D") + #MIS1
  annotate(geom="rect", xmin=-Inf, xmax=Inf, ymin=70000, ymax=135000, alpha=0.2, fill="#F8766D") + #MIS5
  annotate(geom="rect", xmin=-Inf, xmax=Inf, ymin=191000, ymax=243000, alpha=0.2, fill="#F8766D") + #MIS7
  labs(x = NULL, y = NULL) 
clim_rec_plt  

# Combine with d18O global stack plot
cowplot::plot_grid(d18O_plt, Penon_d18O_plt, insol_plt, clim_rec_plt,
          nrow = 1, rel_widths = c(1,1,1,3),
          align = "h", axis = "bt", labels = "auto")

# Combine with d18O global stack plot
paleoclim_plt <- cowplot::plot_grid(d18O_plt, chinese_d18O_plt, Penon_d18O_plt, insol_plt,
                   nrow = 1, rel_widths = c(1,1,1,1),
                   align = "h", axis = "bt")

ggsave("outputs/for_paper/paleoclim_comparison.png",
       plot = last_plot(), bg="white",
       width = 18,
       height=8,
       units="in",
       dpi = 300)



myPalette <- c("chartreuse3", "goldenrod", "#483D8B")

#Read in the data
F_diatoms <- read.csv("data/Fuquene_diatoms_counts.csv") %>%
  filter(!core=="M1") %>%
  #filter(core=="F7") %>%
  drop_na(upper_age)

#Matrix of new taxa names and groups
nms_diat <- read.csv("data/nms_changes_diatoms_M1_F7_F9.csv", sep=";")

agedepth <- F_diatoms[, names(F_diatoms) %in% c("depth", "upper_age", "lower_age", "core")] 
diat <- F_diatoms[, !names(F_diatoms) %in% c("depth", "upper_age", "lower_age", "core")]
diat[is.na(diat)] <- 0

names(diat) <- nms_diat[,1]

#Select most common species 
criteria <- 0.25 #% of the total samples

n.occur <- apply(diat>0, 2, sum)
# diat <- diat[, n.occur > 0] #
# rich <- specnumber(diat)
# summary(rich)

diat_red <- diat[, n.occur > (dim(diat)[1])*criteria] #
diatoms <- cbind(agedepth, diat_red)

diat_data <- diatoms %>% 
  gather(key = taxa, value = count, -depth, -upper_age, -lower_age, -core) %>%
  mutate(taxa = plyr::mapvalues(taxa, from = nms_diat[,1], to = nms_diat[,2])) %>%
  mutate(taxa_agg = plyr::mapvalues(taxa, from = nms_diat[,2], to=nms_diat[,4])) %>%
  mutate(group = plyr::mapvalues(taxa_agg, from = nms_diat[,4], to = nms_diat[,5])) %>%
  group_by(depth, upper_age, lower_age) %>%
  mutate(total_sample = sum(count)) %>% 
  filter(!total_sample == "0") %>% #this is to remove empty samples
  mutate(log_total_counts = log10(total_sample+1)) %>%
  mutate(relative_abundance_percent = count / sum(count) * 100) %>%
  mutate(negAge = -upper_age) %>%
  #mutate(AgeCE = upper_age*(-1)+1950) %>%
  #filter(AgeCE >= "0") %>%
  mutate(elapsedTime = round(abs(upper_age - lower_age),0)) %>%
  ungroup() %>%
  filter(!is.na(elapsedTime)) %>%
  mutate(spp = factor(taxa)) %>%
  mutate(spp_agg = factor(taxa_agg)) %>%
  mutate(group = factor(group))

levels(diat_data$spp)
levels(diat_data$spp_agg)
levels(diat_data$group)

diat_data <- diat_data %>%
  mutate(group=fct_recode(group, "Planktonic"="non-motile planktonic",
                          "Tychoplanktonic"="non-motile tychoplanktonic",
                          "Slightly weakly motile "="Slightly_weakly_mod motile "))

theme_set(theme_bw(base_size=14))
diat_plot <- diat_data %>%
  filter(upper_age<240000) %>%
  ggplot(., aes(x = relative_abundance_percent, y = upper_age)) +
  geom_col_segsh() +
  #geom_hline(yintercept = zones, col = "black", lty = 2) +
  #scale_y_reverse() +
  geom_hline(yintercept = 23730, col = "forestgreen", lty = 2, linewidth = 1) +
  geom_hline(yintercept = 33244, col = "blue", lty = 2, linewidth = 1) +
  facet_abundanceh(vars(group), rotate_facet_labels = 0) +
  #scale_x_continuous(breaks=seq(0, 100, by=10)) +
  scale_y_reverse(breaks=seq(0, 250000, by=20000)) +
  scale_colour_manual(values = myPalette) +
  annotate(geom="rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=11800, alpha=0.2, fill="#F8766D") + #MIS1
  annotate(geom="rect", xmin=-Inf, xmax=Inf, ymin=70000, ymax=135000, alpha=0.2, fill="#F8766D") + #MIS5
  annotate(geom="rect", xmin=-Inf, xmax=Inf, ymin=191000, ymax=243000, alpha=0.2, fill="#F8766D") + #MIS7
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_text(size=12)) +
  labs(y = "Age (cal yr BP)") 
diat_plot

ggsave("outputs/for_paper/diat_stratiplot_taxa.png",
       plot = diat_plot,
       width = 10,
       height=8,
       units="in",
       dpi = 300)

#### POLLEN
## Read in change nms table
nms_pollen <- read.csv("data/nms_changes_pollen_F7_F9_2RA.csv", sep = ",") #F7,F9
head(nms_pollen)
#colnames(nms_pollen)[1] <- "taxa"

## Change pollen names and groups
#pollen_df <- read.csv("outputs/pollen_M1_F7_F9_counts_2RA.csv")[-1]
pollen_df <- read.csv("outputs/pollen_F7_F9_counts_2RA.csv")[-1] 

pollen_comp_long <- pollen_df %>% #remove column core
  gather(key = taxa, value = count, -depth, -upper_age, -lower_age, -core) %>%
  mutate(taxa = plyr::mapvalues(taxa, from = nms_pollen[,1], to = nms_pollen$new)) %>%
  mutate(group = plyr::mapvalues(taxa, from = nms_pollen[,2], to = nms_pollen$habitat)) %>%
  mutate(habitat = plyr::mapvalues(taxa, from = nms_pollen[,2], to = nms_pollen$aquatic_habitat)) %>%
  mutate(terrestrial_biome = plyr::mapvalues(taxa, from = nms_pollen[,2], to = nms_pollen$terrestrial_group)) %>%
  mutate(group = replace(group, str_detect(group, "Botryococcus"), "Algae")) %>% #something weird with Botryococcus
  mutate(habitat = replace(habitat, str_detect(habitat, "Botryococcus"), "floating")) %>% #something weird with Botryococcus
  mutate(terrestrial_biome = replace(terrestrial_biome, str_detect(terrestrial_biome, "Botryococcus"), "aquatic")) %>% #something weird with Botryococcus
  group_by(depth, upper_age, lower_age) %>%
  mutate(relative_abundance_percent = count / sum(count) * 100,
         sqrt_pc = sqrt(relative_abundance_percent),
         Norm_sqrt_pc=sqrt_pc/sum(sqrt_pc)*100) %>%
  filter(relative_abundance_percent>5) %>%
  ungroup() %>%
  mutate(group=factor(group)) %>%
  mutate(habitat=factor(habitat)) %>%
  mutate(taxa=factor(taxa)) %>%
  mutate(terrestrial_biome=factor(terrestrial_biome)) %>%
  select(-sqrt_pc) 

levels(pollen_comp_long$group)
levels(pollen_comp_long$habitat)
levels(pollen_comp_long$taxa)
levels(pollen_comp_long$terrestrial_biome)

range(pollen_comp_long$upper_age)
median(diff(pollen_comp_long$upper_age))

## Plot the data
theme_set(theme_bw(base_size=14))
pollen_plot <- pollen_comp_long %>%
  filter(!habitat=="unkown") %>%
  filter(!group=="unknown") %>%
  #filter(relative_abundance_percent > 5) %>%
  #filter(!str_detect(taxa,"Botryococcus|Pediastrum")) %>%
  filter(!terrestrial_biome=="") %>%
  filter(upper_age<240000) %>%
  ggplot(., aes(x = relative_abundance_percent, y = upper_age)) +
  geom_col_segsh() +
  geom_hline(yintercept = 23610, col = "forestgreen", lty = 2, linewidth = 1) + #F7 ends
  geom_hline(yintercept = 26969, col = "darkblue", lty = 2, linewidth = 1) + #F9 begins
  #geom_hline(yintercept = zones, col = "black", lty = 2) +
  #scale_y_reverse() +
  facet_abundanceh(vars(habitat), rotate_facet_labels = 0) +
  labs(x = "Relative abundance (%)", y = "Cal yr BP") +
  #scale_x_continuous(breaks=seq(0, 100, by=20)) +
  scale_y_reverse(breaks=seq(0, 240000, by=20000)) +
  #scale_colour_manual(values = c("royalblue", "olivedrab", "khaki3", "sandybrown")) +
  annotate(geom="rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=11800, alpha=0.2, fill="#F8766D") + #MIS1
  annotate(geom="rect", xmin=-Inf, xmax=Inf, ymin=70000, ymax=135000, alpha=0.2, fill="#F8766D") + #MIS5
  annotate(geom="rect", xmin=-Inf, xmax=Inf, ymin=191000, ymax=243000, alpha=0.2, fill="#F8766D") + #MIS7
  theme(legend.position = "bottom",
        strip.text.x=element_text(angle=0, size=12, vjust=0, hjust=0),
        legend.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size=12)) + labs(y = "Age (cal yr BP)") 
pollen_plot

ggsave("outputs/for_paper/pollen_stratiplot.png",
       plot = pollen_plot,
       width = 10,
       height=8,
       units="in",
       dpi = 300)

## Read in interpolated data
F9_geochem_grainsize <- read.csv("outputs/F9_geochem_grainsize_i.csv", sep = ",")[-1] 

# Prepare data for modeling: grainsize
geochem_grainsize_long <- F9_geochem_grainsize %>%
  #select(-9) %>%
  gather(key = element, value = value, -depth, -upper_age, -lower_age) %>%
  mutate(negAge = -upper_age) %>%
  mutate(elapsedTime = round(abs(upper_age - lower_age),0)) %>%
  filter(elapsedTime < 40000) %>%
  filter(!is.na(as.numeric(value))) %>%
  filter(!is.na(elapsedTime)) %>%
  mutate(element = factor(element)) %>%
  mutate(value_raw=value) %>%
  mutate_at(vars(value), sqrt) %>% #
  rename(value_sqrt=value)

levels(geochem_grainsize_long$element)
range(geochem_grainsize_long$elapsedTime)
range(geochem_grainsize_long$upper_age)

geochem_grainsize_long <- geochem_grainsize_long %>%
  mutate(element=fct_recode(element, "Clay"="clay",
                            "Fine silt"="finesilt",
                            "Coarse silt"="coarsesilt",
                            "Sand" = "sand",
                            "%C"="C_percent",
                            "%N"="N_percent",
                            "d13C"="expression(italic(delta)^13*C",
                            "TC/TN"="TC_TN",
                            "%TIC"="TIC_percent",
                            "%TOC"="TOC_percent",
                            "LOI375"="LOI375")) %>%
  mutate(element=fct_relevel(element,"Clay","Fine silt","Coarse silt","Sand","%C","%N","d13C","TC/TN",
                             "%TIC","%TOC","LOI375"))

levels(geochem_grainsize_long$element)

## Quick plot the trends
theme_set(theme_bw(14))
geochem_plot <- geochem_grainsize_long %>%
  filter(upper_age<240000) %>%
  filter(!element=="LOI375") %>%
  ggplot(., aes(x = value_raw, y = upper_age)) +
  geom_lineh() +
  geom_hline(yintercept = 26969, col = "darkblue", lty = 2, linewidth = 1) + #F9 begins
  #scale_y_reverse() +
  #geom_point() +
  scale_y_reverse(breaks=seq(0, 240000, by=20000)) +
  facet_geochem_gridh(vars(element), rotate_axis_labels =0) +
  labs(x = NULL) +
  labs(y = NULL) +
  annotate(geom="rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=11800, alpha=0.2, fill="#F8766D") + #MIS1
  annotate(geom="rect", xmin=-Inf, xmax=Inf, ymin=70000, ymax=135000, alpha=0.2, fill="#F8766D") + #MIS5
  annotate(geom="rect", xmin=-Inf, xmax=Inf, ymin=191000, ymax=243000, alpha=0.2, fill="#F8766D") + #MIS7
  theme(strip.text.x=element_text(angle=0, size=12, vjust=0, hjust=0),
        strip.background = element_rect(fill = "white", colour = "white"),
        legend.position = "bottom",
        axis.text.y = element_blank(),
        axis.ticks.y.left = element_blank()) +
  ylab("")
geochem_plot

# Prepare data for C/N vs d13C ratio
geochem_grainsize_long <- geochem_grainsize_long %>%
  mutate(MIS = case_when(between(upper_age, 191000, 243000) ~ "MIS7",
                       between(upper_age, 70000,135000) ~ "MIS5",
                       between(upper_age, 26000,60000) ~ "MIS3",
                       between(upper_age, 0,11800) ~ "MIS1")) %>%
  mutate(MIS = replace_na(MIS, "Glacial"))            

TC_TN_d13C <- geochem_grainsize_long %>%
  filter(element %in% c("TC/TN", "d13C")) %>%
  select(upper_age,element,MIS,value_raw) %>%
  spread(key=element, value=value_raw) %>%
  rename(TC_TN="TC/TN") %>%
  ggplot(., aes(TC_TN, d13C)) +
  geom_point(aes(colour=MIS), size=5) +
  scale_colour_viridis_d() +
  #scale_colour_viridis_c() +
  xlab("TC/TN (%)") + ylab(expression(delta^13*C ~ "\U2030")) +
  theme_article() +
  theme(legend.text = element_text(size=13),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=13),
        axis.title.x = element_text(size=13),
        legend.title = element_blank())
TC_TN_d13C

ggsave("outputs/for_paper/CN_d13C.png",
       plot = TC_TN_d13C, bg="white",
       width = 12,
       height=8,
       units="in",
       dpi = 300)

# Composite geochem, grainsize and biota plots
library(cowplot)
# plot_composite <- plot_grid(diat_plot, pollen_plot, grainsize_plot,
#                             nrow = 1, labels = c("Diatoms", "Pollen", "Geochemistry and granulometry"),
#                             rel_widths = c(0.5,0.5,0.5), align = "hv") 

plot_composite <- plot_grid(diat_plot, pollen_plot, geochem_plot, clim_rec_plt,
                            nrow = 1, 
                            labels = c("Diatoms", "Pollen", "Geochemistry & granulometry", "Climate reconstructions"), label_size = 12,
                            rel_widths = c(0.5,0.5,1,0.3,0,2), align = "hv") 

# Combine with d18O global stack plot
library(patchwork)

# plot_composite <- diat_plot + pollen_plot + grainsize_plot + clim_rec_plt + paleoclim_plt +
#   patchwork::plot_layout(widths = c(1,1,1,1,1)) + plot_annotation(tag_levels = "A") 
# plot_composite

plot_composite <- cowplot::plot_grid(diat_plot, pollen_plot, geochem_plot, clim_rec_plt, chinese_d18O_plt, insol_plt,
                   nrow = 1, rel_widths = c(1,1,1.5,1,0.3,0.3),
                   align = "h", axis = "bt", labels = "auto")
plot_composite

ggsave("outputs/for_paper/diat_pollen_geochem_climrec_6.svg",
       plot = plot_composite, bg="white",
       width =22,
       height=8,
       units="in",
       dpi = 300)

