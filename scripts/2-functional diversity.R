## Step 2: Functional diversity analyses

source("scripts/functions_custom.R")
library(doParallel)
library(cluster)
library(FD)

## Functional diversity of diatom assemblages
diat_comm <- core_counts_wide_diatoms[,-c(1:3)]
diat_traits <- nms_diat[nms_diat$new_2_groups %in% names(diat_comm), ]
length(unique(diat_traits$new_2_groups))

diat_traits <- diat_traits %>%
  distinct(new_2_groups, .keep_all = TRUE) %>%
  select(-c(1,2)) %>%
  column_to_rownames(var = "new_2_groups") %>%
  mutate_if(is.character, as.factor)

diat_comm[is.na(diat_comm)] <- 0
diat_comm <- diat_comm[,order(colnames(diat_comm))]
diat_traits <- diat_traits[order(row.names(diat_traits)),]


# Get standardized effect size of the functional richness  of every assemblage
SesFric <- getSesFric(diat_comm, # Observed diatom assemblages
                      diat_traits, # diatom traits
                      Runs = 99,  # 99x functional richness according to a Null model of community assembly
                      Ncores = 15) # Generate Null communities in parallel on XX cores

#write.table(SesFric, "outputs/SesFric.txt")

# calculate rarefied species richness
rich <- core_counts_wide_diatoms %>%
  replace(is.na(.), 0) %>%
  dplyr::select(-c(3)) %>%
  rename("Depth"=`depth`, "Age_BP"=`upper_age`) %>%
  summarise(Depth, Age_BP, richness=rarefy(round(.[, -(1:2)]), 100)) 

#Sr <- rowSums(ifelse(diat_comm > 0, 1, 0)) # Species richness of the 
# Linear regression model for the relationship between SES Functional richness (response) and species richness (predictor)
FricSrLM <- lm(SesFric ~ rich$richness) 
summary(FricSrLM)

# Display regression result in nice table for the html document
#kable(tidy(FricSrLM)) 

Fric_rich <- read.table("outputs/Fric_rich.txt") %>%
  mutate(MIS = case_when(between(Age_BP, 191000, 243000) ~ "MIS7",
                         between(Age_BP, 70000,135000) ~ "MIS5",
                         between(Age_BP, 26000,60000) ~ "MIS3",
                         between(Age_BP, 0,11800) ~ "MIS1")) %>%
  mutate(MIS = replace_na(MIS, "Glacial"))                         

glacial_samples <- Fric_rich %>%
  filter(MIS=="Glacial")
interg_samples <- Fric_rich %>%
  filter(!MIS=="Glacial")

# plot Frich vs species richness in ggplot
# Fric_rich_taxa <- cbind(Fric_rich, diat_comm$`Aulacoseira ambigua`, diat_comm$`Small Staurosirella`)
# colnames(Fric_rich_taxa)[c(5,6)] <- c("Aulacoseira_ambigua", "Small_Staurosirella")
#write.table(Fric_rich, "outputs/Fric_rich.txt")

library(psych)
library(egg)
corr.test(Fric_rich$SesFric, Fric_rich$richness)$p
corr.test(Fric_rich$SesFric, Fric_rich$richness)$r

Fric_rich_plt <- Fric_rich %>%
  filter(str_detect(MIS, "MIS7|MIS5|MIS3|MIS1")) %>%
  ggplot(., aes(richness, SesFric, shape=MIS, colour=Age_BP)) +
  geom_point(size=5) +
  geom_hline(yintercept = 0, colour="grey", linetype=2) +
  scale_color_distiller(palette = "PuOr", direction=1)+
  geom_point(data=glacial_samples, aes(richness, SesFric), size=2, colour="blue") +
  #geom_smooth(method = "lm") +
  xlab("Rarefied species richness") + ylab("SES Functional richness") +
  theme_article() +
  #guides(colour=guide_legend("cal yrs BP")) +
  theme(axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11),
        axis.title.y = element_text(size=12),
        axis.title.x = element_text(size=12))
Fric_rich_plt

Fric_rich_plt <- Fric_rich %>%
  #filter(str_detect(MIS, "MIS7|MIS5|MIS3|MIS1")) %>%
  ggplot(., aes(richness, SesFric)) +
  geom_point(aes(colour=MIS), size=5) +
  geom_hline(yintercept = 0, colour="grey", linetype=2) +
  scale_colour_viridis_d() +
  #scale_color_discrete(palette = "PuOr", direction=1)+
  #geom_point(data=glacial_samples, aes(richness, SesFric, shape=("Glacial")), size=2) +
  #scale_shape_manual(values = c(1), name="")+
  geom_smooth(method = "lm") +
  xlab("Rarefied species richness") + ylab("SES Functional richness") +
  theme_article() +
  #geom_mark_ellipse(aes(fill = MIS, filter=MIS=="MIS5"))+
  #guides(colour=guide_legend("cal yrs BP")) +
  theme(legend.text = element_text(size=12),
        legend.title = element_blank(),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=13),
        axis.title.x = element_text(size=13))
Fric_rich_plt

Fric_age_plt <- ggplot(Fric_rich, aes(Age_BP, SesFric)) +
  geom_line(size=0.8) +
  #geom_point(aes(colour=richness), size=3) +
  #scale_color_distiller(palette = "YIOrBr", direction = 1)+
  geom_smooth() +
  #scale_x_continuous(breaks=seq(0, 245000, by=2000)) +
  xlab("Calibrated years BP") + ylab("SES Functional richness") +
  theme_article() +
  theme(legend.position="none",
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=13),
        axis.title.x = element_text(size=13)) +
  geom_hline(yintercept = 0, colour="grey", linetype=2) +
  annotate(geom="rect", ymin=-Inf, ymax=Inf, xmin=0, xmax=11800, alpha=0.2, fill="#F8766D") + #MIS1
  annotate(geom="rect", ymin=-Inf, ymax=Inf, xmin=70000, xmax=135000, alpha=0.2, fill="#F8766D") + #MIS5
  annotate(geom="rect", ymin=-Inf, ymax=Inf, xmin=191000, xmax=243000, alpha=0.2, fill="#F8766D")  #MIS7

rich_age_plt <- ggplot(Fric_rich, aes(Age_BP, richness)) +
  geom_line(size=0.8) +
  #geom_point(aes(colour=richness), size=3) +
  #scale_color_distiller(palette = "YIOrBr", direction = 1)+
  geom_smooth() +
  #scale_x_continuous(breaks=seq(0, 245000, by=2000)) +
  xlab("Cal yrs BP") + ylab("Rarefied species richness") +
  theme_article() +
  theme(legend.position="none",
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=11),
        axis.title.y = element_text(size=12),
        axis.title.x = element_blank()) +
  geom_hline(yintercept = 0, colour="grey", linetype=2) +
  annotate(geom="rect", ymin=-Inf, ymax=Inf, xmin=0, xmax=11800, alpha=0.2, fill="#F8766D") + #MIS1
  annotate(geom="rect", ymin=-Inf, ymax=Inf, xmin=70000, xmax=135000, alpha=0.2, fill="#F8766D") + #MIS5
  annotate(geom="rect", ymin=-Inf, ymax=Inf, xmin=191000, xmax=243000, alpha=0.2, fill="#F8766D")  #MIS7


plot_grid(Fric_age_plt, Fric_rich_plt,
          nrow = 1, rel_widths = c(1,1),
          align = "hv", labels = "auto") 

plot_grid(rich_age_plt, Fric_age_plt,
          ncol = 1, rel_widths = c(1,1),
          align = "hv", labels = "auto") 


ggsave <- function(..., bg = 'white') ggplot2::ggsave(..., bg = bg)

ggsave("outputs/for_paper/Frich_rich_age_plots_MIS_def.png",
       plot = last_plot(),
       width = 12,
       height=8,
       units="in",
       dpi = 300)
