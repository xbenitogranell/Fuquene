######################
#### STRATIPLOTS #####
######################

## To work with interpolated dataset (yearly)

source("scripts/functions_custom.R") #functions to lag datasets and model (Blas) and binning (Seddon)
library(tidyverse)
library(tidypaleo)

## Elegant way to read multiple excel sheets per Excel file
# load names of excel files 
files <- list.files(path = "data/", full.names = TRUE, pattern = ".xlsx")

# create function (see t that transpose the dataframe for later gather)
read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  sapply(sheets, function(f) as.data.frame(readxl::read_excel(filename, sheet = f, col_names = TRUE)), 
         simplify = FALSE)
}

# execute function for all excel files in "files"
all_data <- lapply(files, read_excel_allsheets)

#List names
names(all_data[[1]])

# extract datasets
aquatics <- all_data[[1]]$IntCln_aquatics 
aquatics <- aquatics[, !(names(aquatics) %in% c("age", "depth"))]
totals <- apply(aquatics, 1, sum) #suma les abundancies
aquatics <- aquatics / totals * 100 #calcul proporcions

pollen <- all_data[[1]]$IntClen_pollen[-1]
age <- pollen[,"age"] #all datasets will share this variable

pollen <- pollen[,!(names(pollen) %in% c("age", "depth"))]
pollen <- pollen[, !(names(pollen) %in% names(aquatics))]

totals <- apply(pollen, 1, sum) #suma les abund?ncies
pollen <- pollen / totals * 100 #calcul proporcions


diatoms <- all_data[[1]]$IntCln_diatgrp %>%
  rename(depth=Depth) %>%
  select(-c("depth", "Age"))
totals <- apply(diatoms, 1, sum) #suma les abund?ncies
diatoms <- diatoms / totals * 100 #calcul proporcions

geochem <- all_data[[1]]$IntClen_geochem %>%
  rename(depth=Depth) %>%
  select(-c("depth", "Age"))

range(age)

#make the binning - 100-yr
Pollen_bin <- binFunc(as.data.frame(pollen), as.numeric(age), 100, min(age), max(age)) 
diatoms_bin <- binFunc(as.data.frame(diatoms), as.numeric(age), 100, min(age), max(age)) 
aquatics_bin <- binFunc(as.data.frame(aquatics), as.numeric(age), 100, min(age), max(age)) 
geochem_bin <- binFunc(as.data.frame(geochem), as.numeric(age), 100, min(age), max(age)) 

# cbind dataframes
combinedData <- cbind(Pollen_bin, diatoms_bin, aquatics_bin, geochem_bin)

# make name vectors to later replace taxon names to assemblage group
nms_list <- data.frame(taxa = c(colnames(diatoms), colnames(pollen), colnames(aquatics), colnames(geochem)), 
                       group = c(rep("Diatoms",length(diatoms)), rep("Pollen",length(pollen)), 
                                 rep("Aquatics",length(aquatics)),
                                 rep("Geochemistry", length(geochem))))

nms_aquatics <- data.frame(taxa = c(colnames(aquatics), colnames(diatoms)),
                           group = c("floating", "emergent", "submerged", "emergent", "floating", "submerged",
                                     "submerged", "unknown", "emergent", "unknown", "floating", "floating", 
                                     rep("Diatoms", length(diatoms))))

nms_aquatics <- data.frame(taxa = c(colnames(aquatics)),
                           group = c("floating", "emergent", "submerged", "emergent", "floating", "submerged",
                                     "submerged", "unknown", "emergent", "unknown", "floating", "floating"))


#gather dataset for plotting
df_long <- combinedData %>% 
  as.data.frame() %>%
  mutate(age=row.names(combinedData)) %>%
  mutate(age=as.numeric(age)) %>%
  gather(key=taxa, value=relative_abundance_percent, -age) %>%
  mutate(assemblage = plyr::mapvalues(taxa, from = nms_list$taxa, to = nms_list$group)) %>%
  #mutate(assemblage2 = plyr::mapvalues(taxa, from = nms_aquatics$taxa, to = nms_aquatics$group)) %>%
  mutate(taxa=factor(taxa)) %>%
  mutate(taxa=reorder(taxa,relative_abundance_percent)) %>%
  group_by(assemblage,taxa) %>%
  arrange(relative_abundance_percent) %>%
  ungroup() %>%
  drop_na(relative_abundance_percent)

levels(df_long$taxa)

## Using tidypaleo R package (https://fishandwhistle.net/post/2018/stratigraphic-diagrams-with-tidypaleo-ggplot2/)
theme_set(theme_bw(12))
#theme_set(theme_paleo())

## subset biotic assemblages
assemblage_data <- df_long %>% 
  filter(assemblage!="Geochemistry") %>%
  filter(assemblage!="Pollen") %>%
  #filter(relative_abundance_percent>4) %>%
  mutate(assemblage=factor(assemblage)) %>%
  group_by(age, taxa) %>%
  mutate(taxa = fct_relevel(taxa, "Non.motile.planktonic", "Non.motile.tychoplanktonic", "Non.motile.benthic", "Slightly...weakly..mod.motile",
                            "Highly.motile", "Cyperaceae", "Debarya", "Hydrocotyle", "Isoetes", "Ludwigia", "Myriophyllum", "Polygonum",
                            "Potamogeton", "Spirogyra", "Umbelliferae", "Typha", "Zygnema")) %>%
  mutate(assemblage = plyr::mapvalues(taxa, from = nms_aquatics$taxa, to = nms_aquatics$group)) %>%
  
  # mutate(floating=sum(relative_abundance_percent[taxa=="Ludwigia" | taxa=="Hydrocotyle"| taxa=="Spirogyra"| taxa=="Zygnema"])) %>%
  # mutate(emergent=sum(relative_abundance_percent[taxa=="Polygonum" | taxa=="Typha" | taxa=="Cyperaceae"])) %>%
  # mutate(submerged=sum(relative_abundance_percent[taxa=="Potamogeton" | taxa=="Isoetes" | taxa=="Myriophyllum"])) %>%
  droplevels()
  # mutate(taxa = fct_relevel(taxa, "Alnus", "Aragoa", "Asteraceae.Asteroidea", "Asteraceae.Liguliflora", "Borreria",
  #                           "Caryophyllaceae", "Cecropia", "Dodonea", "Ericaceae", "Geranium", "Hedyosmum", "Hypericum",
  #                           "Ilex", "Lycopodium", "Lycopodium.Foveolate", "Miconia", "Myrica", "Poaceae", "Podocarpus", "Polylepis.Acaena",
  #                           "Puya", "Quercus", "Rumex", "Sordariaceae", "Sporormiella", "Tot.LocTot.100", "Weinmannia",
  #                           "Cyperaceae", "Debarya", "Hydrocotyle", "Isoetes", "Ludwigia", "Myriophyllum", "Polygonum",
  #                           "Potamogeton", "Spirogyra", "Umbelliferae", "Typha", "Zygnema",
  #                           "Non.motile.planktonic", "Non.motile.tychoplanktonic", "Non.motile.benthic", "Slightly...weakly..mod.motile",
  #                           "Highly.motile")) %>%
  #mutate(assemblage = fct_relevel(assemblage, "Pollen", "Aquatics", "Diatoms"))

## Plot
assemblage_plot <- ggplot(assemblage_data, aes(x = relative_abundance_percent, y = age)) +
  geom_col_segsh(size=1.3) +
  scale_y_reverse() +
  geom_areah() +
  facet_abundanceh(vars(assemblage), rotate_facet_labels = 70) +
  #facet_abundanceh(vars(taxa), grouping = vars(assemblage), rotate_facet_labels = 70) +
  #geom_lineh_exaggerate(exaggerate_x = 5, col = "grey70", lty = 2) +
  labs(x = "Relative abundance (%)", y = "cal years BP", colour="Assemblage") +
  #ggtitle("Assemblages") +
  theme(axis.text.x = element_text(size = 9),
        legend.position = "bottom",
        legend.title = element_blank()) 
assemblage_plot

# Subset geochemistry
geochem_data <- df_long %>% 
  filter(assemblage=="Geochemistry") %>%
  rename(value=relative_abundance_percent) %>%
  rename(variable=taxa)

# Plot geochemistry using tidypaleo R package
theme_set(theme_bw(12))

geochem_plot <- ggplot(geochem_data, aes(x = value, y = age)) +
  geom_lineh() +
  #geom_smooth()+
  scale_y_reverse() +
  geom_point() +
  facet_geochem_gridh(vars(variable)) +
  labs(x = NULL) +
  labs(y = NULL)
#ggtitle("XRF")
geochem_plot

## combine plots
library(patchwork)
strat_plt2 <- wrap_plots(
  assemblage_plot + 
    theme(strip.text.y = element_text(size=10),
          strip.text.x = element_text(size=10)),
  #panel.grid = element_blank()),
  geochem_plot +
    theme(axis.text.x = element_text(size = 10),
          strip.text.x = element_text(size = 10, angle = 90), 
          axis.ticks.y.left = element_blank(),
          axis.text.y=element_blank(),  #remove y axis labels
          axis.ticks.y=element_blank()),
  nrow = 1,
  widths = c(6, 2)
)
strat_plt2

ggsave("outputs/aquatic_proxies_stratplot_Fuquene.png", strat_plt2, height = 11, width = 15, dpi=300)



