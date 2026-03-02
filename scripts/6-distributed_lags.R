## Step 4- Distributed Hierarchical Generalized Additive Models (DHGAM)

library(mvgam)
library(tidyverse)
library(egg)

## Prepare the dataset for modeling
df <- core_counts_wide_diatoms %>%
  select(depth,`Aulacoseira ambigua`,`Discostella stelligera`,`Small Fragilarioids`,`Cymbella minor`,`Gomphonema angustum`) %>%
  bind_cols(prc_F7F9_df_i) %>%
  mutate(time=1:NROW(prc_F7F9_df_i)) %>%
  rename(aulac='Aulacoseira ambigua',
         cymbm='Cymbella minor',
         disc='Discostella stelligera',
         gomphn='Gomphonema angustum',
         fragil='Small Fragilarioids')

#write.csv(df,'outputs/df_for_distributedlags.csv')

# Read in the previously saved dataset
df <- read.csv("outputs/df_for_distributedlags.csv")[-1]
prc_F7F9_df_i <- read.csv("outputs/F7_F9/prc_geochem_F7F9_df_i_April2025.csv")[-1]

#plot.ts(diff(prc_F7F9_df_i$upper_age)) #check ellapsed time between samples

df_long <- df %>%
  gather(series, counts, -upper_age, -lower_age,-time, -depth, -c(names(prc_F7F9_df_i))) %>%
  select(-c("diatoms")) #diatoms are the PrC
  
df_long_ts <- df_long %>%
  dplyr::mutate(series = as.factor(series)) %>%
  mutate(MIS = case_when(between(upper_age, 191000, 243000) ~ "MIS7",
                         between(upper_age, 70000,135000) ~ "MIS5",
                         between(upper_age, 26000,60000) ~ "MIS3",
                         between(upper_age, 0,11800) ~ "MIS1")) %>%
  mutate(MIS = replace_na(MIS, "Glacial"))     


dplyr::glimpse(df_long_ts)
levels(df_long_ts$series)
max(df_long_ts$time)

plot_mvgam_series(data = df_long_ts, y = 'counts', series = 'all')
plot_mvgam_series(data = df_long_ts, y = 'counts', series = 5)


unique(df_long_ts$series)
spp <- "gomphn"
  
df_long_ts %>% 
  dplyr::filter(str_detect(series, spp)) %>%
  #ggplot(aes(x = terrestrial, y = log(counts + 1), colour=upper_age/1000)) +
  ggplot(aes(x = terrestrial, y = log(counts + 1), colour=MIS)) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 12),
              col = 'darkred', fill = "#A25050") +
  labs(title = paste0(spp),
       y = "log(counts)", 
       x = 'PrC terrestrial pollen') +
  theme_article() + guides(colour=guide_legend(title="Age (kyrs BP)")) +

  df_long_ts %>% 
  dplyr::filter(str_detect(series, spp)) %>%
  #ggplot(aes(x = aquatic, y = log(counts + 1), colour=upper_age/1000)) +
  ggplot(aes(x = aquatic, y = log(counts + 1), colour=MIS)) +
  
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 12),
              col = 'darkred', fill = "#A25050") +
  labs(y = NULL, 
       x = 'PrC aquatic pollen') +
  theme_article() + theme(legend.position = "none") +

  df_long_ts %>% 
  #dplyr::filter(str_detect(series, spp)) %>%
  #ggplot(aes(x = map_i, y = log(counts + 1),colour=upper_age/1000)) +
  ggplot(aes(x = map_i, y = log(counts + 1), colour=MIS)) +
  
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 12),
              col = 'darkred', fill = "#A25050") +
  labs(y = NULL, 
       x = 'map_i') +
  theme_article() +theme(legend.position = "none") +
  
  df_long_ts %>% 
  dplyr::filter(str_detect(series, spp)) %>%
  #ggplot(aes(x = mat_i, y = log(counts + 1),colour=upper_age/1000)) +
  ggplot(aes(x = mat_i, y = log(counts + 1), colour=MIS)) +
  
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 12),
              col = 'darkred', fill = "#A25050") +
  labs(y = NULL, 
       x = 'mat_i') +
  theme_article() +theme(legend.position = "none") +
  
  df_long_ts %>% 
  dplyr::filter(str_detect(series, spp)) %>%
  #ggplot(aes(x = ai_i, y = log(counts + 1), colour=upper_age/1000)) +
  ggplot(aes(x = ai_i, y = log(counts + 1), colour=MIS)) +
  
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 12),
              col = 'darkred', fill = "#A25050") +
  labs(y = NULL, 
       x = 'ai_i') +
  theme_article() +theme(legend.position = "none") +

df_long_ts %>% 
  dplyr::filter(str_detect(series, spp)) %>%
  #ggplot(aes(x = d13C, y = log(counts + 1),colour=upper_age/1000)) +
  ggplot(aes(x = d13C, y = log(counts + 1), colour=MIS)) +
  
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 12),
              col = 'darkred', fill = "#A25050") +
  labs(y = NULL, 
       x = 'd13C') +
  theme_article() +theme(legend.position = "none") +
  
df_long_ts %>% 
  dplyr::filter(str_detect(series, spp)) %>%
  #ggplot(aes(x = C_percent, y = log(counts + 1),colour=upper_age/1000)) +
  ggplot(aes(x = C_percent, y = log(counts + 1), colour=MIS)) +
  
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 12),
              col = 'darkred', fill = "#A25050") +
  labs(y = NULL, 
       x = 'C_percent') +
  theme_article() + theme(legend.position = "none") +
  
df_long_ts %>% 
  dplyr::filter(str_detect(series, spp)) %>%
  #ggplot(aes(x = TC_TN, y = log(counts + 1),colour=upper_age/1000)) +
  ggplot(aes(x = TC_TN, y = log(counts + 1), colour=MIS)) +
  
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 12),
              col = 'darkred', fill = "#A25050") +
  labs(y = NULL, 
       x = 'TC_TN') +
  theme_article() + theme(legend.position = "none") +

df_long_ts %>% 
  dplyr::filter(str_detect(series, spp)) %>%
  #ggplot(aes(x = sand, y = log(counts + 1),colour=upper_age/1000)) +
  ggplot(aes(x = sand, y = log(counts + 1), colour=MIS)) +
  
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 12),
              col = 'darkred', fill = "#A25050") +
  labs(y = NULL, 
       x = 'sand') +
  theme_article() + theme(legend.position = "none") +


df_long_ts %>% 
  dplyr::filter(str_detect(series, spp)) %>%
  #ggplot(aes(x = coarsesilt, y = log(counts + 1),colour=upper_age/1000)) +
  ggplot(aes(x = coarsesilt, y = log(counts + 1), colour=MIS)) +
  
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 12),
              col = 'darkred', fill = "#A25050") +
  labs(y = NULL, 
       x = 'coarsesilt') +
  theme_article() + theme(legend.position = "none") +


df_long_ts %>% 
  dplyr::filter(str_detect(series, spp)) %>%
  #ggplot(aes(x = finesilt, y = log(counts + 1),colour=upper_age/1000)) +
  ggplot(aes(x = finesilt, y = log(counts + 1), colour=MIS)) +
  
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 12),
              col = 'darkred', fill = "#A25050") +
  labs(y = NULL, 
       x = 'finesilt') +
  theme_article() + theme(legend.position = "none") +

df_long_ts %>% 
  dplyr::filter(str_detect(series, spp)) %>%
  #ggplot(aes(x = clay, y = log(counts + 1),colour=upper_age/1000)) +
  ggplot(aes(x = clay, y = log(counts + 1), colour=MIS)) +
  
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 12),
              col = 'darkred', fill = "#A25050") +
  labs(y = NULL, 
       x = 'clay') +
  theme_article()  + theme(legend.position = "none")

ggsave(paste0("outputs/for paper/",spp,"_drivers.png"),
       plot = last_plot(),
       width = 12,
       height=8,
       units="in",
       dpi = 300)


## function to lag
lagard <- function(x, n_lag = 6) {
  n <- length(x)
  X <- matrix(NA, n, n_lag)
  for (i in 1:n_lag){
    X[i:n, i] <- x[i:n - i + 1]
  } 
  X
}

# here select which variables to lag
names(df_long_ts)
driver <- "terrestrial"
#covariate <- "aquatic"

lag_driver <- do.call(rbind, lapply(seq_along(levels(df_long_ts$series)), function(x){
  df_long_ts %>%
    dplyr::filter(series == levels(df_long_ts$series)[x]) %>%
    dplyr::arrange(time) %>%
    #dplyr::select(terrestrial, time) %>%
    dplyr::select(driver, time) %>%
    #dplyr::pull(terrestrial) -> tempdat
    dplyr::pull(driver) -> tempdat
  
  lag_driver <- lagard(tempdat, 6)
  tail(lag_driver, NROW(lag_driver) - 5)
}))

#dim(terrestrial)[1]
dim(lag_driver)[1]

#head(terrestrial, 10)
head(lag_driver, 10)

df_long_ts %>%
  dplyr::arrange(series, time) %>%
  dplyr::filter(time > 5) -> df_long_ts

data_all <- list(
  lag = matrix(0:5, nrow(df_long_ts), 6, byrow = TRUE),
  counts = df_long_ts$counts,
  #aquatic = df_long_ts$aquatic,
  #covariate = df_long_ts[,names(df_long_ts) %in% covariate],
  time = df_long_ts$time,
  series = df_long_ts$series,
  elapsedTime = df_long_ts$lower_age - df_long_ts$upper_age,
  MIS = factor(df_long_ts$MIS), #here make it numeric instead of factor?
  driver= df_long_ts[,names(df_long_ts) %in% driver])

# check temp resolution of lags
lags<-1:6
median(data_all$elapsedTime) * lags

#data_all$terrestrial <- terrestrial
data_all$lag_driver <- lag_driver

#The dimensions of all the objects need to match up
dim(data_all$lag)
## [1] 375   6
#dim(data_all$terrestrial)
dim(data_all$lag_driver)

## [1] 375   6
length(data_all$time)
## [1] 375

# aulac <- which(data_all$series == 'aulac')
# data_aulac <- lapply(data_all, function(x){
#   if(is.matrix(x)){
#     x[aulac, ]
#   } else {
#     x[aulac]
#   }
# })
# 
# data_aulac$lag_driver <- lag_driver
# 
# # 
# mod1 <- gam(counts ~
#               te(terrestrial, lag, k = 6) +
#               s(aquatic),
#             family = poisson(),
#             data = data_aulac,
#             method = 'REML')
# 
# summary(mod1)
# draw(mod1, select = 1)
# appraise(mod1)

unique(data_all$series)

weights_aulac <- weights_cymbm <- weights_disc <- weights_gomphn <- weights_fragil <- matrix(1, ncol = ncol(data_all$lag), nrow = nrow(data_all$lag))

weights_aulac[!(data_all$series == 'aulac'), ] <- 0
weights_cymbm[!(data_all$series == 'cymbm'), ] <- 0
weights_disc[!(data_all$series == 'disc'), ] <- 0
weights_gomphn[!(data_all$series == 'fragil'), ] <- 0
weights_fragil[!(data_all$series == 'gomphn'), ] <- 0

data_all$weights_aulac <- weights_aulac
data_all$weights_cymbm <- weights_cymbm
data_all$weights_disc <- weights_disc
data_all$weights_gomphn <- weights_gomphn
data_all$weights_fragil <- weights_fragil

# lagged model
mod_lag <- gam(counts ~ 
              # Hierarchical intercepts
              s(series, bs = 're') +
              # Smooths of time to try and capture autocorrelation
              s(time, by = series, bs="fs", k = 20) +
              
              # Smooths of aquatic pollen effects with random intercepts
              #s(covariate, by = series, k = 5) +
              
              # Distributed lags of terrestrial pollen
              te(lag_driver, lag, k = 4, by = weights_aulac) +
              te(lag_driver, lag, k = 4, by = weights_cymbm) +
              te(lag_driver, lag, k = 4, by = weights_disc) +
              te(lag_driver, lag, k = 4, by = weights_gomphn)+
              te(lag_driver, lag, k = 4, by = weights_fragil),
              
            weights = elapsedTime/mean(elapsedTime), #this is the "years mud slide"

            family = nb(),
            data = data_all,
            control = list(nthreads = 6),
            method = 'REML')


appraise(mod_lag)
summary(mod_lag)
k.check(mod_lag)
draw(mod_lag, select = 1:3)
# draw(mod_lag)


#no-lagged model
mod_nolag <- gam(counts ~ 
                 # Hierarchical intercepts
                 s(series, bs = 're') +
                 # Smooths of time to try and capture autocorrelation
                 s(time, by = series, bs="fs", k = 20) +
                   
               # Smooths of covariate effects with random intercepts
               #s(covariate, by = series, k = 5) +
                 
               # Smooths of driver effects with random intercepts
               s(driver, by = series, k = 5),

               weights = elapsedTime/mean(elapsedTime), #this is the "years mud slide"
               
               family = nb(),
               data = data_all,
               control = list(nthreads = 6),
               method = 'REML')

#appraise(mod_lag)
summary(mod_nolag)
appraise(mod_nolag)
draw(mod_nolag)
k.check(mod_nolag)


### it can be argued that lag effects can be stronger at certain MIS. driver-lag effect should then change with MIS
# knots <- list(time = c(6,24,36,37,42,56,73,80)) #end points of MIS interglacials
# knots <- list(time = c(6,42,56,73,80)) #end points of MIS interglacials

system.time(
  mod_temp_lag <- gam(counts ~ 
                 # Hierarchical intercepts
                 s(series, bs = 're') +
                 # Smooths of time to try and capture autocorrelation
                 s(time, by = series, bs="fs", k=20) + 
                   
                 # Smooths of aquatic pollen effects with random intercepts
                 #s(covariate, by = series, k = 5) +
                 
                 # Distributed lags of terrestrial pollen
                 te(lag_driver, lag, k = 4, by = weights_aulac) +
                 te(lag_driver, lag, k = 4, by = weights_cymbm) +
                 te(lag_driver, lag, k = 4, by = weights_disc) +
                 te(lag_driver, lag, k = 4, by = weights_gomphn)+
                 te(lag_driver, lag, k = 4, by = weights_fragil) +
                 
                 ti(lag_driver, lag, time) + MIS + s(time,by=MIS),
               
               weights = elapsedTime/mean(elapsedTime), #this is the "years mud slide"
               family = nb(),
               data = data_all,
               control = list(nthreads = 6),
               method = 'REML'))

summary(mod_temp_lag)
mod_temp_lag$aic
k.check(mod_temp_lag)
draw(mod_temp_lag, select = 6:10)

#Compare different model fits using AIC
AIC_table <- AIC(mod_lag, mod_nolag, mod_temp_lag) %>%
  rownames_to_column(var= paste0(driver,"_Model"))%>%
  mutate(data_source = rep(c("diatom_data")))%>%
  group_by(data_source)%>%
  mutate(deltaAIC = AIC - min(AIC))%>%
  ungroup()%>%
  dplyr::select(-data_source)%>%
  mutate_at(.vars = vars(df,AIC, deltaAIC), 
            .funs = funs(round,.args = list(digits=0)))

AIC_table

# Extract Adjusted R-Squared Values 
adj_r_squared_1 <- summary(mod_lag)$r.sq
adj_r_squared_2 <- summary(mod_nolag)$r.sq
adj_r_squared_3 <- summary(mod_temp_lag)$r.sq

results_data_temp <- data.frame(
  Model = c("mod_lag", "mod_nolag", "mod_temp_lag"),
  Adjusted_R_Squared = c(adj_r_squared_1, adj_r_squared_2, adj_r_squared_3),
  AIC = c(AIC_table$AIC[1],AIC_table$AIC[2],AIC_table$AIC[3]),
  driver = driver
)
results_data_temp

results_data <- results_data_temp

terrestrial_aquatic <- rbind(results_data, results_data_temp)
terrestrial_aquatic_map <- rbind(terrestrial_aquatic, results_data)
terrestrial_aquatic_map_ai <- rbind(terrestrial_aquatic_map, results_data_temp)
terrestrial_aquatic_map_ai_mat <- rbind(terrestrial_aquatic_map_ai, results_data_temp)

#write.csv(terrestrial_aquatic_map_ai_mat, "outputs/for paper/models_results_AIC.csv")

spp_labels <- c(
  'aulac'="Aulacoseira ambigua",
  'cymbm'="Cymbella minor",
  'disc'="Discostella stelligera",
  'fragil'="Small Fragilarioids",
  "gomphn" = "Gomphonema spp"
)

spp_labeller <- as_labeller(spp_labels)

lag_labels <- c(round(median(data_all$elapsedTime) * lags),0)

plot_dist_lags = function(model, data_all){
  
  require(viridis)
  require(egg)
  
  all_species <- levels(data_all$series)
  
  # Loop across species to create the effect plot dataframe
  sp_plot_dat <- do.call(rbind, lapply(all_species, function(sp){
    
    # Zero out all predictors to start the newdata
    newdata <- lapply(data_all, function(x){
      if(is.matrix(x)){
        matrix(0, nrow = nrow(x), ncol = ncol(x))
      } else {
        rep(0, length(x))
      }
    })
    
    # Modify to only focus on the species of interest
    newdata$series <- rep(sp, nrow(data_all$lag))
    newdata$lag <- data_all$lag
    which_weightmat <- grep(paste0('weights_', 
                                   tolower(sp)), 
                            names(newdata))
    newdata[[which_weightmat]] <- matrix(1, nrow = nrow(newdata[[which_weightmat]]),
                                         ncol = ncol(newdata[[which_weightmat]]))
    
    # Calculate predictions for when driver is zero to find the baseline
    # value for centring the plot
    if(inherits(model, 'mvgam')){
      preds <- predict(model, newdata = newdata, 
                       type = 'link', process_error = FALSE)
      preds <- apply(preds, 2, median)
    } else {
      preds <- predict(model, newdata = newdata, type = 'link')
    }
    
    offset <- mean(preds)
    plot_dat <- do.call(rbind, lapply(seq(1:6), function(lag){
      # Set up prediction matrix for mintemp; 
      # use a sequence of values across the full range of observed values
      newdata$lag_driver <- matrix(0, ncol = ncol(newdata$lag),
                                    nrow = nrow(newdata$lag))
      

      newdata$lag_driver[,lag] <- seq(min(data_all$lag_driver),
                                       max(data_all$lag_driver),
                                       length.out = length(newdata$time))
      
      
      # Predict on the link scale and shift by the offset 
      # so that values are roughly centred at zero
      if(inherits(model, 'mvgam')){
        preds <- predict(model, newdata = newdata, 
                         type = 'link', process_error = FALSE) 
        preds <- apply(preds, 2, median)
      } else {
        preds <- predict(model, newdata = newdata,
                         type = 'link') 
      }
      preds <- preds - offset
      
      # data.frame(lag = lag,
      #            preds = preds,
      #            terrestrial = seq(min(data_all$terrestrial),
      #                          max(data_all$terrestrial),
      #                          length.out = length(newdata$time)))
      
      data.frame(lag = lag,
                 preds = preds,
                 lag_driver = seq(min(data_all$lag_driver),
                                   max(data_all$lag_driver),
                                   length.out = length(newdata$time)))
      
      
    }))
    plot_dat$species <- sp
    
    plot_dat
  }))
  
  # Build the facetted distributed lag plot
  ggplot(data = sp_plot_dat %>%
           dplyr::mutate(lag = as.factor(lag)),
         aes(x = lag_driver, y = preds, 
             colour = lag, fill = lag)) +
    facet_wrap(~ species, scales = 'free', labeller = spp_labeller, nrow = 2) +
    geom_hline(yintercept = 0) +
    # Use geom_smooth, though beware these uncertainty
    # intervals aren't necessarily correct
    geom_smooth() +
    #scale_fill_manual(labels = c("A", "B","C","D","E","G")) +
    scale_fill_viridis(discrete = TRUE, labels= lag_labels, name = "Lag (cal yr)") +
    scale_colour_viridis(discrete = TRUE, labels= lag_labels, name = "Lag (cal yr)") +
    labs(x = paste0("PrC ",driver," pollen (z-scored)"),
         y = 'Partial effect') +
    theme_article() +
    theme(axis.text = element_text(size=14),
          axis.title = element_text(size=14),
          legend.text = element_text(size=13),
          legend.title = element_text(size=14),
          strip.text = element_text(size=14, face = "italic"))
  }

plot_dist_lags(mod_lag, data_all)

ggsave <- function(..., bg = 'white') ggplot2::ggsave(..., bg = bg)

ggsave(paste0("outputs/for_paper/HGAM_distributed_lags_",driver,".png"),
       plot = last_plot(),
       width = 12,
       height=8,
       units="in",
       dpi = 300)


