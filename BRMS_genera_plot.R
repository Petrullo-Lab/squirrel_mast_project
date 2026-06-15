asv_list <- readRDS("jaccard_aligned_by_genus_filtered_new.RDS")
library(tidybayes)

library(marginaleffects)
library(dplyr)
library(ggplot2)
library(tidyverse)

mmmm = readRDS("/Users/gabri/Desktop/mouse/brms_model_Faecalibacterium.rds")
#View(mmmm)




# Option 1: use your actual dataframe name (check with ls())
# If it's called e.g. `df` or `squirrel_data`, use that

library(tidyverse)
library(marginaleffects)

model_dir <- "/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/Squirrels/masting_project/BRMS_data_again/model_results_genera/"

rds_files <- list.files(model_dir, pattern = "\\.rds$", full.names = TRUE)
rds_files[1:10]
vars <- c("mother_offspring", "masting", "same_season", 
          "year_dist", "same_year", "dist_m", "same_sex", 
          "run", "age_dist")

run_model <- function(filepath) {
  
  genus <- gsub("brms_model_|\\.rds$", "", basename(filepath))
  message("Processing: ", genus)
  
  mmmm <- readRDS(filepath)
  
  newdata_grid <- datagrid(
    model            = mmmm,
    mother_offspring = c(0, 1),
    masting          = c(0, 1),
    same_season      = mean(mmmm$data$same_season),
    year_dist        = mean(mmmm$data$year_dist),
    same_year        = mean(mmmm$data$same_year),
    dist_m           = mean(mmmm$data$dist_m),
    same_sex         = mean(mmmm$data$same_sex),
    run              = mean(mmmm$data$run),
    age_dist         = mean(mmmm$data$age_dist)
  )
  
  main1 <- lapply(vars, function(v) {
    avg_slopes(mmmm,
               variables  = v,
               newdata    = newdata_grid,
               re_formula = NA,
               ndraws     = 500)
  }) %>% bind_rows() %>%
    select(term, estimate, conf.low, conf.high)
  
  int1 <- avg_comparisons(mmmm,
                          variables  = list(mother_offspring = c(0,1), masting = c(0,1)),
                          cross      = TRUE,
                          re_formula = NA,
                          ndraws     = 500) %>%
    mutate(term = "mother_offspring × masting") %>%
    select(term, estimate, conf.low, conf.high)
  
  int2 <- avg_comparisons(mmmm,
                          variables  = list(mother_offspring = c(0,1), same_season = c(0,1)),
                          cross      = TRUE,
                          re_formula = NA,
                          ndraws     = 500) %>%
    mutate(term = "mother_offspring × same_season") %>%
    select(term, estimate, conf.low, conf.high)
  
  bind_rows(main1, int1, int2) %>%
    mutate(genus = genus)
}

all_results <- map_dfr(rds_files, run_model)
#saveRDS(all_results,"/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/Squirrels/masting_project/BRMS_data_again/all_results.RDS" )
all_results = readRDS("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/Squirrels/masting_project/BRMS_data_again/all_results.RDS" )






all_results1 = as.data.frame(all_results)

all_results1 <- all_results1 %>%
  mutate(sig = case_when(
    conf.low > 0 & abs(estimate) > 0.05  ~ "positive",
    conf.high < 0 & abs(estimate) > 0.05 ~ "negative",
    TRUE                                  ~ "ns"
  ))

ggplot(all_results1, aes(x = genus, y = estimate, color = sig)) +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
  scale_color_manual(
    values = c("positive" = "steelblue", "negative" = "firebrick", "ns" = "grey60"),
    labels = c("positive" = "CI > 0", "negative" = "CI < 0", "ns" = "ns"),
    name = NULL
  ) +
  facet_wrap(~ term, ncol = 1, scales = "free_y") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  ) +
  labs(x = "Genus", y = "Estimate", title = "Model estimates by genus and term")
allres


unique(all_results1$term)

all_results1 = all_results1 %>%
  filter(sig != "ns", term != "year_dist") %>%
  mutate(term = recode(term,
                       "mother_offspring"               = "MO Pairs",
                       "masting"                        = "Masting",
                       "same_season"                    = "Same season",
                       "same_year"                      = "Same year",
                       "dist_m"                         = "Geo dist",
                       "same_sex"                       = "Same sex",
                       "run"                            = "Run",
                       "age_dist"                       = "Age dist",
                       "mother_offspring × masting"     = "MO x mast",
                       "mother_offspring × same_season" = "MO x season"
  ))



  a= ggplot(all_results1, aes(x = genus, y = estimate, color = sig)) +
  geom_point() +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("positive" = "steelblue", "negative" = "firebrick")) +
  facet_wrap(~ term, scales = "free_x", nrow = 1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "Genus", y = "Difference in Similarity", color = NULL,
       title = "Jaccard")
### can I only also plot genera with signficant results
a





### BC similairty ---------
model_dir <- "/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/Squirrels/masting_project/BRMS_data_again/model_results_genera_BC/"

rds_files <- list.files(model_dir, pattern = "\\.rds$", full.names = TRUE)
rds_files[1:10]
vars <- c("mother_offspring", "masting", "same_season", 
          "year_dist", "same_year", "dist_m", "same_sex", 
          "run", "age_dist")

run_model <- function(filepath) {
  
  genus <- gsub("brms_model_|\\.rds$", "", basename(filepath))
  message("Processing: ", genus)
  
  mmmm <- readRDS(filepath)
  
  newdata_grid <- datagrid(
    model            = mmmm,
    mother_offspring = c(0, 1),
    masting          = c(0, 1),
    same_season      = mean(mmmm$data$same_season),
    year_dist        = mean(mmmm$data$year_dist),
    same_year        = mean(mmmm$data$same_year),
    dist_m           = mean(mmmm$data$dist_m),
    same_sex         = mean(mmmm$data$same_sex),
    run              = mean(mmmm$data$run),
    age_dist         = mean(mmmm$data$age_dist)
  )
  
  main1 <- lapply(vars, function(v) {
    avg_slopes(mmmm,
               variables  = v,
               newdata    = newdata_grid,
               re_formula = NA,
               ndraws     = 500)
  }) %>% bind_rows() %>%
    select(term, estimate, conf.low, conf.high)
  
  int1 <- avg_comparisons(mmmm,
                          variables  = list(mother_offspring = c(0,1), masting = c(0,1)),
                          cross      = TRUE,
                          re_formula = NA,
                          ndraws     = 500) %>%
    mutate(term = "mother_offspring × masting") %>%
    select(term, estimate, conf.low, conf.high)
  
  int2 <- avg_comparisons(mmmm,
                          variables  = list(mother_offspring = c(0,1), same_season = c(0,1)),
                          cross      = TRUE,
                          re_formula = NA,
                          ndraws     = 500) %>%
    mutate(term = "mother_offspring × same_season") %>%
    select(term, estimate, conf.low, conf.high)
  
  bind_rows(main1, int1, int2) %>%
    mutate(genus = genus)
}

all_results <- map_dfr(rds_files, run_model)
saveRDS(all_results,"/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/Squirrels/masting_project/BRMS_data_again/all_results_BC.RDS" )
all_results = readRDS("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/Squirrels/masting_project/BRMS_data_again/all_results_BC.RDS" )

all_results_BC = as.data.frame(all_results)

all_results_BC <- all_results_BC %>%
  mutate(sig = case_when(
    conf.low > 0 & abs(estimate) > 0.05  ~ "positive",
    conf.high < 0 & abs(estimate) > 0.05 ~ "negative",
    TRUE                                  ~ "ns"
  ))

ggplot(all_results_BC, aes(x = genus, y = estimate, color = sig)) +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
  scale_color_manual(
    values = c("positive" = "steelblue", "negative" = "firebrick", "ns" = "grey60"),
    labels = c("positive" = "CI > 0", "negative" = "CI < 0", "ns" = "ns"),
    name = NULL
  ) +
  facet_wrap(~ term, ncol = 1, scales = "free_y") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  ) +
  labs(x = "Genus", y = "Estimate", title = "Model estimates by genus and term")
allres


unique(all_results1$term)

all_results_BC1 = all_results_BC %>%
  filter(sig != "ns", term != "year_dist") %>%
  mutate(term = recode(term,
                       "mother_offspring"               = "MO Pairs",
                       "masting"                        = "Masting",
                       "same_season"                    = "Same season",
                       "same_year"                      = "Same year",
                       "dist_m"                         = "Geo dist",
                       "same_sex"                       = "Same sex",
                       "run"                            = "Run",
                       "age_dist"                       = "Age dist",
                       "mother_offspring × masting"     = "MO x mast",
                       "mother_offspring × same_season" = "MO x season"
  ))
all_results_BC1 <- all_results_BC1[all_results_BC1$genus != "BC_NA", ]
all_results_BC1$genus = gsub("BC_","",all_results_BC1$genus)
  b = ggplot(all_results_BC1, aes(x = genus, y = estimate, color = sig)) +
  geom_point() +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("positive" = "steelblue", "negative" = "firebrick")) +
  facet_wrap(~ term, scales = "free_x", nrow = 1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "Genus",y = "Difference in Similiarity", color = NULL,
       title = "Bray-Curtis")

  b
  ### can I only also plot genera with signficant results
  
  library(ggpubr)
  ggarrange(b,a, ncol = 1, common.legend = TRUE)
  