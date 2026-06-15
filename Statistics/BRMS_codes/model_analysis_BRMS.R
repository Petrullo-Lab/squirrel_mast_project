library(emmeans)
library(tidybayes)
library(ggplot2)
library(brms)
library(tibble)
library(dplyr)


#### Only year distance august ------
model_bc =readRDS("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/Squirrels/masting_project/BRMS_data_again/rarefied_dataset/model1_BC_rarefied.RDS")
summary(model_bc)


post <- as_draws_df(model_bc)

# look at posterior correlations
cor(post[, c(
  "b_same_year",
  "b_year_dist",
  "b_masting",
  "b_run",
  "b_mother_offspring:masting",
  "b_mother_offspring",
  "b_same_season"
)])

colnames(post)

# extract fixed effects summary and tidy column names
fe <- as.data.frame(brms::fixef(model_bc))
fe <- rownames_to_column(fe, var = "term")
fe
fe <- fe %>% dplyr::rename(estimate = Estimate,
                           se = Est.Error,
                           lower = `Q2.5`,
                           upper = `Q97.5`)
# reorder terms for plotting (largest -> smallest)
fe$term <- factor(fe$term, levels = rev(fe$term))

fe <- fe %>%
  mutate(sig = ifelse(lower > 0 | upper < 0, "significant", "not sigificant"))




a = ggplot(fe, aes(x = term, y = estimate, color = sig)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25) +
  coord_flip() +
  scale_color_manual(values = c("significant" = "red", "ns" = "black")) +
  labs(x = NULL, y = "Estimate (logit scale)", color = "") +
  theme_classic(base_size = 13)+ ggtitle("Bray Curtis")

a
performance::check_collinearity(model_bc)

# Get EMMs for the interaction
emm <- emmeans(model_bc, ~ mother_offspring * same_season)
# Compare masting levels (within each mother_offspring group)
pairs(emm, by = "same_season")
pairs(emm, by = "same_season",level = 0.99)  # For **
pairs(emm,by = "same_season", level = 0.999) # For ***
emm_masting <- emmeans(model_bc, ~ same_season, type = "response")
pairs(emm_masting)
pairs(emm_masting, level = 0.99)  # For **
pairs(emm_masting, level = 0.999) # For ***

# Get EMMs for the interaction
emm <- emmeans(model_bc, ~ mother_offspring * masting)
# Compare masting levels (within each mother_offspring group)
pairs(emm, by = "masting")
pairs(emm, by = "masting",level = 0.99)  # For **
pairs(emm,by = "masting", level = 0.999) # For ***
emm_masting <- emmeans(model_bc, ~ masting, type = "response")
pairs(emm_masting)
pairs(emm_masting, level = 0.99)  # For **
pairs(emm_masting, level = 0.999) # For ***






fe <- as.data.frame(fixef(model_bc)) %>%
  rownames_to_column("term") %>%
  rename(estimate = Estimate,
         lower = `Q2.5`,
         upper = `Q97.5`) %>%
  mutate(
    intercept = estimate[term == "Intercept"],
    # transform from intercept to intercept+β on logit scale
    delta = plogis(intercept + estimate) - plogis(intercept),
    delta_low = plogis(intercept + lower) - plogis(intercept),
    delta_high = plogis(intercept + upper) - plogis(intercept),
    term = factor(term, levels = rev(term))
  )

fe <- fe %>%
  filter(!term %in% c("Intercept", "read_dist", "run")) %>%
  mutate(term = recode(term,
                       "dist_m" = "Spatial dist (m)",
                       "julian_day_dist" = "Dist Julian day",
                       "mother_offspring" = "Mother–offspring pairs",
                       "age_dist" = "Age distance (years)",
                       "masting_factor2_mast" = "Mast year",
                       "masting_factor3_mixed" = "Mixed mast/non-mast",
                       "same_season" = "Same season",
                       "same_year" = "Same year",
                       "year_dist" = "Distance in years",
                       "same_sex" = "Same sex",
                       "mother_offspring:masting_factor2_mast" = "Mother–offspring × Mast",
                       "mother_offspring:masting_factor3_mixed" = "Mother–offspring × Mixed Mast",
                       "mother_offspring:same_season" = "Mother–offspring × same Season"
  ))


fe <- fe %>%
  mutate(direction = case_when(
    delta_low > 0  ~ "increase",
    delta_high < 0 ~ "decrease",
    TRUE           ~ "not sigificant"
  ))



bc = ggplot(fe, aes(x = term, y = delta, color = direction)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = delta_low, ymax = delta_high), width = 0.25) +
  coord_flip() +
  theme_classic() +
  scale_color_manual(values = c(
    "increase" = "red",
    "decrease" = "blue",
    "ns" = "black"
  )) +
  labs(x = NULL,
       y = "Difference in BC Similarity",
       title = "Bray-Curtis",
       color = "")

bc



post <- model_bc %>%
  gather_draws(b_mother_offspring, 
               b_masting, 
               b_same_season,
               `b_mother_offspring:masting`,
               `b_mother_offspring:same_season`) %>%
  ungroup() %>%
  mutate(
    # Transform to response scale (probability difference)
    delta = plogis(.value) - plogis(0),
    # Clean up term names
    .variable = recode(.variable,
                       "b_mother_offspring" = "M/O pairs",
                       "b_masting" = "Mast year",
                       "b_same_season" = "Same season",
                       "b_mother_offspring:masting" = "M-O × Mast",
                       "b_mother_offspring:same_season"= "M-O × Same season"),
    # Set order
    .variable = factor(.variable, levels = c(
      "M/O pairs",
      "Mast year",
      "Same season",
      "M-O × Mast",
      "M-O × Same season"
    ))
  )

# Plot distributions
bc2 = ggplot(post, aes(x = delta)) +
  geom_density(fill = "grey70", color = "grey30", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~ .variable, scales = "free_y", nrow = 1) +
  theme_classic() +
  labs(x = "Difference in Bray-Curtis Similarity",
       y = "Posterior Density", title = "Bray-Curtis") +
  theme(legend.position = "none",
        strip.text = element_text(size = 9))

bc2
### level2
model_jacc = readRDS("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/Squirrels/masting_project/BRMS_data_again/rarefied_dataset/model1_jacc_rarefied.RDS")
summary(model_jacc)

# extract fixed effects summary and tidy column names
fe <- as.data.frame(brms::fixef(model_jacc))
fe <- rownames_to_column(fe, var = "term")
fe
fe <- fe %>% dplyr::rename(estimate = Estimate,
                           se = Est.Error,
                           lower = `Q2.5`,
                           upper = `Q97.5`)
# reorder terms for plotting (largest -> smallest)
fe$term <- factor(fe$term, levels = rev(fe$term))

fe <- fe %>%
  mutate(sig = ifelse(lower > 0 | upper < 0, "significant", "not sigificant"))

b = ggplot(fe, aes(x = term, y = estimate, color = sig)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25) +
  coord_flip() +
  scale_color_manual(values = c("significant" = "red", "ns" = "black")) +
  labs(x = NULL, y = "Estimate (logit scale)", color = "") +
  theme_classic(base_size = 13) + ggtitle("Jaccard")
b
library(ggpubr)
ggarrange(a,b, common.legend = TRUE)

performance::check_collinearity(model_jacc) # collinearity ok

# Get EMMs for the interaction
emm <- emmeans(model_jacc, ~ mother_offspring * same_season)
# Compare masting levels (within each mother_offspring group)
pairs(emm, by = "same_season")
pairs(emm, by = "same_season",level = 0.99)  # For **
pairs(emm,by = "same_season", level = 0.999) # For ***
emm_masting <- emmeans(model_jacc, ~ same_season, type = "response")
pairs(emm_masting)
pairs(emm_masting, level = 0.99)  # For **
pairs(emm_masting, level = 0.999) # For ***

# Get EMMs for the interaction
emm <- emmeans(model_bc, ~ mother_offspring * masting_factor)
# Compare masting levels (within each mother_offspring group)
pairs(emm, by = "masting_factor")
pairs(emm, by = "masting_factor",level = 0.99)  # For **
pairs(emm,by = "masting_factor", level = 0.999) # For ***
emm_masting <- emmeans(model_bc, ~ masting_factor, type = "response")
pairs(emm_masting)
pairs(emm_masting, level = 0.99)  # For **
pairs(emm_masting, level = 0.999) # For ***




fe <- as.data.frame(fixef(model_jacc)) %>%
  rownames_to_column("term") %>%
  rename(estimate = Estimate,
         lower = `Q2.5`,
         upper = `Q97.5`) %>%
  mutate(
    intercept = estimate[term == "Intercept"],
    # transform from intercept to intercept+β on logit scale
    delta = plogis(intercept + estimate) - plogis(intercept),
    delta_low = plogis(intercept + lower) - plogis(intercept),
    delta_high = plogis(intercept + upper) - plogis(intercept),
    term = factor(term, levels = rev(term))
  )



fe <- fe %>%
  filter(!term %in% c("Intercept", "read_dist", "run")) %>%
  mutate(term = recode(term,
                       "dist_m" = "Spatial dist (m)",
                       "julian_day_dist" = "Dist Julian day",
                       "mother_offspring" = "Mother–offspring pairs",
                       "age_dist" = "Age distance (years)",
                       "masting_factor2_mast" = "Mast year",
                       "masting_factor3_mixed" = "Mixed mast/non-mast",
                       "same_season" = "Same season",
                       "same_year" = "Same year",
                       "year_dist" = "Distance in years",
                       "same_sex" = "Same sex",
                       "mother_offspring:masting_factor2_mast" = "Mother–offspring × Mast",
                       "mother_offspring:masting_factor3_mixed" = "Mother–offspring × Mixed Mast",
                       "mother_offspring:same_season" = "Mother–offspring × same Season"
  ))


fe <- fe %>%
  mutate(direction = case_when(
    delta_low > 0  ~ "increase",
    delta_high < 0 ~ "decrease",
    TRUE           ~ "not sigificant"
  ))





fe <- fe %>%
  mutate(direction = case_when(
    delta_low > 0  ~ "increase",
    delta_high < 0 ~ "decrease",
    TRUE           ~ "not sigificant"
  ))


jacc = ggplot(fe, aes(x = term, y = delta, color = direction)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = delta_low, ymax = delta_high), width = 0.25) +
  coord_flip() +
  theme_classic() +
  scale_color_manual(values = c(
    "increase" = "red",
    "decrease" = "blue",
    "ns" = "black"
  )) +
  labs(x = NULL,
       y = "Difference in Jaccard Similarity",
       title = "Jaccard",
       color = "")

jacc



post1 <- model_jacc %>%
  gather_draws(b_mother_offspring, 
               b_masting_factor2_mast, 
               b_masting_factor3_mixed,
               b_same_season,
               `b_mother_offspring:masting_factor2_mast`,
               `b_mother_offspring:masting_factor3_mixed`,
               `b_mother_offspring:same_season`)%>%
  ungroup() %>%
  mutate(
    # Transform to response scale (probability difference)
    delta = plogis(.value) - plogis(0),
    # Clean up term names
    .variable = recode(.variable,
                       "b_mother_offspring" = "M/O pairs",
                       "b_masting_factor2_mast" = "Mast year",
                       "b_masting_factor3_mixed" = "Mixed mast",
                       "b_same_season" = "Same season",
                       "b_mother_offspring:masting_factor2_mast" = "M-O × Mast",
                       "b_mother_offspring:masting_factor3_mixed" = "M-O × mixed Mast",
                       "b_mother_offspring:same_season"= "M-O × Same season"),
    # Set order
    .variable = factor(.variable, levels = c(
      "M/O pairs",
      "Mast year",
      "Mixed mast",
      "Same season",
      "M-O × Mast",
      "M-O × mixed Mast",
      "M-O × Same season"
    ))
  )
# Plot distributions
jacc2 = ggplot(post1, aes(x = delta)) +
  geom_density(fill = "grey70", color = "grey30", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~ .variable, scales = "free_y", nrow = 1) +
  theme_classic() +
  labs(x = "Difference in Bray-Curtis Similarity",
       y = "Posterior Density", title = "Jaccard") +
  theme(legend.position = "none",
        strip.text = element_text(size = 9))

jacc2

ggarrange(bc, jacc, ncol = 2, common.legend = TRUE)
ggarrange(bc2, jacc2, ncol = 1)



