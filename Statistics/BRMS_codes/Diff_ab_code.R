

#### diff ab -----------------
table_table = generate.tax.summary.modified(ASVtab_16s,as.data.frame(Taxtable))

#genus_tab = data.frame(resss$ASVtab)
genus_tab = data.frame(t(table_table$tax7))
genus_tab_rel = decostand(genus_tab, method = "total")
rowSums(genus_tab_rel)
genus_tab1 = colSums(genus_tab_rel)
list_NAMES = names(sort(genus_tab1, decreasing = TRUE)[1:50])
saveRDS(list_NAMES,"list_top50_genera.RDS")

list_NAMES
results_list <- list()
emmeans_list <- list()
for (asv in list_NAMES) {
  #asv = "Ruminococcus"
  metadata_table$genus <- genus_tab[[asv]]
  result <- tryCatch({
    model <- glmmTMB(
      genus ~  season_cat*mast + (1 | squirrel_id_letter) + (1|year_factor)  + (1|run.x) + offset(log(read_number)), # change here
      data = metadata_table,
      family = nbinom2
    )
    
    # Save fixed effect summary
    tidy_out <- broom.mixed::tidy(model)
    tidy_out$ASV <- asv
    
    # get emmeans for mast within each season
    emm <- emmeans(model, ~ mast | season_cat)
    
    # get mast vs non-mast contrasts within each season
    contrasts <- pairs(emm) %>%        # same as contrast(emm, method = "pairwise")
      as.data.frame() %>%
      mutate(ASV = asv)
    
    # <<--- ADD THESE LINES: convert from ln-fold-change to log2-fold-change --->
    contrasts <- contrasts %>%
      mutate(
        estimate_log2 = estimate / log(2),    # log2 fold change
        SE_log2      = SE / log(2)            # standard error on log2 scale
      )
    
    
    
    
    list(tidy = tidy_out, emm = contrasts)
    
  }, error = function(e) {
    message(paste("Model failed for", asv, ":", e$message))
    list(
      tidy = data.frame(term = NA, estimate = NA, std.error = NA, statistic = NA, p.value = NA, ASV = asv),
      emm = data.frame(season_cat = NA, contrast = NA, estimate = NA, SE = NA, df = NA, z.ratio = NA, p.value = NA, ASV = asv)
    )
  })
  
  results_list[[asv]] <- result$tidy
  emmeans_list[[asv]] <- result$emm
}

# Combine
results_df  <- bind_rows(results_list)
emmeans_df  <- bind_rows(emmeans_list)
emmeans_df$qvalue = p.adjust(emmeans_df$p.value, method = "fdr")
alpha <- 0.05

df_emm_plot <- emmeans_df %>%
  #  filter(!is.na(season_cat)) %>%
  mutate(
    signif = qvalue < alpha,
    show_coef = sprintf("%.2f", estimate_log2),
    fill_col = ifelse(signif, estimate_log2, NA)
  )

a <- ggplot(df_emm_plot, aes(x = season_cat, y = ASV, fill = fill_col)) +
  geom_tile(color = "grey90") +
  geom_text(aes(label = show_coef), size = 3) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-2, 2),
    breaks = c(-2, 0, 2),
    labels = c("< -2", "0", "> 2"),
    oob = scales::squish,
    na.value = "white",
    name = "Log(fold change)"
  ) +
  labs(
    title = "Mast across Seasons",
    x = "Season",
    y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ theme(axis.text.y = element_blank())

a





#### diff ab
#table_table = generate.tax.summary.modified(ASVtab_16s.rar,as.data.frame(Taxtable))
results_list <- list()
emmeans_list <- list()

for (asv in list_NAMES) {
  # add genus abundance to metadata
  metadata_table$genus <- genus_tab[[asv]]
  
  result <- tryCatch({
    # --- MODEL: only test for seasonal effects ---
    model <- glmmTMB(
      genus ~ season_cat + (1 | squirrel_id) + (1 | year_factor) + (1 | run.x)+ offset(log(read_number)),
      data = metadata_table,
      family = nbinom2
    )
    
    # extract tidy summary of fixed effects
    tidy_out <- broom.mixed::tidy(model)
    tidy_out$ASV <- asv
    
    # get estimated marginal means for each season
    emm <- emmeans(model, ~ season_cat)
    
    # pairwise contrasts between seasons
    contrasts <- pairs(emm) %>%
      as.data.frame() %>%
      mutate(
        ASV = asv,
        estimate_log2 = estimate / log(2),  # log2 fold change
        SE_log2 = SE / log(2)
      )
    
    list(tidy = tidy_out, emm = contrasts)
    
  }, error = function(e) {
    message(paste("Model failed for", asv, ":", e$message))
    list(
      tidy = data.frame(term = NA, estimate = NA, std.error = NA,
                        statistic = NA, p.value = NA, ASV = asv),
      emm = data.frame(season_cat = NA, contrast = NA, estimate = NA,
                       SE = NA, df = NA, z.ratio = NA, p.value = NA, ASV = asv)
    )
  })
  
  results_list[[asv]] <- result$tidy
  emmeans_list[[asv]] <- result$emm
}

# combine results
results_df <- bind_rows(results_list)
emmeans_df <- bind_rows(emmeans_list)

emmeans_df$qvalue = p.adjust(emmeans_df$p.value, method = "fdr")
alpha <- 0.05

df_emm_plot <- emmeans_df %>%
  #  filter(!is.na(season_cat)) %>%
  mutate(
    signif = qvalue < alpha,
    show_coef = sprintf("%.2f", estimate_log2),
    fill_col = ifelse(signif, estimate_log2, NA)
  )

b= ggplot(df_emm_plot, aes(x = contrast, y = ASV, fill = fill_col)) +
  geom_tile(color = "grey90") +
  geom_text(aes(label = show_coef), size = 3) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-2, 2),
    breaks = c(-2, 0, 2),
    labels = c("< -2", "0", "> 2"),
    oob = scales::squish,
    na.value = "white",
    name = "Log(fold change)"
  ) +
  labs(
    title = "Seasonal differences",
    x = "Season",
    y = "Genus"
  ) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

b

ggarrange(b,a, common.legend = TRUE, align = "h", 
          widths  = c(1.9, 1))


### example plots
#### Plot 
###boxplot
samdata_subset1 = metadata_table
list_NAMES
cor.test(genus_tab[["Prevotellaceae.UCG.001"]] ,genus_tab[["Bacteroides"]] )
samdata_subset1$Genus = genus_tab[["Ruminococcus"]]
#samdata_subset1$season_cat = droplevels(samdata_subset1$food_avail)
ggplot(samdata_subset1, aes(year_factor,Genus, fill =season_cat , color = season_cat)) +  
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, alpha = 0.7) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), 
              size = 1, alpha = 0.6) +
  theme_classic() +
  labs(y = "Counts", x = "Season", fill = "Year", title = "Prevotella") +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_x_discrete(labels = c(
    "Early Spring" = "Day < 120",
    "Late Spring" = "120 < Day <  180",
    "Summer" = "Day > 180"
  ))

samdata_subset2019_2022 = samdata_subset1[samdata_subset1$years_cat %in% c("2019", "2022"),]
genus_tab_2019_2022 = genus_tab[rownames(genus_tab) %in% samdata_subset2019_2022$name,]
sort(colSums(genus_tab_2019_2022), decreasing = TRUE)[1:10]
genus_tab_others = genus_tab[rownames(genus_tab) %nin% samdata_subset2019_2022$name,]
sort(colSums(genus_tab_others), decreasing = TRUE)[1:10]
sort(colSums(genus_tab), decreasing = TRUE)[1:10]


# Step 0: Fit the model
model <- glmmTMB(
  Genus ~  season_cat*mast + (1 | squirrel_id_letter) + (1|year_factor)  + (1|run.x),
  data = samdata_subset1,
  family = nbinom2
)
summary(model)

# Step 1: Find all squirrels with multiple seasons within the same year
squirrels_multiple_seasons <- samdata_subset1 %>%
  group_by(squirrel_id_letter, year_factor) %>%
  summarise(
    n_seasons = n_distinct(season_cat),
    n_obs = n(),
    .groups = "drop"
  ) %>%
  filter(n_seasons >= 2)

# Step 2: Filter your data to only these squirrel-year combinations
data_to_plot <- samdata_subset1 %>%
  semi_join(squirrels_multiple_seasons, by = c("squirrel_id_letter", "year_factor"))

# Step 3: Get predicted trends with confidence intervals using emmeans
library(emmeans)

emm <- emmeans(model, ~ season_cat | mast, type = "response")
pred_trend <- as.data.frame(emm)

# Check what columns we have
print(names(pred_trend))
print(head(pred_trend))

# Step 4: Plot (I'll adjust the column names based on what we see)
b = data_to_plot %>%
  mutate(squirrel_year = paste(squirrel_id_letter, year_factor, sep = " | ")) %>%
  ggplot(aes(x = season_cat, y = Genus, group = squirrel_year, color = year_factor)) +
  geom_point(size = 2, alpha = 0.2) +
  geom_line(linewidth = 0.8, alpha = 0.2) +
  # Add confidence ribbon - using different possible column names
  geom_ribbon(data = pred_trend, 
              aes(x = season_cat, ymin = asymp.LCL, ymax = asymp.UCL, group = 1),
              fill = "gray30", alpha = 0.3, inherit.aes = FALSE) +
  # Add bold prediction line
  geom_line(data = pred_trend, 
            aes(x = season_cat, y = response, group = 1), 
            color = "black", linewidth = 1.5) +
  facet_wrap(~mast) +
  theme_bw() +
  labs(title = "Ruminococcus",
       #subtitle = "Thin lines = individual squirrel-years | Bold line = predicted trend | Gray ribbon = 95% CI",
       x = "Season", 
       y = "Rarefied genus abundance",
       color = "Year")
b
ggarrange(a,b,c, ncol = 1, common.legend = TRUE)
