### start with analysis
### libraries
library(vegan)
library(ggplot2)
library(dplyr)
library(lubridate)
library(ggpubr)
library(rcompanion)
library(lmerTest)
library(glmmTMB)
library(emmeans)
library(splines)
library(lomb )
library(dada2)

source("~/Library/CloudStorage/OneDrive-UniversityofArizona/helpful_r_functions.R")

### read tabales in
setwd("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/Squirrels/analysis_data_oct_2025/16s/comb_ASV_tabs")
ASVtable = readRDS("new_ASVtable_all_samples.RDS")
Taxtable = readRDS("Taxtab_total_silva.RDS")
metadata_table = readRDS("new_meta_all_samples_bdate.RDS")



unique(rownames(ASVtable)  == metadata_table$name) # doh
metadata_table <- metadata_table[match(rownames(ASVtable), metadata_table$name), ]
unique(rownames(ASVtable)  == metadata_table$name) # yay
unique(rownames(Taxtable)  == colnames(ASVtable)) # yay
sort(rowSums(ASVtable))
ncol(ASVtable)




mypalette2  <-c("#40004b","#ffffbf","#762a83","#de77ae","#f46d43","#f7f7f7","#d9f0d3","#a6dba0","#a6cee3","#1f78b4","#b2df8a","#33a02c","#4d9221",
                "#2166ac","#5aae61","#1b7837","#d73027","#00441b","#543005","#8c510a","#bf812d","#dfc27d","#f6e8c3","#b2182b","#de77ae","#bc80bd","#f5f5f5","#c7eae5","#80cdc1","#35978f","#003c30","#8e0152","#fb8072","#c51b7d","#f1b6da","#fde0ef","#fdb462","#b3de69","#7fbc41","#276419","#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#ff7f00","#cab2d6","#6a3d9a","#ffff99")



taxaDTDB.sum = generate.tax.summary.modified(data.frame(ASVtable), data.frame(Taxtable)) # generates list
tax_phyla = taxaDTDB.sum$tax3
tax_phyla_rel = decostand(t(tax_phyla), method = "total")
sort(colMeans(tax_phyla_rel))
tax_phyla_rel = data.frame(t(tax_phyla_rel))
top10 = names(sort(rowSums(tax_phyla_rel), decreasing = TRUE)[1:10])#change this number here to change the number of groups displayed
tax_phyla_rel$phylum = rownames(tax_phyla_rel)
ASVtab_rel_melt = reshape2::melt(tax_phyla_rel)
ASVtab_rel_melt$phylum = ifelse(ASVtab_rel_melt$phylum %in% top10, ASVtab_rel_melt$phylum, "others")
metadata_table$name1 = gsub("-",".",metadata_table$name)

# Stacked + percent
# name
ASVtab_rel_melt$day  = as.factor(metadata_table$day_of_year[match(ASVtab_rel_melt$variable,metadata_table$name1 )])
ASVtab_rel_melt$year  = metadata_table$year[match(ASVtab_rel_melt$variable,metadata_table$name1 )]
ASVtab_rel_melt$day <- droplevels(ASVtab_rel_melt$day)

ggplot(ASVtab_rel_melt, aes(fill = phylum, y = value, x = day)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = mypalette2) +
  #facet_grid(. ~ year, scales = "free_x", space = "free") +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  )

tax_fam = taxaDTDB.sum$tax6
tax_fam_rel = decostand(t(tax_fam), method = "total")
sort(colMeans(tax_fam_rel))
tax_fam_rel = data.frame(t(tax_fam_rel))
top10 = names(sort(rowSums(tax_fam_rel), decreasing = TRUE)[1:25])#change this number here to change the number of groups displayed
tax_fam_rel$family = rownames(tax_fam_rel)
ASVtab_rel_melt = reshape2::melt(tax_fam_rel)
ASVtab_rel_melt$family = ifelse(ASVtab_rel_melt$family %in% top10, ASVtab_rel_melt$family, "others")
metadata_table$name1 = gsub("-",".",metadata_table$name)

# Stacked + percent
# name
ASVtab_rel_melt$day  = as.factor(metadata_table$day_of_year[match(ASVtab_rel_melt$variable,metadata_table$name1 )])
ASVtab_rel_melt$year  = metadata_table$year[match(ASVtab_rel_melt$variable,metadata_table$name1 )]
ASVtab_rel_melt$day <- droplevels(ASVtab_rel_melt$day)

ggplot(ASVtab_rel_melt, aes(fill = family, y = value, x = day)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = mypalette2) +
  #facet_grid(. ~ year, scales = "free_x", space = "free") +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  )




metadata_table = metadata_table[metadata_table$day_of_year < 225,]
nrow(metadata_table)
### put year and mast
metadata_table$year_factor = as.factor(metadata_table$year)
metadata_table$mast = ifelse(metadata_table$year %in% c(2010, 2014, 2019, 2022), "mast", "non-mast")
table(metadata_table$gr)
metadata_table = metadata_table[metadata_table$gr != "BT",]
ASVtable = ASVtable[rownames(ASVtable) %in% metadata_table$name,]



metadata_table$season_cat <- cut(
  metadata_table$day_of_year,
  breaks = c(min(metadata_table$day_of_year)-1, 120, 180, max(metadata_table$day_of_year)+1),
  labels = c("Early Spring", "Late Spring", "Summer"),
  right = FALSE # so 60 is included in Early Spring, 120 in Late Spring, etc.
)


min(metadata_table$day_of_year)
max(metadata_table$day_of_year)


# Count number of samples per masting status, season, and year
counts <- metadata_table %>%
  count(mast, season_cat, year_factor)

# Plot
ggplot(counts, aes(x = year_factor, y = n, fill = year_factor)) +
  geom_col() +
  facet_grid(mast ~ season_cat) +
  labs(
    x = "Year",
    y = "Number of samples",
    title = "Sample counts by season, year, and masting status"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "bold")
  )


season_bg <- data.frame(
  ymin = c(0, 120, 180),
  ymax = c(120, 180, Inf),
  season = c("Early Spring", "Late Spring", "Summer")
)

ggplot(metadata_table,
       aes(x = factor(year),
           y = day_of_year,
           color = mast)) +
  
  geom_rect(data = season_bg,
            aes(ymin = ymin,
                ymax = ymax,
                xmin = -Inf,
                xmax = Inf,
                fill = season),
            inherit.aes = FALSE,
            alpha = 0.15) +
  
  geom_jitter(width = 0.15,
              alpha = 0.7,
              size = 2) +
  
  scale_fill_manual(values = c(
    "Early Spring" = "steelblue",
    "Late Spring" = "orange",
    "Summer" = "firebrick"
  )) +
  
  labs(
    x = "Year",
    y = "Day of year sampled",
    color = "Mast status",
    fill = "Season"
  ) +
  
  theme_bw()



### alpha diversity -------------
sort(rowSums(ASVtable))
metadata_table$read_number = rowSums(ASVtable)
metadata_table$Richness<-specnumber(ASVtable)
metadata_table$Shannon<-vegan::diversity(ASVtable, index="shannon")#Chang
metadata_table$Richness.rar<-specnumber(rrarefy(ASVtable, 10000)) #Change this number according to the lowest sample
metadata_table$Shannon.rar<-vegan::diversity(rrarefy(ASVtable, 10000), index="shannon")#Change this number according to the lowest sample
metadata_table = metadata_table[metadata_table$Shannon.rar >3,]
#saveRDS(metadata_table,"metadata_table_sub.RDS")
#saveRDS(ASVtable,"ASVtable_sub.RDS")



metadata_table$tukey_shannon = transformTukey(metadata_table$Shannon)
metadata_table$squirrel_id_letter = paste("a_",metadata_table$squirrel_id, sep = "")
metadata_table$day_scaled <- scale(metadata_table$day_of_year, center = TRUE, scale = TRUE)



a <- ggplot(metadata_table, aes(x = season_cat, y = Richness.rar, fill = mast, color = mast)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, alpha = 0.7) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), 
              size = 1, alpha = 0.6) +
  theme_classic() +
  labs(y = "Number of ASVs", x = NULL, fill = "Year") +
  scale_fill_manual(values = c("#1b9e77", "#7570b3")) +
  scale_color_manual(values = c("#1b9e77", "#7570b3")) +
  theme(axis.text.x = element_blank()) + guides(color = "none")

a
b = ggplot(metadata_table, aes(x = season_cat, y = Shannon.rar, fill = mast, color = mast)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, alpha = 0.7) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), 
              size = 1, alpha = 0.6) +
  theme_classic() +
  labs(y = "Shannon H'", x = "Season", fill = "Year") +
  scale_fill_manual(values = c("#1b9e77", "#7570b3")) +
  scale_color_manual(values = c("#1b9e77", "#7570b3")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_x_discrete(labels = c(
    "Early Spring" = "Early Spring",
    "Late Spring" = "Late Spring",
    "Summer" = "Summer"
  ))+ guides(color = "none")

b

ggarrange(a,b, common.legend = TRUE, ncol =1, legend = "top", heights = c(1, 1.5))



shannon_season_no_slopes = lmer(Richness ~ season_cat*mast+ age+sex+gr+  read_number+
                               (1 | squirrel_id_letter)+ (1|year_factor) +(1|run.x) , data = metadata_table_nores)


residuals <- residuals(shannon_season_no_slopes)
qqnorm(residuals) # o

model_formula <- formula(shannon_season_no_slopes)
model_formula
anova(shannon_season_no_slopes)

#tab_model(shannon_season_no_slopes,
#          #show.anova = TRUE,
#          file = "/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/Squirrels/masting_project/output_models/model_results.doc")

summary(shannon_season_no_slopes)

emm_season_mast <- emmeans(shannon_season_no_slopes, ~ season_cat)
pairs(emm_season_mast,  adjust = "tukey")
emm_season_mast <- emmeans(shannon_season_no_slopes, ~ mast|season_cat)
pairs(emm_season_mast,  adjust = "tukey")

# Average values for each squirrel-year-season combination (for ALL squirrels)
data_to_plot <- metadata_table_nores %>%
  group_by(squirrel_id_letter, year_factor, season_cat, mast) %>%
  summarise(
    Richness = mean(Richness, na.rm = TRUE),
    # Add any other variables you need
    .groups = "drop"
  )

# Get predicted trends with confidence intervals using emmeans
library(emmeans)
emm <- emmeans(shannon_season_no_slopes, ~ season_cat | mast, type = "response")
pred_trend <- as.data.frame(emm)

library(dplyr)
library(ggplot2)

a <- data_to_plot %>%
  mutate(squirrel_year = paste(squirrel_id_letter, year_factor, sep = " | ")) %>%
  ggplot(aes(x = season_cat, y = Richness, group = squirrel_year, color = year_factor)) +
  
  # thin faint points & subject-level lines (colored by year)
  geom_point(size = 2, alpha = 0.2) +
  geom_line(linewidth = 0.8, alpha = 0.2) +
  
  # ribbon (prediction interval) — uses pred_trend data; facetting will match by mast
  geom_ribbon(
    data = pred_trend,
    aes(x = season_cat, ymin = lower.CL, ymax = upper.CL, group = 1),
    fill = "gray30",
    alpha = 0.3,
    inherit.aes = FALSE
  ) +
  
  # two thicker trend lines: one for mast (orange), one for non-mast (green)
  geom_line(
    data = subset(pred_trend, mast == "mast"),
    aes(x = season_cat, y = emmean, group = 1),
    color = "#1b9e77",   # green
    linewidth = 1.5,
    inherit.aes = FALSE
  ) +
  geom_line(
    data = subset(pred_trend, mast == "non-mast"),
    aes(x = season_cat, y = emmean, group = 1),
    color = "#7570b3",   # blue/purple
    linewidth = 1.5,
    inherit.aes = FALSE
  ) +
  
  # facets, x-axis spacing, legend title, theme, and rotated x labels
  facet_wrap(~mast) +
  scale_x_discrete(expand = c(0.05, 0.05)) +
  labs(
    x = "Season",
    y = "Number of ASVs",
    color = "Year"          # legend title for the year mapping
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

# print
a


shannon_season_no_slopes = lmer(tukey_shannon ~ season_cat*mast+ age+sex+gr+  read_number+
                                  (1 | squirrel_id_letter)+ (1|year_factor) +(1|run.x) , data = metadata_table_nores)


residuals <- residuals(shannon_season_no_slopes)
qqnorm(residuals) # o
#tab_model(shannon_season_no_slopes,
#          #show.anova = TRUE,
#          file = "/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/Squirrels/masting_project/output_models/model_results_shannon.doc")


anova(shannon_season_no_slopes)
summary(shannon_season_no_slopes)

emm_season_mast <- emmeans(shannon_season_no_slopes, ~ season_cat)
pairs(emm_season_mast,  adjust = "tukey")
emm_season_mast <- emmeans(shannon_season_no_slopes, ~ mast|season_cat)
pairs(emm_season_mast,  adjust = "tukey")



# Average values for each squirrel-year-season combination (for ALL squirrels)
data_to_plot <- metadata_table_nores %>%
  group_by(squirrel_id_letter, year_factor, season_cat, mast) %>%
  summarise(
    tukey_shannon = mean(tukey_shannon, na.rm = TRUE),
    # Add any other variables you need
    .groups = "drop"
  )

# Get predicted trends with confidence intervals using emmeans
library(emmeans)
emm <- emmeans(shannon_season_no_slopes, ~ season_cat | mast, type = "response")
pred_trend <- as.data.frame(emm)

# Plot
b <- data_to_plot %>%
  mutate(squirrel_year = paste(squirrel_id_letter, year_factor, sep = " | ")) %>%
  ggplot(aes(x = season_cat, y = tukey_shannon, group = squirrel_year, color = year_factor)) +
  
  # thin faint points & subject-level lines (colored by year)
  geom_point(size = 2, alpha = 0.2) +
  geom_line(linewidth = 0.8, alpha = 0.2) +
  
  # ribbon (prediction interval) — uses pred_trend data; facetting will match by mast
  geom_ribbon(
    data = pred_trend,
    aes(x = season_cat, ymin = lower.CL, ymax = upper.CL, group = 1),
    fill = "gray30",
    alpha = 0.3,
    inherit.aes = FALSE
  ) +
  
  # two thicker trend lines: one for mast (orange), one for non-mast (green)
  geom_line(
    data = subset(pred_trend, mast == "mast"),
    aes(x = season_cat, y = emmean, group = 1),
    color = "#1b9e77",   # green
    linewidth = 1.5,
    inherit.aes = FALSE
  ) +
  geom_line(
    data = subset(pred_trend, mast == "non-mast"),
    aes(x = season_cat, y = emmean, group = 1),
    color = "#7570b3",   # blue/purple
    linewidth = 1.5,
    inherit.aes = FALSE
  ) +
  
  # facets, x-axis spacing, legend title, theme, and rotated x labels
  facet_wrap(~mast) +
  scale_x_discrete(expand = c(0.05, 0.05)) +
  labs(
    x = "Season",
    y = "Tukey Transformed Shannon H'",
    color = "Year"          # legend title for the year mapping
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

# print
b
 a = a+ theme(
   axis.text.x = element_blank(),axis.title.x = element_blank()
 )
 


a <- a + guides(color = guide_legend(nrow = 2)) +theme(
  legend.title = element_text(size = 6),
  legend.text  = element_text(size = 6),
  legend.key.size = unit(0.3, "cm"),
  legend.spacing.y = unit(0.10, "cm")
)
b <- b + guides(color = guide_legend(nrow = 2)) +theme(
  legend.title = element_text(size = 6),
  legend.text  = element_text(size = 6),
  legend.key.size = unit(0.3, "cm"),
  legend.spacing.y = unit(0.10, "cm")
)

# Then combine
ggarrange(a, b, ncol = 1, heights = c(0.7,1),common.legend = TRUE, align = "v")
 
 










 #### Beta diversity ---------
metadata_table <- metadata_table %>%
  mutate(years_cat = if_else(mast == "non-mast", "non-mast", as.character(year_factor)))
#metadata_table_non_mast = metadata_table[metadata_table$years_cat =="non-mast",]

ASVtab_16s.rar = rrarefy(ASVtable[rownames(ASVtable) %in% metadata_table$name,], 10000) 
ASVtab_16s = ASVtable[rownames(ASVtable) %in% metadata_table$name,] 

#saveRDS(ASVtab_16s.rar,"ASVtab_16s.rar.RDS")
#table_table <- generate.tax.summary.modified(ASVtab_16s.rar, as.data.frame(Taxtable))
set.seed(12)
#ASVtab_16s.rar1 = t(table_table$tax6)

bac.dist_non_rar<-vegdist(ASVtab_16s, method="bray")
bac.dist.jac<-vegdist(ASVtab_16s.rar, method="jaccard", binary = TRUE)
bac_dist_jack_not_rar <-vegdist(ASVtab_16s, method="jaccard", binary = TRUE)
bac.dist<-vegdist(ASVtab_16s.rar, method="bray")

#saveRDS(bac.dist.jac, "bac_dist_jac_rar.RDS")
#saveRDS(bac_dist_jack_not_rar, "bac_dist_jac_not_rar.RDS")
#saveRDS(bac.dist_non_rar, "bac_dist_BC_not_rar.RDS")
#saveRDS(metadata_table, "metadata_mast.RDS")



#plot(bac.dist,bac.dist.jac )

set.seed(12)
bac.nmds<-metaMDS(bac.dist, k=2, try=999)
bac.nmds$stress #0.21
stress= paste("stress = ", round(bac.nmds$stress, digits = 2))

metadata_table$Axis01<-bac.nmds$points[,1]
metadata_table$Axis02<-bac.nmds$points[,2]


metadata_table_non_mast = metadata_table[metadata_table$years_cat =="non-mast",]



library(viridis)

bac_nmds <- ggplot(metadata_table, aes(x = Axis01, y = Axis02)) +
  geom_point(aes(color =season_cat), size = 2, alpha = 0.8) +
  stat_ellipse(aes( fill = season_cat), 
               geom = "polygon", alpha = 0.5, level = 0.68, show.legend = FALSE) +
  facet_wrap(~years_cat) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1",labels = c("Early Spring", "Late Spring", "Summer")) +
  theme_classic() +
  labs(
    x = "NMDS1",
    y = "NMDS2",
    color  = "Season")
   
bac_nmds

a = ggplot(metadata_table, aes(x = season_cat, y = Axis01, color = mast)) +
  geom_boxplot() +
  scale_color_manual(values = c("mast" = "#1b9e77", 
                                "non-mast" = "#7570b3")) +
  theme_classic() +
  labs(
    x = "Season",
    y = "NMDS1",
    color = "Year",
    title = "NMDS1"
  )


b = ggplot(metadata_table, aes(x = season_cat, y = Axis02, color = mast)) +
  geom_boxplot() +
  scale_color_manual(values = c("mast" = "#1b9e77", 
                                "non-mast" = "#7570b3")) +
  theme_classic() +
  labs(
    x = "Season",
    y = "NMDS2",
    color = "Year",
    title = "NMDS2"
  )

ggarrange(a,b, common.legend = TRUE)

b <- ggplot(subset(metadata_table, season_cat == "Summer"),
            aes(x = mast, y = Axis01, fill = mast)) +
  
  # Violin
  geom_violin(trim = FALSE, alpha = 0.6) +
  
  # Boxplot inside violin
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8) +
  
  # Jitter points (same color as fill)
  geom_jitter(aes(color = mast),
              width = 0.1, size = 1.5, alpha = 0.6) +
  
  scale_fill_manual(values = c("mast" = "#1b9e77", 
                               "non-mast" = "#7570b3")) +
  scale_color_manual(values = c("mast" = "#1b9e77", 
                                "non-mast" = "#7570b3")) +
  
  theme_classic() +
  labs(
    x = NULL,
    y = "NMDS1",
    fill = "Year",
    color = "Year",
    title = "NMDS1 (Summer only)"
  ) +
  
  # Remove duplicate legend
  guides(color = "none")
b

set.seed(12)
permanova_res = adonis2(bac.dist~run.x+sex+gr+age+year_factor+mast+season_cat*mast ,  
                        by = "terms",data = metadata_table, strata = metadata_table$squirrel_id_letter)
permanova_res


metadata_table1 <- metadata_table1 %>%
  # ensure columns exist and classes are ok
  mutate(
    date = as.Date(date), 
    squirrel_id_letter = as.character(squirrel_id_letter),
    year = as.integer(year),           # if year is not integer already
    season_cat = as.character(season_cat)
  ) %>%
  # counts within same (squirrel, year, season)
  group_by(squirrel_id_letter, year, season_cat) %>%
  mutate(
    n_yr_season = n()                      # includes the focal sample
  ) %>%
  ungroup() %>%
  
  # counts within same (squirrel, year) regardless of season
  group_by(squirrel_id_letter, year) %>%
  mutate(
    n_year = n()
  ) %>%
  ungroup() %>%
  
  # counts within same (squirrel, season) across years
  group_by(squirrel_id_letter, season_cat) %>%
  mutate(
    n_season = n()
  ) %>%
  ungroup() %>%
  
  # add "other" counts (exclude the focal sample)
  mutate(
    n_other_yr_season = if_else(is.na(n_yr_season), NA_integer_, pmax(0L, n_yr_season - 1L)),
    n_other_year       = if_else(is.na(n_year),       NA_integer_, pmax(0L, n_year - 1L)),
    n_other_season     = if_else(is.na(n_season),     NA_integer_, pmax(0L, n_season - 1L))
  )

#### diff ab -----------------
table_table = generate.tax.summary.modified(ASVtab_16s.rar,as.data.frame(Taxtable))

#genus_tab = data.frame(resss$ASVtab)
genus_tab = data.frame(t(table_table$tax7))
genus_tab_rel = decostand(genus_tab, method = "total")
rowSums(genus_tab_rel)
genus_tab1 = colSums(genus_tab_rel)
list_NAMES = names(sort(genus_tab1, decreasing = TRUE)[1:50])
#saveRDS(list_NAMES,"list_top50_genera.RDS")

list_NAMES
results_list <- list()
emmeans_list <- list()
for (asv in list_NAMES) {
  #asv = "Ruminococcus"
  metadata_table$genus <- genus_tab[[asv]]
  result <- tryCatch({
    model <- glmmTMB(
      genus ~  season_cat*mast + (1 | squirrel_id_letter) + (1|year_factor)  + (1|run.x) , # change here
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
      genus ~ season_cat + (1 | squirrel_id) + (1 | year_factor) + (1 | run.x),
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
library(stringr)

df_emm_plot2 <- df_emm_plot %>%
  mutate(
    comp1 = str_trim(str_split_fixed(contrast, " - ", 2)[, 1]),
    comp2 = str_trim(str_split_fixed(contrast, " - ", 2)[, 2]),
    
    # use "vs" instead of "-"
    contrast_flipped = paste(comp2, "vs", comp1),
    
    estimate_log2_flipped = -estimate_log2,
    show_coef_flipped = round(estimate_log2_flipped, 2),
    fill_col_flipped = estimate_log2_flipped
  ) %>%
  select(-comp1, -comp2)


b <- ggplot(df_emm_plot2, aes(x = contrast_flipped, y = ASV, fill = fill_col_flipped)) +
  geom_tile(color = "grey90") +
  geom_text(aes(label = show_coef_flipped), size = 3) +
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

b
### example plots
#### Plot 
###boxplot
samdata_subset1 = metadata_table
list_NAMES
cor.test(genus_tab[["Prevotellaceae.UCG.001"]] ,genus_tab[["Bacteroides"]] )




samdata_subset1$Genus = genus_tab[["Marvinbryantia"]]
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

#b



# Plot
c <- data_to_plot %>%
  mutate(squirrel_year = paste(squirrel_id_letter, year_factor, sep = " | ")) %>%
  ggplot(aes(x = season_cat, y = Genus, group = squirrel_year, color = year_factor)) +
  
  # thin faint points & subject-level lines (colored by year)
  geom_point(size = 2, alpha = 0.2) +
  geom_line(linewidth = 0.8, alpha = 0.2) +
  
  # ribbon (prediction interval) — uses pred_trend data; facetting will match by mast
  geom_ribbon(data = pred_trend, 
              aes(x = season_cat, ymin = asymp.LCL, ymax = asymp.UCL, group = 1),
              fill = "gray30", alpha = 0.3, inherit.aes = FALSE)  +
  
  # two thicker trend lines: one for mast (orange), one for non-mast (green)
  geom_line(
    data = subset(pred_trend, mast == "mast"),
    aes(x = season_cat, y = response, group = 1),
    color = "#1b9e77",
    linewidth = 1.5,
    inherit.aes = FALSE
  ) +
  geom_line(
    data = subset(pred_trend, mast == "non-mast"),
    aes(x = season_cat, y = response, group = 1),
    color = "#7570b3",
    linewidth = 1.5,
    inherit.aes = FALSE
  ) +
  
  # facets, x-axis spacing, legend title, theme, and rotated x labels
  facet_wrap(~mast) +
  scale_x_discrete(expand = c(0.05, 0.05)) +
  labs(
    title = "Marvinbryantia",
    x = "Season",
    y = "Rarefied counts",
    color = "Year"          # legend title for the year mapping
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
c

a = a + theme(
  axis.text.x = element_blank(),axis.title.x = element_blank()
)
b = b + theme(
  axis.text.x = element_blank(),axis.title.x = element_blank()
)

ggarrange(a,b,c,heights = c(0.7, 0.7,1), ncol = 1, common.legend = TRUE)
