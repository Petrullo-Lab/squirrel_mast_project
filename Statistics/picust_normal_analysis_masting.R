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
library(tidyr)





source("~/Library/CloudStorage/OneDrive-UniversityofArizona/helpful_r_functions.R")

### read tabales in
setwd("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/Squirrels/analysis_data_oct_2025/16s/comb_ASV_tabs")
metadata_table = readRDS("new_meta_all_samples.RDS")
unique(metadata_table$sex)
picrust_tab = read.delim("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/Squirrels/analysis_data_oct_2025/16s/picrust/path_abun_unstrat.tsv")
rownames(picrust_tab) = picrust_tab$pathway
picrust_tab = picrust_tab[,-1]
picrust_tab = data.frame(t(picrust_tab))

metadata_table = readRDS("new_meta_all_samples_bdate.RDS")

### test
rownames(picrust_tab) 
metadata_table$name
rownames(picrust_tab) = gsub("\\.", "-", rownames(picrust_tab))
unique(rownames(picrust_tab)  == metadata_table$name) # doh

metadata_table <- metadata_table[match(rownames(picrust_tab), metadata_table$name), ]
unique(rownames(picrust_tab)  == metadata_table$name) # yay
sort(rowSums(picrust_tab))


metadata_table = metadata_table[metadata_table$day_of_year < 227,]
nrow(metadata_table)
### put year and mast
metadata_table$year_factor = as.factor(metadata_table$year)
metadata_table$mast = ifelse(metadata_table$year %in% c(2010, 2014, 2019, 2022), "mast", "non-mast")
table(metadata_table$gr)
metadata_table = metadata_table[metadata_table$gr != "BT",]
picrust_tab = picrust_tab[rownames(picrust_tab) %in% metadata_table$name,]
metadata_table$cones = ifelse(metadata_table$day_of_year > 200, "cone_available", "cone unavailable")
metadata_table$sept = ifelse(metadata_table$day_of_year > 200, "post", "pre")
table(metadata_table$sept, metadata_table$mast)

metadata_table$squirrel_id_letter = paste("a_",metadata_table$squirrel_id, sep = "")
metadata_table$day_scaled <- scale(metadata_table$day_of_year, center = TRUE, scale = TRUE)


metadata_table$season_cat <- cut(
  metadata_table$day_of_year,
  breaks = c(min(metadata_table$day_of_year)-1, 120, 180, max(metadata_table$day_of_year)+1),
  labels = c("Early Spring", "Late Spring", "Summer"),
  right = FALSE # so 60 is included in Early Spring, 120 in Late Spring, etc.
)





#genus_tab = data.frame(resss$ASVtab)
genus_tab = picrust_tab
genus_tab_rel = genus_tab
list_NAMES = colnames(genus_tab_rel)

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
alpha <- 0.10

#write.csv(emmeans_df,"/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/Squirrels/analysis_data_oct_2025/16s/picrust/emmeans_df_picrust.csv")

df_emm_plot <- emmeans_df %>%
  #  filter(!is.na(season_cat)) %>%
  mutate(
    signif = qvalue < alpha,
    show_coef = sprintf("%.2f", estimate_log2),
    fill_col = ifelse(signif, estimate_log2, NA)
  )
df_emm_plot = df_emm_plot[which(df_emm_plot$qvalue < 0.10),]
df_emm_plot_mast = df_emm_plot
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
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

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
alpha <- 0.10

#write.csv(emmeans_df,"/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/Squirrels/analysis_data_oct_2025/16s/picrust/emmeans_df_picrust_seasons.csv")
#emmeans_df = read.csv("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/Squirrels/analysis_data_oct_2025/16s/picrust/emmeans_df_picrust_seasons.csv")

df_emm_plot <- emmeans_df %>%
  #  filter(!is.na(season_cat)) %>%
  mutate(
    signif = qvalue < alpha,
    show_coef = sprintf("%.2f", estimate_log2),
    fill_col = ifelse(signif, estimate_log2, NA)
  )
df_emm_plot = df_emm_plot[which(df_emm_plot$qvalue < 0.1),]
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
 
### ok way to many. Let's just plot 10 most abudant changes per season comparisons
### calculate most abudant (averaged per sample)
read_picrust_names = read.delim("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/Squirrels/analysis_data_oct_2025/16s/picrust/map_metacyc-pwy_name.txt",
                                header = FALSE)
read_picrust_cat = read.delim("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/Squirrels/analysis_data_oct_2025/16s/picrust/map_metacyc-pwy_lineage.tsv",
                                header = FALSE)

read_picrust_cat_sep <- read_picrust_cat %>%
  separate(
    V2,
    into = c("Level1", "Level2", "Level3", "Level4"),
    sep = "\\|",
    fill = "right"
  )

head(read_picrust_cat_sep)
head(df_emm_plot)
library(dplyr)

df_emm_plot_annotated <- df_emm_plot %>%
  mutate(
    V1 = gsub("\\.", "-", sub("^X", "", ASV))
  ) %>%
  left_join(read_picrust_cat_sep, by = "V1")


avg_ab = sort(colMeans(picrust_tab), decreasing = TRUE)
df_emm_plot
df_emm_plot_spring_summer = df_emm_plot_annotated[df_emm_plot_annotated$contrast == "Early Spring - Summer",]


df_div <- df_emm_plot_spring_summer %>%
  mutate(
    Level1 = ifelse(is.na(Level1), "Non-categorized", Level1),
    Direction = case_when(
      show_coef > 0 ~ "Increasing",
      show_coef < 0 ~ "Decreasing"
    )
  ) %>%
  filter(!is.na(Direction)) %>%
  group_by(Level1, Direction) %>%
  summarise(Count = n(), .groups = "drop") %>%
  tidyr::complete(
    Level1,
    Direction = c("Increasing", "Decreasing"),
    fill = list(Count = 0)
  ) %>%
  mutate(
    Count_signed = ifelse(Direction == "Decreasing", -Count, Count)
  )

# ---- FIX: compute Level1 order safely ----
level_order <- df_div %>%
  group_by(Level1) %>%
  summarise(Total = sum(abs(Count_signed)), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  pull(Level1)

df_div$Level1 <- factor(df_div$Level1, levels = level_order)

# ensure Label exists (hide zeros)
df_div <- df_div %>% mutate(Label = ifelse(abs(Count_signed) == 0, "", as.character(abs(Count_signed))))

# split for reliable labeling
df_pos <- df_div %>% filter(Count_signed > 0)
df_neg <- df_div %>% filter(Count_signed < 0)

# ---- Plot (replace your previous earl_vs_summ block) ----
earl_vs_summ <- ggplot(df_div, aes(x = Count_signed, y = Level1, fill = Direction)) +
  geom_col(width = 0.7) +
  
  # positive labels: always placed to the right of bar
  geom_text(
    data = df_pos,
    aes(label = Label),
    hjust = -0.05,    # move label slightly right of the positive bar; tweak if needed
    size = 3
  ) +
  
  # negative labels: always placed to the left of bar
  geom_text(
    data = df_neg,
    aes(label = Label),
    hjust = 1.05,     # move label slightly left of the negative bar; tweak if needed
    size = 3
  ) +
  
  geom_vline(xintercept = 0, linewidth = 0.6) +
  scale_fill_manual(values = c("Increasing" = "red", "Decreasing" = "blue")) +
  labs(
    title = "Early spring vs Summer",
    x = "Number of pathways",
    y = NULL,
    fill = "Direction"
  ) +
  theme_bw(base_size = 12) +
  scale_x_continuous(expand = expansion(mult = c(0.30, 0.30))) +  # give ample room for labels
  coord_cartesian(clip = "off")                                    # allow labels outside plot

earl_vs_summ


### other contrast ---------
df_emm_plot_late_spring_summer = df_emm_plot_annotated[df_emm_plot_annotated$contrast == "Late Spring - Summer",]

df_div <- df_emm_plot_late_spring_summer %>%
  mutate(
    Level1 = ifelse(is.na(Level1), "Non-categorized", Level1),
    Direction = case_when(
      show_coef > 0 ~ "Increasing",
      show_coef < 0 ~ "Decreasing"
    )
  ) %>%
  filter(!is.na(Direction)) %>%
  group_by(Level1, Direction) %>%
  summarise(Count = n(), .groups = "drop") %>%
  tidyr::complete(
    Level1,
    Direction = c("Increasing", "Decreasing"),
    fill = list(Count = 0)
  ) %>%
  mutate(
    Count_signed = ifelse(Direction == "Decreasing", -Count, Count)
  )

# ---- FIX: compute Level1 order safely ----
level_order <- df_div %>%
  group_by(Level1) %>%
  summarise(Total = sum(abs(Count_signed)), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  pull(Level1)

df_div$Level1 <- factor(df_div$Level1, levels = level_order)

# ---- Plot ----
lat_vs_summer = # ensure Label exists (hide zeros)
  df_div <- df_div %>% mutate(Label = ifelse(abs(Count_signed) == 0, "", as.character(abs(Count_signed))))

# split for reliable labeling
df_pos <- df_div %>% filter(Count_signed > 0)
df_neg <- df_div %>% filter(Count_signed < 0)

# ---- Plot (replace your previous earl_vs_summ block) ----
lat_vs_summer <- ggplot(df_div, aes(x = Count_signed, y = Level1, fill = Direction)) +
  geom_col(width = 0.7) +
  
  # positive labels: always placed to the right of bar
  geom_text(
    data = df_pos,
    aes(label = Label),
    hjust = -0.05,    # move label slightly right of the positive bar; tweak if needed
    size = 3
  ) +
  
  # negative labels: always placed to the left of bar
  geom_text(
    data = df_neg,
    aes(label = Label),
    hjust = 1.05,     # move label slightly left of the negative bar; tweak if needed
    size = 3
  ) +
  
  geom_vline(xintercept = 0, linewidth = 0.6) +
  scale_fill_manual(values = c("Increasing" = "red", "Decreasing" = "blue")) +
  labs(
    title = "Late spring vs Summer",
    x = "Number of pathways",
    y = NULL,
    fill = "Direction"
  ) +
  theme_bw(base_size = 12) +
  scale_x_continuous(expand = expansion(mult = c(0.30, 0.30))) +  # give ample room for labels
  coord_cartesian(clip = "off")                                    # allow labels outside plot

lat_vs_summer





### other contrast ---------
df_emm_plot_late_spring_early = df_emm_plot_annotated[df_emm_plot_annotated$contrast == "Early Spring - Late Spring",]
unique(df_emm_plot_annotated$contrast)


df_div <- df_emm_plot_late_spring_early %>%
  mutate(
    Level1 = ifelse(is.na(Level1), "Non-categorized", Level1),
    Direction = case_when(
      show_coef > 0 ~ "Increasing",
      show_coef < 0 ~ "Decreasing"
    )
  ) %>%
  filter(!is.na(Direction)) %>%
  group_by(Level1, Direction) %>%
  summarise(Count = n(), .groups = "drop") %>%
  tidyr::complete(
    Level1,
    Direction = c("Increasing", "Decreasing"),
    fill = list(Count = 0)
  ) %>%
  mutate(
    Count_signed = ifelse(Direction == "Decreasing", -Count, Count)
  )

# ---- FIX: compute Level1 order safely ----
level_order <- df_div %>%
  group_by(Level1) %>%
  summarise(Total = sum(abs(Count_signed)), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  pull(Level1)

df_div$Level1 <- factor(df_div$Level1, levels = level_order)

# ---- Plot ----
earl_vs_lat = # ensure Label exists (hide zeros)
  df_div <- df_div %>% mutate(Label = ifelse(abs(Count_signed) == 0, "", as.character(abs(Count_signed))))

# split for reliable labeling
df_pos <- df_div %>% filter(Count_signed > 0)
df_neg <- df_div %>% filter(Count_signed < 0)

# ---- Plot (replace your previous earl_vs_summ block) ----
earl_vs_lat <- ggplot(df_div, aes(x = Count_signed, y = Level1, fill = Direction)) +
  geom_col(width = 0.7) +
  
  # positive labels: always placed to the right of bar
  geom_text(
    data = df_pos,
    aes(label = Label),
    hjust = -0.05,    # move label slightly right of the positive bar; tweak if needed
    size = 3
  ) +
  
  # negative labels: always placed to the left of bar
  geom_text(
    data = df_neg,
    aes(label = Label),
    hjust = 1.05,     # move label slightly left of the negative bar; tweak if needed
    size = 3
  ) +
  
  geom_vline(xintercept = 0, linewidth = 0.6) +
  scale_fill_manual(values = c("Increasing" = "red", "Decreasing" = "blue")) +
  labs(
    title = "Early spring vs Late spring",
    x = "Number of pathways",
    y = NULL,
    fill = "Direction"
  ) +
  theme_bw(base_size = 12) +
  scale_x_continuous(expand = expansion(mult = c(0.30, 0.30))) +  # give ample room for labels
  coord_cartesian(clip = "off")                                    # allow labels outside plot

earl_vs_lat


ggarrange(earl_vs_lat,earl_vs_summ, lat_vs_summer, common.legend = TRUE, ncol = 1 , align = "v")
#write.csv(df_emm_plot_annotated, "/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/Squirrels/analysis_data_oct_2025/16s/picrust/df_emm_plot_annotated.csv")


### Let's add some examples: 
# SCFA pathways - all significant!
scfa_pathways <- c("PWY.6588", "P108.PWY", "PWY.5677")

# Plant degradation - 3/4 significant!
plant_pathways <- c("PWY.5517", "PWY.6486", "PWY0.321")

# Define pathways
scfa_pathways <- c("PWY.6588", "P108.PWY", "PWY.5677", "PWY0.42")
plant_pathways <- c("PWY.5517", "PWY.6486", "PWY0.321")
all_pathways <- c(scfa_pathways, plant_pathways)
df_emm_plot_annotated$ASV1 = gsub("\\.","-",df_emm_plot_annotated$ASV)
df_emm_plot_annotated$name = read_picrust_names$V2[match(df_emm_plot_annotated$ASV1,read_picrust_names$V1)]

df_emm_plot_annotated_to_print = df_emm_plot_annotated[df_emm_plot_annotated$ASV %in%all_pathways, ]
df_emm_plot_annotated_to_print$name_to_print = paste(df_emm_plot_annotated_to_print$ASV1, df_emm_plot_annotated_to_print$name)

ggplot(df_emm_plot_annotated_to_print, aes(as.numeric(show_coef), contrast)) + geom_point() + facet_wrap(~name_to_print)

df_emm_plot_annotated_to_print <- df_emm_plot_annotated_to_print %>%
  mutate(
    signif_label = case_when(
      qvalue < 0.001 ~ "***",
      qvalue < 0.01 ~ "**",
      qvalue < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    signif_shape = ifelse(qvalue < 0.05, 19, 1)  # filled vs open circle
  )
ggplot(df_emm_plot_annotated_to_print, 
       aes(x = -as.numeric(show_coef), 
           y = name_to_print,
           color = as.numeric(show_coef) > 0)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(aes(shape = signif_label, size = signif_label), alpha = 0.8) +
  facet_wrap(~ contrast, ncol = 1,
             labeller = as_labeller(c(
               "Early Spring - Late Spring" = "Late vs Early Spring",
               "Early Spring - Summer"  = "Summer vs Early Spring",
               "Late Spring - Summer"  = "Summer vs Late Spring"
             ))
  )+
  scale_color_manual(values = c("FALSE" = "#E41A1C", "TRUE" = "#377EB8"),
                     labels = c("Enriched", "Depleted"),
                     name = "Direction") +
  scale_shape_manual(values = c("***" = 19, "**" = 19, "*" = 19, "ns" = 1),
                     name = "Significance") +
  scale_size_manual(values = c("***" = 4, "**" = 3.5, "*" = 3, "ns" = 2),
                    name = "Significance") +
  labs(
    x = "Log2 Fold Change",
    y = ""
  ) +
  theme_classic() +
  theme(
    legend.position = "top",
    legend.justification = "right",
    legend.box = "horizontal",
    legend.box.spacing = unit(0.2, "cm"),
    legend.spacing.x = unit(0.1, "cm"),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10)
  ) +
  guides(
    shape = guide_legend(nrow = 2),
    size = guide_legend(nrow = 2),
    color = guide_legend(nrow = 2)
  )






library(forcats)

# Combine all contrasts into one summarized table ----
df_all <- df_emm_plot_annotated %>%
  # recode NA level and direction
  mutate(
    Level1 = ifelse(is.na(Level1), "Non-categorized", Level1),
    Direction = case_when(
      show_coef > 0  ~ "Decreasing",
      show_coef < 0  ~ "Increasing",
      TRUE           ~ NA_character_    # drop no-change
    )
  ) %>%
  filter(!is.na(Direction)) %>%
  # count per contrast x Level1 x Direction
  group_by(contrast, Level1, Direction) %>%
  summarise(Count = n(), .groups = "drop") %>%
  # ensure zero rows exist (so stacking is stable)
  complete(contrast, Level1, Direction = c("Increasing", "Decreasing"), fill = list(Count = 0)) 

# Optionally: order Level1 by overall total across all contrasts (largest first)
level_order <- df_all %>%
  group_by(Level1) %>%
  summarise(Total = sum(Count), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  pull(Level1)

df_all <- df_all %>% mutate(Level1 = factor(Level1, levels = level_order))
df_all$Direction = factor(df_all$Direction, levels = c("Increasing", "Decreasing"))
# Add label column (hide zero labels)
df_all <- df_all %>% mutate(Label = ifelse(Count == 0, "", as.character(Count)))
unique(df_all$contrast)
p <- ggplot(df_all, aes(x = Direction, y = Count, fill = Level1)) +
  geom_col(width = 0.75, position = "stack") +
  facet_wrap(~ contrast, ncol = 3,
             labeller = as_labeller(c(
               "Early Spring - Late Spring" = "Late vs Early Spring",
               "Early Spring - Summer"  = "Summer vs Early Spring",
               "Late Spring - Summer"  = "Summer vs Late Spring"
             ))
  ) +
  scale_fill_brewer(palette = "Set3", name = "MetaCyc category") +
  labs(
    y = "Number of pathways"
  ) +
  theme_classic() 
  
p

df_emm_plot

df_summary <- df_all %>%
  group_by(contrast, Direction) %>%
  summarise(Total = sum(Count), .groups = "drop")

df_summary
read_picrust_names = read.delim("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/Squirrels/analysis_data_oct_2025/16s/picrust/map_metacyc-pwy_name.txt",
                                header = FALSE)
read_picrust_cat = read.delim("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/Squirrels/analysis_data_oct_2025/16s/picrust/map_metacyc-pwy_lineage.tsv",
                              header = FALSE)

read_picrust_cat_sep <- read_picrust_cat %>%
  separate(
    V2,
    into = c("Level1", "Level2", "Level3", "Level4"),
    sep = "\\|",
    fill = "right"
  )

head(read_picrust_cat_sep)
head(df_emm_plot)
library(dplyr)

df_emm_plot_annotated_mast <- df_emm_plot_mast %>%
  mutate(
    V1 = gsub("\\.", "-", sub("^X", "", ASV))
  ) %>%
  left_join(read_picrust_cat_sep, by = "V1")



df_emm_plot_annotated_mast$ASV1 = gsub("\\.","-",df_emm_plot_annotated_mast$ASV)
df_emm_plot_annotated_mast$name = read_picrust_names$V2[match(df_emm_plot_annotated_mast$ASV1,read_picrust_names$V1)]

df_emm_plot_annotated_to_print = df_emm_plot_annotated_mast
df_emm_plot_annotated_to_print$name_to_print = paste(df_emm_plot_annotated_to_print$ASV1, df_emm_plot_annotated_to_print$name)

ggplot(df_emm_plot_annotated_to_print, aes(as.numeric(show_coef), contrast)) + geom_point() + facet_wrap(~name_to_print)

df_emm_plot_annotated_to_print$contrast = "Mast - Non Mast (Summer)"
df_emm_plot_annotated_to_print <- df_emm_plot_annotated_to_print %>%
  mutate(
    signif_label = case_when(
      qvalue < 0.001 ~ "***",
      qvalue < 0.01 ~ "**",
      qvalue < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    signif_shape = ifelse(qvalue < 0.05, 19, 1)  # filled vs open circle
  )
df_emm_plot_annotated_to_print$name_to_print


df_emm_plot_annotated_to_print <- df_emm_plot_annotated_to_print %>%
  mutate(name_to_print = case_when(
    ASV1 == "PWY-7090" ~ "PWY-7090 LPS sugar biosynthesis",
    TRUE ~ name_to_print
  ))
ggplot(df_emm_plot_annotated_to_print, 
       aes(y = as.numeric(show_coef), 
           x = name_to_print,
           color = as.numeric(show_coef) > 0)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(aes(shape = signif_label, size = signif_label), alpha = 0.8) +
  facet_wrap(~contrast, ncol = 1) +
  scale_color_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "#377EB8"),
                     labels = c("Depleted", "Enriched"),
                     name = "Direction") +
  scale_shape_manual(values = c("***" = 19, "**" = 19, "*" = 19, "ns" = 1),
                     name = "Significance") +
  scale_size_manual(values = c("***" = 4, "**" = 3.5, "*" = 3, "ns" = 2),
                    name = "Significance") +
  labs(
    y = "Log2 Fold Change",
    x = ""
  ) +
  theme_classic() +
  theme(
    legend.position = "top",
    legend.justification = "right",
    legend.box = "horizontal",
    legend.box.spacing = unit(0.2, "cm"),
    legend.spacing.x = unit(0.1, "cm"),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  guides(
    shape = guide_legend(nrow = 2),
    size = guide_legend(nrow = 2),
    color = guide_legend(nrow = 2)
  )
