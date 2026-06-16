#### get food data #####
library(dplyr)
library(ggplot2)
library(lubridate)
library(scales)
library(krsp)
library(tidyr)
library(tibble)

setwd("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/Squirrels/analysis_data_oct_2025/16s/comb_ASV_tabs")
source("~/Library/CloudStorage/OneDrive-UniversityofArizona/helpful_r_functions.R")


### behaviour
con <- krsp_connect (host = "krsp.cepb5cjvqban.us-east-2.rds.amazonaws.com",
                     dbname ="krsp",
                     username = Sys.getenv("krsp_user"),
                     password = Sys.getenv("krsp_password")
)
con
# Pull the trapping table
behaviour = tbl(con, "behaviour") %>%
  collect()


behaviour = behaviour[behaviour$grid %in% c("SU" ,"KL" ),]
feed = behaviour[behaviour$behaviour == 1,]
feed$detail = as.factor(feed$detail)

sort(table(feed$detail), decreasing = TRUE)

feed$category <- dplyr::case_when(
  feed$detail %in% c("1") ~ "animal material",
  
  # regular cones
  feed$detail %in% c("2") ~ "cone",
  
  # new cones + subcodes
  feed$detail %in% c("6", "6A", "6B", "6C", "6D",
                     "6E", "6F", "6G", "6H") ~ "new cone",
  
  # buds
  feed$detail %in% c("3") ~ "buds",
  
  feed$detail %in% c("13") ~ "needle",
  feed$detail %in% c("19") ~ "witches' broom",
  # mushrooms / truffles
  feed$detail %in% c("4", "04", "31", "32") ~ "mushroom_truffle",
  
  TRUE ~ "other"
)


table(feed$category)





feed2 <- feed %>%
  mutate(
    date = ymd(date, quiet = TRUE),
    year = year(date),
    week = isoweek(date)
  ) %>%
  filter(!is.na(date), !is.na(category))

feed_subset = feed2[feed2$year %in% c(2008:2023),]
feed_subset <- feed_subset[
  feed_subset$year != 2013 &
    feed_subset$year != 2020 &
    feed_subset$year != 2021,
]


weekly_rel <- feed_subset %>%
  count(year, week, category) %>%
  group_by(year, week) %>%
  mutate(
    rel_abundance = n / sum(n),
    week_date = as.Date(paste0(year, "-01-01")) + (week - 1) * 7,
    doy = yday(week_date),
    year_week = paste0(year, "_W", sprintf("%02d", week))
  ) %>%
  ungroup()

ggplot(weekly_rel, aes(x = doy, y = rel_abundance, fill = category)) +
  geom_col(width = 6) +
  facet_wrap(~ year) +
  geom_vline(xintercept = c(120, 180), linetype = "dashed") +
  scale_y_continuous(labels = percent_format()) +
  scale_x_continuous(breaks = seq(1, 365, by = 30)) +
  labs(x = "Day of year", y = "% total observation / week", fill = "Category") +
  theme_bw() +
  scale_fill_brewer(palette = "Set1")



weekly_counts <- feed_subset %>%
  count(year, week)

overall_summary <- weekly_counts %>%
  summarise(
    mean_n = mean(n, na.rm = TRUE),
    sd_n   = sd(n, na.rm = TRUE)
  )

overall_summary





### mantel test with ASV composition -----------
metadata_table = readRDS("metadata_table_sub.RDS")
ASVtable = readRDS("ASVtable_sub.RDS")
rownames(ASVtable)
head(metadata_table)


metadata_table <- metadata_table %>%
  mutate(
    week = isoweek(date),
    year_week = paste0(year(date), "_W", sprintf("%02d", week))
  )


#### Mantel test with genus composition
library(dplyr)
library(tibble)

# change `name1` to the sample-ID column in your metadata if needed
ASVtable = ASVtable[rownames(ASVtable) %in% metadata_table$name,]
unique(rownames(ASVtable) %in% metadata_table$name)
# change `year_week` to the name of the grouping column you created
asv_by_week <- ASVtable %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  left_join(
    metadata_table %>%
      mutate(sample_id = name) %>%
      select(sample_id, year_week),
    by = "sample_id"
  ) %>%
  group_by(year_week) %>%
  summarise(across(where(is.numeric), sum), .groups = "drop") %>%
  column_to_rownames("year_week") %>%
  as.matrix()



weekly_rel_mat <- weekly_rel %>%
  select(year_week, category, rel_abundance) %>%
  pivot_wider(
    names_from = category,
    values_from = rel_abundance,
    values_fill = 0
  ) %>%
  column_to_rownames("year_week") %>%
  as.matrix()


weekly_rel_mat = weekly_rel_mat[rownames(weekly_rel_mat) %in% rownames(asv_by_week),]
asv_by_week = asv_by_week[rownames(asv_by_week) %in% rownames(weekly_rel_mat) ,]
rownames(asv_by_week) == rownames(weekly_rel_mat) #ok

library(vegan)


asv_by_week_rel <- sweep(
  asv_by_week,
  1,
  rowSums(asv_by_week),
  FUN = "/"
)


weekly_rel_mat_rel <- sweep(
  weekly_rel_mat,
  1,
  rowSums(weekly_rel_mat),
  FUN = "/"
)


# Bray-Curtis distance matrices
asv_dist  <- vegdist(asv_by_week_rel, method = "bray")
feed_dist <- vegdist(weekly_rel_mat_rel, method = "bray")
dim(asv_dist)
dim(feed_dist)
# Mantel test
mantel_result <- mantel(
  asv_dist,
  feed_dist,
  method = "spearman",   # or "pearson"
  permutations = 9999
)

mantel_result$statistic
mantel_result$signif

# convert to vectors
dd <- data.frame(
  microbiome_dist = as.vector(asv_dist),
  diet_dist       = as.vector(feed_dist)
)

# plot
ggplot(dd, aes(x = diet_dist, y = microbiome_dist)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_classic() +
  labs(
    x = "Diet community dissimilarity (Bray-Curtis)",
    y = "Microbiome dissimilarity (Bray-Curtis)",
    title = "Microbiome~Diet relationship"
  )









set.seed(12)
asv_by_week_int = asv_by_week
asv_by_week_rar = rrarefy(asv_by_week_int, min(rowSums(asv_by_week_int)))
richness = specnumber(asv_by_week_rar)
library(ggplot2)
library(ggpubr)

df <- data.frame(
  richness = richness,
  old_cones = weekly_rel_mat_rel[,1],
  new_cones = weekly_rel_mat_rel[,6]
)

# OLD CONES
p1 <- ggplot(df, aes(x = old_cones, y = richness)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  stat_cor(method = "pearson") +
  theme_classic() +
  labs(
    title = "Old cones",
    x = "Fraction of old cone feeding events",
    y = "ASV Richness"
  )

# NEW CONES
p2 <- ggplot(df, aes(x = new_cones, y = richness)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  stat_cor(method = "pearson") +
  theme_classic() +
  labs(
    title = "New cones",
    x = "Fraction of new cone feeding events",
    y = "ASV Richness"
  )

p1
p2
ggarrange(p1, p2)










### mantel test with genus tables
Taxtable = readRDS("Taxtab_total_silva.RDS")
genus = generate.tax.summary.modified(asv_by_week,as.data.frame(Taxtable))
genus = genus$tax7

genus <- genus %>%
  as.data.frame() %>%
  na.omit() %>%                     # remove rows containing NA
  as.matrix()


asv_dist  <- vegdist(t(genus), method = "bray")
dim(asv_dist)
dim(feed_dist)
# Mantel test
mantel_result <- mantel(
  asv_dist,
  feed_dist,
  method = "spearman",   # or "pearson"
  permutations = 9999
)

mantel_result

# convert to vectors
dd <- data.frame(
  microbiome_dist = as.vector(asv_dist),
  diet_dist       = as.vector(feed_dist)
)

# plot
ggplot(dd, aes(x = diet_dist, y = microbiome_dist)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_classic() +
  labs(
    x = "Diet community dissimilarity (Bray-Curtis)",
    y = "Microbiome dissimilarity (Bray-Curtis)",
    title = "Microbiome~Diet relationship at genus level"
  )



### Picrust results ---------
picrust_tab = readRDS("picrust_tab.RDS")
metadata_table = readRDS("metadata_table_sub.RDS")



metadata_table <- metadata_table %>%
  mutate(
    week = isoweek(date),
    year_week = paste0(year(date), "_W", sprintf("%02d", week))
  )


# change `name1` to the sample-ID column in your metadata if needed
picrust_tab = picrust_tab[rownames(picrust_tab) %in% metadata_table$name,]
unique(rownames(ASVtable) %in% metadata_table$name)
# change `year_week` to the name of the grouping column you created
asv_by_week <- picrust_tab %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  left_join(
    metadata_table %>%
      mutate(sample_id = name) %>%
      select(sample_id, year_week),
    by = "sample_id"
  ) %>%
  group_by(year_week) %>%
  summarise(across(where(is.numeric), sum), .groups = "drop") %>%
  column_to_rownames("year_week") %>%
  as.matrix()



weekly_rel_mat <- weekly_rel %>%
  select(year_week, category, rel_abundance) %>%
  pivot_wider(
    names_from = category,
    values_from = rel_abundance,
    values_fill = 0
  ) %>%
  column_to_rownames("year_week") %>%
  as.matrix()


weekly_rel_mat = weekly_rel_mat[rownames(weekly_rel_mat) %in% rownames(asv_by_week),]
asv_by_week = asv_by_week[rownames(asv_by_week) %in% rownames(weekly_rel_mat) ,]
rownames(asv_by_week) == rownames(weekly_rel_mat) #ok

library(vegan)


asv_by_week_rel <- sweep(
  asv_by_week,
  1,
  rowSums(asv_by_week),
  FUN = "/"
)


weekly_rel_mat_rel <- sweep(
  weekly_rel_mat,
  1,
  rowSums(weekly_rel_mat),
  FUN = "/"
)


# Bray-Curtis distance matrices
asv_dist  <- vegdist(asv_by_week_rel, method = "bray")
feed_dist <- vegdist(weekly_rel_mat, method = "bray")
dim(asv_dist)
dim(feed_dist)
# Mantel test
mantel_result <- mantel(
  asv_dist,
  feed_dist,
  method = "spearman",   # or "pearson"
  permutations = 9999
)

mantel_result$statistic
mantel_result$signif


### correlation with funtional richness:
set.seed(12)
asv_by_week_int = round(asv_by_week)
asv_by_week_rar = rrarefy(asv_by_week_int, min(rowSums(asv_by_week_int)))
richness = specnumber(asv_by_week_rar)
library(ggplot2)
library(ggpubr)

df <- data.frame(
  richness = richness,
  old_cones = weekly_rel_mat_rel[,1],
  new_cones = weekly_rel_mat_rel[,6]
)

# OLD CONES
p1 <- ggplot(df, aes(x = old_cones, y = richness)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  stat_cor(method = "pearson") +
  theme_classic() +
  labs(
    title = "Old cones",
    x = "Fraction of old cone feeding events",
    y = "Pathway Richness"
  )

# NEW CONES
p2 <- ggplot(df, aes(x = new_cones, y = richness)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  stat_cor(method = "pearson") +
  theme_classic() +
  labs(
    title = "New cones",
    x = "Fraction of new cone feeding events",
    y = "Pathway Richness"
  )

p1
p2
ggarrange(p1, p2)
### loop within a loop
results_loop = data.frame()
for(i in 1:ncol(weekly_rel_mat_rel))
{
  for(j in 1:ncol(asv_by_week_rel))
  {
    temp_pi = asv_by_week_rel[,j]
    temp_food = weekly_rel_mat_rel[,i]
    test = cor.test(temp_pi,temp_food, method = "spearman")
    res_temp = c(colnames(asv_by_week_rel)[j],colnames(weekly_rel_mat_rel)[i] , test$statistic, test$parameter, test$p.value, test$estimate)
    results_loop = rbind(res_temp,results_loop)
}
}

colnames(results_loop) = c("patwhay", "food", "stat", "p", "rho")
results_loop$p_adj = p.adjust(results_loop$p, method = "fdr")
results_loop_sig = results_loop[results_loop$p_adj < 0.05,]

unique(results_loop_sig$food )
### new cone 
results_loop_sig_new_cone = results_loop_sig[results_loop_sig$food == "new cone",]
results_loop_sig_new_cone = results_loop_sig_new_cone[results_loop_sig_new_cone$rho > 0,]



### old cone 
results_loop_sig_old_cone = results_loop_sig[results_loop_sig$food == "cone",]
results_loop_sig_old_cone = results_loop_sig_old_cone[results_loop_sig_old_cone$rho > 0,]
results_loop_sig_old_cone

### new cone 
results_loop_sig_new_cone = results_loop_sig[results_loop_sig$food == "new cone",]
results_loop_sig_new_cone = results_loop_sig_new_cone[results_loop_sig_new_cone$rho > 0,]



### old cone 
results_loop_sig_buds = results_loop_sig[results_loop_sig$food == "buds",]
results_loop_sig_buds = results_loop_sig_buds[results_loop_sig_buds$rho > 0,]
results_loop_sig_buds

### old cone 
results_loop_sig_broom = results_loop_sig[results_loop_sig$food == "witches' broom",]
results_loop_sig_broom = results_loop_sig_broom[results_loop_sig_broom$rho > 0,]
results_loop_sig_broom

results_loop_sig_mushroom = results_loop_sig[results_loop_sig$food == "mushroom_truffle",]
results_loop_sig_mushroom = results_loop_sig_mushroom[results_loop_sig_mushroom$rho > 0,]
results_loop_sig_mushroom



library(forcats)
library(tidytext)   # for reorder_within()

# top 5 positive correlations per food category
foods_keep = c("mushroom_truffle", "witches' broom", "buds", "new cone","cone")
top5_paths <- results_loop_sig %>%
  filter(food %in% foods_keep,
         rho > 0)

top5_paths

top5_paths$rho <- as.numeric(top5_paths$rho)

read_picrust_names = read.delim("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/Squirrels/analysis_data_oct_2025/16s/picrust/map_metacyc-pwy_name.txt",
                                header = FALSE)

read_picrust_names
top5_paths_annotated <- top5_paths %>%
  mutate(
    V1 = gsub("\\.", "-", sub("^X", "", patwhay))
  ) %>%
  left_join(read_picrust_names, by = "V1") 


table(top5_paths_annotated$food)
top5_paths_annotated = top5_paths_annotated %>%
  mutate(
    pathway_label = ifelse(is.na(V2), patwhay, V2),
    pathway_plot = reorder_within(pathway_label, rho, food)
  )

#write.csv(top5_paths_annotated, "cor_pathway_food.csv")
getwd()


emmeans_df = read.csv("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/Squirrels/analysis_data_oct_2025/16s/picrust/emmeans_df_picrust_seasons.csv")
emmeans_df_pwd = emmeans_df[emmeans_df$ASV %in% top5_paths$patwhay,]
top5_paths_pwd = top5_paths[top5_paths$patwhay %in% emmeans_df$ASV,]
length(unique(top5_paths_pwd$patwhay))
length(unique(emmeans_df$ASV))









top5_paths = top5_paths %>% group_by(food) %>%
  slice_max(order_by = rho, n = 10) %>%
  arrange(food, desc(rho)) %>%
  ungroup()

top5_paths$rho <- as.numeric(top5_paths$rho)

# reorder within each facet
top5_paths <- top5_paths %>%
  mutate(pathway_plot = reorder_within(patwhay, rho, food))

# plot
ggplot(top5_paths,
       aes(x = pathway_plot,
           y = rho,
           fill = food)) +
  
  geom_col(width = 0.8) +
  
  coord_flip() +
  
  scale_x_reordered() +
  
  facet_wrap(~food,
             scales = "free_y") +
  
  labs(
    x = "Pathway",
    y = "Spearman rho",
    title = "Top 5 Positive Pathway Correlations by Food Type"
  ) +
  
  theme_bw(base_size = 12) +
  
  theme(
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold"),
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )



#read_picrust_cat = read.delim("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/Squirrels/analysis_data_oct_2025/16s/picrust/map_metacyc-pwy_lineage.tsv",
#                              header = FALSE)

#read_picrust_cat_sep <- read_picrust_cat %>%
#  separate(
#    V2,
#    into = c("Level1", "Level2", "Level3", "Level4"),
##    sep = "\\|",
#    fill = "right"
#  )

#head(read_picrust_cat_sep)
top5_paths_annotated <- top5_paths %>%
  mutate(
    V1 = gsub("\\.", "-", sub("^X", "", patwhay))
  ) %>%
  left_join(read_picrust_names, by = "V1") %>%
  mutate(
    pathway_label = ifelse(is.na(V2), patwhay, V2),
    pathway_plot = reorder_within(pathway_label, rho, food)
  )

ggplot(top5_paths_annotated,
       aes(x = pathway_plot,
           y = rho)) +
  geom_col(width = 0.8) +
  coord_flip() +
  scale_x_reordered() +
  facet_wrap(~food, scales = "free_y", ncol = 1, nrow = 5) +
  labs(
    x = "Pathway",
    y = "Spearman rho",
    title = "Top 10 Positive Pathway Correlations by Food Type"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold"),
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )



