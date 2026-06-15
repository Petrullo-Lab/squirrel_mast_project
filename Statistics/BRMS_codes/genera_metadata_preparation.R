### Load TAXtab --------
library(vegan)
library(reshape2)
library(dplyr)



source("~/Library/CloudStorage/OneDrive-UniversityofArizona/helpful_r_functions.R")
setwd("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/Squirrels/analysis_data_oct_2025/16s/comb_ASV_tabs")


### Read tables in ----
ASVtable = readRDS("new_ASVtable_all_samples.RDS")
Taxtable = as.data.frame(readRDS("Taxtab_total_silva.RDS"))
metadata_table = readRDS("new_meta_all_samples.RDS")
#ASVtab_16s.rar = readRDS("ASVtab_16s.rar.RDS")
#ASVtab_16s.rar_t = data.frame(t(ASVtab_16s.rar))
full_match_meta = readRDS("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/Squirrels/masting_project/BRMS_data_again/Full_data_set_BRMS.RDS")
### Read the metadata_tab
ASVtable = rrarefy(ASVtable, 10000)
rowSums(ASVtable)

colnames(full_match_meta)
### Top 50 genera  ----
top_50_genera = readRDS("list_top50_genera.RDS")
### THis I have to read in!
genera_to_analyze = top_50_genera
head(Taxtable)
Taxtable$Genus[is.na(Taxtable$Genus)] = "unknown"
Taxtable$Genus = gsub(" ",".",Taxtable$Genus)
Taxtable$Genus = gsub("-",".",Taxtable$Genus)
Taxtable$Genus = gsub("\\[","X.",Taxtable$Genus)
Taxtable$Genus = gsub("\\]",".",Taxtable$Genus)
Taxtable$Genus = gsub("28","X28",Taxtable$Genus)
unique(Taxtable$Genus)
genera_to_analyze[which(genera_to_analyze %nin% Taxtable$Genus)]
#Taxtable$Genus
genus_rows <- lapply(top_50_genera, function(g) {
  rownames(Taxtable)[Taxtable$Genus == g]
})
names(genus_rows) <- top_50_genera


ASVtable_t = data.frame(t(ASVtable))
#### STEP 2 create a series of ASVtabs for each genus (maybe avoiding 
# 4. Create list of ASV tables by genus
ASVtabs_by_genus <- list()
#genus = genera_to_analyze[1]
ASVtabs_by_genus <- list()
for (genus in names(genus_rows)) {
  ids <- genus_rows[[genus]]                 # character vector of ASV IDs
  ASVtabs_by_genus[[genus]] <- ASVtable_t[ids, , drop = FALSE]
}
dff = ASVtabs_by_genus$Marvinbryantia


#### STEP 3 write a script that goes through the list 
#### calculates melted ASV tabs for each ASVtab and save it in a list.
full_match_meta <- full_match_meta %>%
  mutate(
    pair_min = pmin(as.character(Sample1), as.character(Sample2)),
    pair_max = pmax(as.character(Sample1), as.character(Sample2)),
    pair_id = paste(pair_min, pair_max, sep = "_")
  )


jaccard_aligned_by_genus <- list()

meta_cols <- c("squirrel1", "squirrel2", "Sample1", "Sample2",
               "age_dist", "run", "read_dist", "year_dist",
               "same_sex", "dist_m", "same_year", "same_season",
               "mother_offspring", "masting_factor")

for (genus in names(jaccard_melted_by_genus)) {
  
  genus_df <- jaccard_melted_by_genus[[genus]] %>%
    mutate(
      pair_min = pmin(as.character(Sample1), as.character(Sample2)),
      pair_max = pmax(as.character(Sample1), as.character(Sample2)),
      pair_id = paste(pair_min, pair_max, sep = "_")
    ) %>%
    select(pair_id, Jaccard_Similarity)
  
  aligned <- full_match_meta %>%
    select(all_of(meta_cols), pair_id) %>%
    inner_join(genus_df, by = "pair_id") %>%
    mutate(Genus = genus)
  
  jaccard_aligned_by_genus[[genus]] <- aligned
}
check = jaccard_aligned_by_genus[[25]]

summary(check$Jaccard_Similarity)


check[789,]
head(full_match_meta)
full_match_meta[full_match_meta$Sample1 == "Sample297_S303_L001"&  full_match_meta$Sample2 =="Sample100_S102_L001",]
saveRDS(jaccard_aligned_by_genus,"/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/Squirrels/masting_project/BRMS_data_again/jaccard_aligned_by_genus_filtered_new.RDS")
nnnnnn = names(jaccard_aligned_by_genus)
saveRDS(nnnnnn, "/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/Squirrels/masting_project/BRMS_data_again/genera_names_new.rds")


