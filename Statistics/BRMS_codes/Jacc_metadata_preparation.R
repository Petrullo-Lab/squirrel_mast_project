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
# masting_factor



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
# List to store cleaned, melted Jacc matrices per genus
Jacc_melted_by_genus <- list()
for (genus in names(ASVtabs_by_genus)) {
  asvtab <- ASVtabs_by_genus[[genus]]
  # Skip genera with fewer than 2 ASVs
  if (nrow(asvtab) < 2) next
  # Transpose to samples x ASVs
  asvtab_t <- t(asvtab)
  # Convert to binary presence/absence
  asvtab_bin <- asvtab_t > 0
  # Compute binary Jacc distance (Sørensen)
  jacc_dist <- vegdist(asvtab_bin, method = "jaccard", binary = TRUE)
  jacc_mat <- as.matrix(jacc_dist)
  jacc_mat = 1 - jacc_mat
  # Melt the distance matrix
  jacc_melt <- melt(jacc_mat, varnames = c("Sample1", "Sample2"), value.name = "Jacc_Similarity")
  # Remove self-comparisons and duplicate pairs
  jacc_melt_filtered <- jacc_melt %>%
    filter(Sample1 != Sample2) %>%
    mutate(pair_min = pmin(as.character(Sample1), as.character(Sample2)),
           pair_max =  pmax(as.character(Sample1), as.character(Sample2))) %>%
    distinct(pair_min, pair_max, .keep_all = TRUE) %>%
    select(Sample1, Sample2, Jacc_Similarity)
  # Add genus column
  jacc_melt_filtered$Genus <- genus
  # Save to list
  Jacc_melted_by_genus[[genus]] <- jacc_melt_filtered
}








#### STEP 4 match the metadata in the columns using the one I already calculated for my Previous table: 
colnames(full_match_meta)
jacc_mat_melt_filtered <- full_match_meta %>% ### It should already good but check if triangle is taken
  mutate(
    pair_min = pmin(as.character(Sample1), as.character(Sample2)),
    pair_max = pmax(as.character(Sample1), as.character(Sample2)),
    pair_id = paste(pair_min, pair_max, sep = "_")
  )
Jacc_aligned_by_genus <- list()

for (genus in names(Jacc_melted_by_genus)) {
  
  genus_df <- Jacc_melted_by_genus[[genus]]
  
  # Ensure Sample1 and Sample2 are character
  genus_df <- genus_df %>%
    mutate(
      pair_min = pmin(as.character(Sample1), as.character(Sample2)),
      pair_max = pmax(as.character(Sample1), as.character(Sample2)),
      pair_id = paste(pair_min, pair_max, sep = "_")
    )
  
  # Join with global table to preserve order and alignment
  aligned <- jacc_mat_melt_filtered %>%
    select(Sample1, Sample2, pair_id) %>%
    left_join(genus_df %>% select(pair_id, Jacc_Similarity), by = "pair_id")# %>%
  #rename(Sample1 = Var1, Sample2 = Var2)
  
  # Add genus column
  aligned$Genus <- genus
  
  # Save aligned table
  Jacc_aligned_by_genus[[genus]] <- aligned
}

as_df = Jacc_aligned_by_genus$Prevotella
unique(as_df$Sample1 == jacc_mat_melt_filtered$Sample1) # Ok 
unique(as_df$Sample2 == jacc_mat_melt_filtered$Sample2) # Ok
hist(as_df$Jacc_Similarity)


### I would add how many ASVs are found in each comparison as in case, we remove the comparisons where 0 asvs where found
for (genus in names(Jacc_aligned_by_genus)) {
  
  # Get the aligned melted distance matrix
  df <- Jacc_aligned_by_genus[[genus]]
  
  # Get the ASV table for this genus (ASVs x Samples)
  asvtab <- ASVtabs_by_genus[[genus]]
  
  # Transpose to Samples x ASVs
  asvtab_t <- t(asvtab)
  
  # Calculate number of ASVs present (count of non-zero ASVs) per sample
  asv_counts <- rowSums(asvtab_t > 0)
  
  # Add columns to the melted distance table
  df$Sample1_ASVcount <- asv_counts[as.character(df$Sample1)]
  df$Sample2_ASVcount <- asv_counts[as.character(df$Sample2)]
  
  # Replace the entry in the list
  Jacc_aligned_by_genus[[genus]] <- df
}


dfff =Jacc_aligned_by_genus$Prevotella

### Append metadata to it: 
all(
  Jacc_aligned_by_genus[[1]]$Sample1 == jacc_mat_melt_filtered$Sample1 &
    Jacc_aligned_by_genus[[1]]$Sample2 == jacc_mat_melt_filtered$Sample2
) # good 

# Columns you want to add
#mother_offspring*masting+mother_offspring*same_season+
#  same_year+dist_m+same_sex+read_dist+run+age_dist+
#  (1|mm(Sample1,Sample2)) + (1|mm(squirrel1,squirrel2)),



meta_cols <- c("squirrel1", "squirrel2", "Sample1", "Sample2", "age_dist", "run",
               "read_dist", "year_dist" ,"same_sex", "dist_m",  "same_year", "same_season","mother_offspring",
               "masting", "masting_factor")

for (genus in names(Jacc_aligned_by_genus)) {
  df <- Jacc_aligned_by_genus[[genus]]
  
  # Add metadata from the same row in jacc_mat_melt_filtered
  df[, meta_cols] <- jacc_mat_melt_filtered[, meta_cols]
  
  # Save back
  Jacc_aligned_by_genus[[genus]] <- df
}


df_check = Jacc_aligned_by_genus$Ruminococcus




#### Remove zeros from table
for (genus in names(Jacc_aligned_by_genus)) {
  df <- Jacc_aligned_by_genus[[genus]]
  
  # Remove rows where both samples have zero ASVs
  df <- df %>%
    filter(!(Sample1_ASVcount == 0 & Sample2_ASVcount == 0))
  
  # Save back
  Jacc_aligned_by_genus[[genus]] <- df
}

saveRDS(Jacc_aligned_by_genus,"/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/Squirrels/masting_project/BRMS_data_again/Jacc_aligned_by_genus_filtered_new.RDS")
jaccccccc = readRDS("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/Squirrels/masting_project/BRMS_data_again/Jacc_aligned_by_genus_filtered_new.RDS")
bccccc = readRDS("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/Squirrels/masting_project/BRMS_data_again/BC_aligned_by_genus_filtered_new.RDS")



nnnnnn = names(Jacc_aligned_by_genus)
#saveRDS(nnnnnn, "/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/Squirrels/masting_project/BRMS_data_again/genera_names_new.rds")
#Jacc_aligned_by_genus = readRDS("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/Squirrels/masting_project/BRMS_data_again/Jacc_aligned_by_genus_filtered_new.RDS")
mv = bccccc$Marvinbryantia
mv$Jacc_Similarity