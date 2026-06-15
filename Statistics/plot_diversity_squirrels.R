### mother offspring pairs
library(dplyr)


### Pedigree and litter download -----------------
#### Use the fil that lauren sent me on Jun 3rd
### mother offpspring pair
setwd("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/Squirrels/analysis_data_oct_2025/16s/comb_ASV_tabs")
metadata_table = readRDS("metadata_mast.RDS")
pedigree_db = read.csv("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/Squirrels/masting_project/pup_summary_new.csv")
metadata_table$dam = pedigree_db$mother_id[match(metadata_table$squirrel_id, pedigree_db$squirrel_id)]
sq = unique(metadata_table$squirrel_id[which(metadata_table$dam %in% metadata_table$squirrel_id)]) # 25 ok
metadata_table$litter = pedigree_db$litter_id[match(metadata_table$squirrel_id, pedigree_db$squirrel_id)]
metadata_table_NA = metadata_table[is.na(metadata_table$dam),]
unique(metadata_table_NA$squirrel_id)

### metadata day of year
#metadata_table = metadata_table[metadata_table$day_of_year < 227,]
hist(metadata_table$day_of_year)

#### Create metadata -------------
#### Read in the jacc_dissimilarity matrix
jacc_dissimilarity = readRDS("bac_dist_BC_rar.RDS")
BC_dissimilarity = readRDS("bac_dist_BC_rar.RDS")
melted_jaccard <- reshape2::melt(as.matrix(jacc_dissimilarity), varnames = c("Sample1", "Sample2"), value.name = "Jaccard_Dissimilarity")
melted_BC <- reshape2::melt(as.matrix(BC_dissimilarity), varnames = c("Sample1", "Sample2"), value.name = "BC_Dissimilarity")
unique(melted_BC$Sample1 == melted_jaccard$Sample1)
unique(melted_BC$Sample2 == melted_jaccard$Sample2) # Ok

#### Reshape 
melted_diffs = melted_jaccard  
melted_diffs$BC_dissimilarity = melted_BC$BC_Dissimilarity 


### convert into similarity
melted_diffs$Jaccard_Similarity = 1 - melted_diffs$Jaccard_Dissimilarity
melted_diffs$BC_Similarity = 1 - melted_diffs$BC_dissimilarity


#### remove half of the dataset
melted_diffs <- melted_diffs %>%
  filter(Sample1 != Sample2) %>%
  mutate(pair_min = pmin(as.character(Sample1), as.character(Sample2)),
         pair_max = pmax(as.character(Sample1), as.character(Sample2))) %>%
  distinct(pair_min, pair_max, .keep_all = TRUE) %>%
  dplyr::select(Sample1, Sample2, Jaccard_Similarity,BC_Similarity )

### Ok I kept only the diagonal

#### Add squirrels -----
melted_diffs$squirrel1 = metadata_table$squirrel_id[match(melted_diffs$Sample1, metadata_table$name)]
melted_diffs$squirrel2 = metadata_table$squirrel_id[match(melted_diffs$Sample2, metadata_table$name)]
melted_diffs$same_squirrel = ifelse(melted_diffs$squirrel1 == melted_diffs$squirrel2, 1,0 )


####  Add mother -----
melted_diffs$mother1 = metadata_table$dam[match(melted_diffs$Sample1, metadata_table$name)]
melted_diffs$mother2 = metadata_table$dam[match(melted_diffs$Sample2, metadata_table$name)]
melted_diffs$same_mother = ifelse(melted_diffs$mother1 == melted_diffs$mother2, 1,0 )
melted_diffs$same_mother = ifelse(is.na(melted_diffs$same_mother),  0,melted_diffs$same_mother )


melted_diffs$mother_offspring = ifelse(melted_diffs$mother1 == melted_diffs$squirrel2, 1,0 )
melted_diffs$mother_offspring = ifelse(melted_diffs$mother2 == melted_diffs$squirrel1, 1,melted_diffs$mother_offspring )
melted_diffs$mother_offspring = ifelse(is.na(melted_diffs$mother_offspring),  0,melted_diffs$mother_offspring )


##### Add litter -----
melted_diffs$litter1 = metadata_table$litter[match(melted_diffs$Sample1, metadata_table$name)]
melted_diffs$litter2 = metadata_table$litter[match(melted_diffs$Sample2, metadata_table$name)]
melted_diffs$same_litter = ifelse(melted_diffs$litter1 == melted_diffs$litter2, 1,0 )
melted_diffs$same_litter = ifelse(is.na(melted_diffs$same_litter),  0,melted_diffs$same_litter )


#### Mast - pup
#### Mast  ------------
melted_diffs$mast1 = metadata_table$mast[match(melted_diffs$Sample1, metadata_table$name)]
melted_diffs$mast2 = metadata_table$mast[match(melted_diffs$Sample2, metadata_table$name)]
melted_diffs =  melted_diffs %>%
  mutate(
    mast_pup = case_when( melted_diffs$mother2 == melted_diffs$squirrel1 & 
                            melted_diffs$mast2  == "mast" ~ "pup_mast",
                          melted_diffs$mother1 == melted_diffs$squirrel2 & 
                            melted_diffs$mast1  == "mast" ~ "pup_mast",
                          melted_diffs$mother2 == melted_diffs$squirrel1 & 
                            melted_diffs$mast2  == "non-mast" ~ "pup_non_mast",
                          melted_diffs$mother1 == melted_diffs$squirrel2 & 
                            melted_diffs$mast1  == "mast" ~ "pup_non_mast", 
                          TRUE ~"other_comparison"# renamed here
                          
    ))

table(melted_diffs$mast_pup)

melted_diffs =  melted_diffs %>%
  mutate(
    masting = case_when( melted_diffs$mast1 == melted_diffs$mast2 & 
                           melted_diffs$mast2 == "mast" & melted_diffs$mast1  == "mast" ~ 1,
                         #melted_diffs$mast2 == "non-mast" & melted_diffs$mast1  == "non-mast" ~ "non-mast",
                         TRUE ~0 # renamed here
                         
    ))

table(melted_diffs$masting)


### Distance in years
melted_diffs$year1 = metadata_table$year[match(melted_diffs$Sample1, metadata_table$name)]
melted_diffs$year2 = metadata_table$year[match(melted_diffs$Sample2, metadata_table$name)]
melted_diffs$year_dist = abs(melted_diffs$year1 - melted_diffs$year2)
table(melted_diffs$year_dist)
melted_diffs$same_year = ifelse(melted_diffs$year1 ==melted_diffs$year2, 1,0 )


### Distance in day of Julyan years (seasonal)
melted_diffs$day_of_year1 = metadata_table$day_of_year[match(melted_diffs$Sample1, metadata_table$name)]
melted_diffs$day_of_year2 = metadata_table$day_of_year[match(melted_diffs$Sample2, metadata_table$name)]
melted_diffs$julian_day_dist <- with(
  melted_diffs,
  pmin(abs(day_of_year1 - day_of_year2),
       365 - abs(day_of_year1 - day_of_year2))
)
hist(melted_diffs$julian_day_dist)

#### Distance in season
melted_diffs$season_1 = metadata_table$season_cat[match(melted_diffs$Sample1, metadata_table$name)]
melted_diffs$season_2 = metadata_table$season_cat[match(melted_diffs$Sample2, metadata_table$name)]
melted_diffs$same_season = ifelse(melted_diffs$season_1 ==melted_diffs$season_2, 1,0 )



### age distance (I have)
melted_diffs$age1 = metadata_table$age[match(melted_diffs$Sample1, metadata_table$name)]
melted_diffs$age2 = metadata_table$age[match(melted_diffs$Sample2, metadata_table$name)]
melted_diffs$age_dist = abs(melted_diffs$age1 - melted_diffs$age2)
table(melted_diffs$age_dist)


### sex 
melted_diffs$sex1 = metadata_table$sex[match(melted_diffs$Sample1, metadata_table$name)]
melted_diffs$sex2 = metadata_table$sex[match(melted_diffs$Sample2, metadata_table$name)]
melted_diffs$same_sex = ifelse(melted_diffs$sex1 == melted_diffs$sex2, 1,0 )
table(melted_diffs$same_sex)

### Geographical distance
grid_stakes_umts = readRDS("grid_stakes_coord_umts.RDS")
metadata_table1 = metadata_table
metadata1_loc <- metadata_table1 %>%
  select(name,locx, locy,gr)


#### X 
metadata1_loc$locx_new = krsp::loc_to_numeric( metadata1_loc$locx )
metadata1_loc$locx_new
metadata1_loc$locx_new <- trunc(metadata1_loc$locx_new)
metadata1_loc$locx_new
#### Y
metadata1_loc$locy_new = krsp::loc_to_numeric( metadata1_loc$locy )
metadata1_loc$locy_new
metadata1_loc$locy_new <- trunc(metadata1_loc$locy_new)
metadata1_loc$locy_new
### Paste togheter
metadata1_loc$locations = paste(metadata1_loc$locx_new, metadata1_loc$locy_new)
metadata1_loc$locations 
metadata1_loc$locations = paste(metadata1_loc$locations, metadata1_loc$gr)
metadata1_loc$locations
grid_stakes_umts$stake_letter_A = krsp::loc_to_numeric( grid_stakes_umts$stake_letter )
grid_stakes_umts$stake_letter_A  ### Ok I can match
grid_stakes_umts$stake_letter_B = krsp::loc_to_numeric( grid_stakes_umts$stake_number )
grid_stakes_umts$stake_letter_B
### Paste togheter
grid_stakes_umts$locations = paste(grid_stakes_umts$stake_letter_A, grid_stakes_umts$stake_letter_B)
grid_stakes_umts$locations 
grid_stakes_umts$locations = paste(grid_stakes_umts$locations, grid_stakes_umts$grid)
grid_stakes_umts$locations
metadata1_loc$east = grid_stakes_umts$easting[match(metadata1_loc$locations  ,grid_stakes_umts$locations )]
metadata1_loc$west = grid_stakes_umts$northing[match(metadata1_loc$locations  ,grid_stakes_umts$locations )]






#### Match distances
### age distance (I have)
melted_diffs$location1_A = metadata1_loc$east[match(melted_diffs$Sample1, metadata1_loc$name)]
melted_diffs$location1_B = metadata1_loc$west[match(melted_diffs$Sample1, metadata1_loc$name)]
melted_diffs$location2_A = metadata1_loc$east[match(melted_diffs$Sample2, metadata1_loc$name)]
melted_diffs$location2_B = metadata1_loc$west[match(melted_diffs$Sample2, metadata1_loc$name)]

melted_diffs$dist_m <- with(
  melted_diffs,
  ifelse(
    is.na(location1_A) | is.na(location1_B) | is.na(location2_A) | is.na(location2_B),
    NA_real_,
    sqrt((location1_A - location2_A)^2 + (location1_B - location2_B)^2)
  )
)

summary(melted_diffs$dist_m ) ### looks about right


### read depth similarity
melted_diffs$depth1 = metadata_table$read_number[match(melted_diffs$Sample1, metadata_table$name)]
melted_diffs$depth2 = metadata_table$read_number[match(melted_diffs$Sample2, metadata_table$name)]
melted_diffs$read_dist = abs(melted_diffs$depth1 - melted_diffs$depth2)



### now lets do grid
melted_diffs$grid1 = metadata_table$gr[match(melted_diffs$Sample1, metadata_table$name)]
melted_diffs$grid2 = metadata_table$gr[match(melted_diffs$Sample2, metadata_table$name)]
melted_diffs$gr = ifelse(melted_diffs$grid1 == melted_diffs$grid2, 1,0 )

tapply(melted_diffs$dist_m, melted_diffs$gr, summary) # Ok great 



### now lets do run
melted_diffs$run1 = metadata_table$run.x[match(melted_diffs$Sample1, metadata_table$name)]
melted_diffs$run2 = metadata_table$run.x[match(melted_diffs$Sample2, metadata_table$name)]
melted_diffs$run = ifelse(melted_diffs$run1 == melted_diffs$run2, 1,0 )



###Normal melted
colnames(melted_diffs)
melted_diffs_selected = melted_diffs 
coluns_to_check_for_NAs = c("Sample1","Sample2", "Jaccard_Similarity", "BC_Similarity", "squirrel1", "squirrel2", "same_squirrel",
                            "mother_offspring", "masting", "year_dist", "same_year", "julian_day_dist", "same_season", "age_dist", "same_sex",
                            "dist_m", "read_dist","run" )
melted_diffs_selected <- melted_diffs_selected[complete.cases(melted_diffs_selected[, coluns_to_check_for_NAs]), ]



colnames(melted_diffs_selected)
#scale all predictors to range between 0-1 if they are not already naturally on that scale
#define scaling function:
range.use <- function(x,min.use,max.use){ (x - min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T)) * (max.use - min.use) + min.use }
scalecols<-c("year_dist", "dist_m", "julian_day_dist", "read_dist", "age_dist")

for(i in 1:ncol(melted_diffs_selected[,which(colnames(melted_diffs_selected)%in%scalecols)])){
  melted_diffs_selected[,which(colnames(melted_diffs_selected)%in%scalecols)][,i]<-range.use(melted_diffs_selected[,which(colnames(melted_diffs_selected)%in%scalecols)][,i],0,1)
}

hist(melted_diffs_selected$same_year)

### save_mother_offspring_pairs
### Make dataset of mother offspring pairs
squirrelID1 = melted_diffs_selected$squirrel1[melted_diffs_selected$mother_offspring == 1]
squirrelID2 = melted_diffs_selected$squirrel2[melted_diffs_selected$mother_offspring == 1]
squirrelID_tot = unique(c(squirrelID1,squirrelID2 ))


melted_diffs_selected_known_mom_pups = melted_diffs_selected[melted_diffs_selected$squirrel1 %in% squirrelID_tot
                                              & melted_diffs_selected$squirrel2 %in% squirrelID_tot,]


melted_diffs_selected1 = melted_diffs_selected
#melted_diffs_selected1 = Melted_diffs_selected_known_mom_pups
#melted_diffs_selected1 = melted_diffs_selected1[which(melted_diffs_selected1$mast1 == melted_diffs_selected1$mast2),]
#melted_diffs_selected1 = melted_diffs_selected1[melted_diffs_selected1$same_year != 1,]
#melted_diffs_selected1 = melted_diffs_selected1[melted_diffs_selected1$same_squirrel != 1,]

table(melted_diffs_selected1$masting, melted_diffs_selected1$same_year)
table(metadata_table$mast, metadata_table$year)


#melted_diffs_selected1 = melted_diffs_selected1[which(melted_diffs_selected1$day_of_year1 < 227  & 
#                                                      melted_diffs_selected1$day_of_year2 < 227),]


#### Can I correct for 
library(ggplot2)




### mast jaccard --------
jaccard_mast <- ggplot(
  melted_diffs_selected1,
  aes(
    x = factor(masting),
    y = Jaccard_Similarity,
    fill = factor(mother_offspring)
  )
) +
  #geom_jitter(
   # aes(color = factor(mother_offspring)),
  #  position = position_jitterdodge(
  #    jitter.width = 0.4,
  #    dodge.width = 0.9
  #  ),
  #  alpha = 0.1,
  #  size = 1
  #) +
  geom_violin(
    position = position_dodge(width = 0.9),
    trim = FALSE,
    alpha = 0.8
  ) +
  geom_boxplot(
    width = 0.15,
    position = position_dodge(width = 0.9),
    alpha = 0.8,
    outlier.shape = NA
  ) +
  scale_x_discrete(
    labels = c(
      "0" = "non-mast",
      "1" = "mast"
    )
  ) +
  scale_fill_discrete(
    labels = c(
      "0" = "not paired",
      "1" = "mother–offspring"
    )
  ) +
  labs(
    x = "Masting",
    y = "Jaccard Similarity",
    fill = "Pair type"
  ) +
  theme_classic() +
  guides(color = "none")+ facet_wrap(~same_year)

jaccard_mast

### mast BC --------
BC_mast <- ggplot(
  melted_diffs_selected1,
  aes(
    x = factor(masting),
    y = BC_Similarity,
    fill = factor(mother_offspring)
  )
) +
  #geom_jitter(
  #  aes(color = factor(mother_offspring)),
  #  position = position_jitterdodge(
  #    jitter.width = 0.4,
  #    dodge.width = 0.9
  #  ),
  #  alpha = 0.1,
  #  size = 1
  #) +
  geom_violin(
    position = position_dodge(width = 0.9),
    trim = FALSE,
    alpha = 0.8
  ) +
  geom_boxplot(
    width = 0.15,
    position = position_dodge(width = 0.9),
    alpha = 0.8,
    outlier.shape = NA
  ) +
  scale_x_discrete(
    labels = c(
      "0" = "non-mast",
      "1" = "mast"
    )
  ) +
  scale_fill_discrete(
    labels = c(
      "0" = "not paired",
      "1" = "mother–offspring"
    )
  ) +
  labs(
    x = "Masting",
    y = "BC Similarity",
    fill = "Pair type"
  ) +
  theme_classic() +
  guides(color = "none") + facet_wrap(~same_year)

BC_mast



### season jaccard--------
jaccard_season <- ggplot(
  melted_diffs_selected1,
  aes(
    x = factor(same_season),
    y = Jaccard_Similarity,
    fill = factor(mother_offspring)
  )
) +
  #geom_jitter(
  #  aes(color = factor(mother_offspring)),
  #  position = position_jitterdodge(
  #    jitter.width = 0.4,
  #    dodge.width = 0.9
  #  ),
  #  alpha = 0.1,
  #  size = 1
  #) +
  geom_violin(
    position = position_dodge(width = 0.9),
    trim = FALSE,
    alpha = 0.8
  ) +
  geom_boxplot(
    width = 0.15,
    position = position_dodge(width = 0.9),
    alpha = 0.8,
    outlier.shape = NA
  ) +
  scale_x_discrete(
    labels = c(
      "0" = "different",
      "1" = "same"
    )
  ) +
  scale_fill_discrete(
    labels = c(
      "0" = "not paired",
      "1" = "mother–offspring"
    )
  ) +
  labs(
    x = "Season",
    y = "Jaccard Similarity",
    fill = "Pair type"
  ) +
  theme_classic() +
  guides(color = "none") #+ facet_wrap(~same_year)

jaccard_season

### BC season ---------------
BC_season <- ggplot(
  melted_diffs_selected1,
  aes(
    x = factor(same_season),
    y = BC_Similarity,
    fill = factor(mother_offspring)
  )
) +
  #geom_jitter(
  #  aes(color = factor(mother_offspring)),
  #  position = position_jitterdodge(
  #    jitter.width = 0.4,
  #    dodge.width = 0.9
  #  ),
  #  alpha = 0.1,
  #  size = 1
  #) +
  geom_violin(
    position = position_dodge(width = 0.9),
    trim = FALSE,
    alpha = 0.8
  ) +
  geom_boxplot(
    width = 0.15,
    position = position_dodge(width = 0.9),
    alpha = 0.8,
    outlier.shape = NA
  ) +
  scale_x_discrete(
    labels = c(
      "0" = "different",
      "1" = "same"
    )
  ) +
  scale_fill_discrete(
    labels = c(
      "0" = "not paired",
      "1" = "mother–offspring"
    )
  ) +
  labs(
    x = "Season",
    y = "BC Similarity",
    fill = "Pair type"
  ) +
  theme_classic() +
  guides(color = "none") + facet_wrap(~same_year)
BC_season

library(ggpubr)
ggarrange(BC_season , jaccard_season,BC_mast, jaccard_mast, common.legend = TRUE)


cor.test(melted_diffs_selected1$Jaccard_Similarity,melted_diffs_selected1$year_dist, method = "pearson")
unique(melted_diffs_selected1$run1)
melted_diffs_selected1$technology1 = ifelse(melted_diffs_selected1$run1 %in% c("Old_run1", "Old_run2"), "Old_run", "New_run")
melted_diffs_selected1$technology2 = ifelse(melted_diffs_selected1$run2 %in% c("Old_run1", "Old_run2"), "Old_run", "New_run")
melted_diffs_selected1$technology = ifelse(melted_diffs_selected1$technology1 ==  melted_diffs_selected1$technology2, 1, 0)
melted_diffs_selected2 = melted_diffs_selected1[melted_diffs_selected1$same_year != 1,]


cor.test(melted_diffs_selected2$BC_Similarity,melted_diffs_selected2$year_dist, method = "pearson")

# Create the plot
ggplot(melted_diffs_selected2, aes(x = year_dist, y = BC_Similarity)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
 
  labs(x = "Year Distance", 
       y = "BC Similarity",
       title = "BC Similarity vs Year Distance, zero escluded") +
  theme_bw() + facet_wrap(~run)


