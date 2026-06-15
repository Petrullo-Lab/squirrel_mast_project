#### Merge ASVtabs
#### Setwd
setwd("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/Squirrels/analysis_data_oct_2025/16s/comb_ASV_tabs")
library(vegan)
library(dada2)
library(dplyr)
library(stringr)
library(lubridate)
library(ggplot2)
library(ggpubr)



### Old ASVtabs
ASVtab_run1 = readRDS("asv_table_old1.RDS")
ASVtab_run2 = readRDS("asv_table_old2.RDS")
merged_old <- mergeSequenceTables(ASVtab_run1, ASVtab_run2)
sort(rowSums(merged_old)) # 665 samples
grep("Water", rownames(merged_old))
length(which(rowSums(merged_old) > 10000)) ### 609 with also 7 blanks... ok



### new data
ASVtab_new = readRDS("asv_table.RDS")
rownames(ASVtab_new)
ASVtab_new_blnaks = ASVtab_new[grep("Blank|blank", rownames(ASVtab_new)),] # four blanks all few reads ok 
ASVtab_new_Zymo = ASVtab_new[grep("zymo|Zymo|zyno|Zymo", rownames(ASVtab_new)),] # four blanks all few reads ok 
ASVtab_new_Zymo_nonzero <- ASVtab_new_Zymo[, colSums(ASVtab_new_Zymo != 0) > 0]
ASVtab_zymo = ASVtab_new[,colnames(ASVtab_new) %in% colnames(ASVtab_new_Zymo_nonzero)]
ASVtab_new_Zymo <- ASVtab_zymo[order(grepl("zymo|Zymo|zyno|Zymo", rownames(ASVtab_zymo), ignore.case = TRUE), decreasing = TRUE), ]

rowSums(ASVtab_new_blnaks) # 237 6103 20 1671 reads
rowSums(ASVtab_new_blnaks != 0)
ASVtab_new <- ASVtab_new[!grepl("WH|MGRS|MG|FUR|Fur|Flea|ABF|PBF|swabs|BF", rownames(ASVtab_new)), ]
rownames(ASVtab_new) ### good
ASVtab_new_blnaks = ASVtab_new[grep("Blank|blank", rownames(ASVtab_new)),] # four blanks all few reads ok 
rowSums(ASVtab_new_blnaks) # 237 6103 20 1671 reads
rownames(ASVtab_new)
nrow(ASVtab_new)
length(which(rowSums(ASVtab_new) > 10000)) ### 148 / 154 have more than 10000 reads --- good
ASVtab_new <- ASVtab_new[, colSums(ASVtab_new) > 0]

### Merge old and new
merged_fin <- mergeSequenceTables(merged_old, ASVtab_new)
ncol(ASVtab_new)
ncol(merged_old)
ncol(merged_fin)
sum(3678,6745) # more or less 2200 ASVs were merged ### not so many
sort(rowSums(merged_fin))
merged_fin = merged_fin[rowSums(merged_fin) > 10000,]
sort(rowSums(merged_fin))

#### create a metadata
meta_tot = data.frame(name = row.names(merged_fin)) %>%
  mutate(
    run = case_when(
      name %in% rownames(ASVtab_run1) ~ "Old_run1",
      name %in% rownames(ASVtab_run2) ~ "Old_run2",
      name %in% rownames(ASVtab_new) ~ "New_run"
    ),
    sample_type = case_when(
      str_detect(name, regex("blank|water", ignore_case = TRUE)) ~ "blank",
      str_detect(name, regex("zymo|zyno", ignore_case = TRUE)) ~ "zymo",
      TRUE ~ "sample"
    )
  )  

unique(meta_tot$sample_type)

meta_tot = meta_tot[meta_tot$sample_type != "zymo",]
merged_fin = merged_fin[rownames(merged_fin) %in% meta_tot$name,]


meta_tot$Richness.rar<-specnumber(rrarefy(merged_fin, 10000)) #Change this number according to the lowest sample
meta_tot$Shannon.rar<-vegan::diversity(rrarefy(merged_fin, 10000), index="shannon")#Change this number according to the lowest sample



a = ggplot(meta_tot, aes(run, Richness.rar))+
  geom_boxplot( alpha=0.3, outlier.shape = NA)+ geom_jitter()+
  theme_classic() +ylab("ASV richness")
a
b = ggplot(meta_tot, aes(run, Shannon.rar))+
  geom_boxplot( alpha=0.3, outlier.shape = NA)+geom_jitter()+
  theme_classic() +ylab("Shannon H'")
b
ggarrange(a,b)



ASVtab_16s.rar = rrarefy(merged_fin, 10000) #78000 is the number decided
bac.dist<-vegdist(ASVtab_16s.rar, method="bray")
set.seed(12)
permanova = adonis(ASVtab_16s.rar ~ run, data= meta_tot)
permanova$aov.tab





bac.nmds<-metaMDS(bac.dist, k=2, try=100)
bac.nmds$stress #0.06
stress= paste("stress = ", round(bac.nmds$stress, digits = 2))

meta_tot$Axis01<-bac.nmds$points[,1]
meta_tot$Axis02<-bac.nmds$points[,2]



bac_nmds<-ggplot(meta_tot, aes(x=Axis01, y=Axis02))+
  geom_point(aes(color=run, shape=sample_type), size = 3)+
  #geom_mark_hull(aes(group=sample_cat, label = sample_cat), concavity=10) +
  scale_color_brewer(palette = "Dark2")+
  theme_classic()+
  annotate("text",x=min(meta_tot$Axis01)+ 0.17,y=max(meta_tot$Axis02),hjust=1,label= stress)+
  ggtitle("Old runs vs new")
bac_nmds



meta_tot$seq_depth = rowSums(merged_fin) 
meta_tot = meta_tot[meta_tot$seq_depth<1000000,]
### Seq depth

a = ggplot(meta_tot, aes(run, seq_depth))+
  geom_boxplot( alpha=0.3, outlier.shape = NA)+ geom_jitter()+
  theme_classic() +ylab("ASV richness")
a
ggplot(meta_tot, aes(x = seq_depth, y = Richness.rar, color = run)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    x = "Sequencing depth",
    y = "Richness (rarefied)",
    color = "Run",
    title = "Correlation between rarefied richness and sequencing depth"
  ) +
  theme_minimal(base_size = 14)
by(meta_tot, meta_tot$run, function(x) {
  cor.test(x$Richness.rar, x$seq_depth, method = "pearson")
})



#### Construct_metadata
# 1st_old samples
samdata_old_run = readRDS("samdata_old_run.RDS")
meta_tot$name_old <- sub("_.*", "", meta_tot$name)
meta_tot$name_old %in% rownames(samdata_old_run)
samdata_old_run  = samdata_old_run[rownames(samdata_old_run) %in%meta_tot$name_old, ]
which(rownames(samdata_old_run) %in% meta_tot$name_old)
meta_old = meta_tot[which(meta_tot$name_old %in%rownames(samdata_old_run)),] ### they are all there
samdata_old_run$name_old = rownames(samdata_old_run)
old_meta_combined <- merge(meta_old, samdata_old_run, by = "name_old", all = TRUE)

colnames(old_meta_combined)
old_meta_combined <- old_meta_combined[, c("name_old", "name", "run.x", "sample_type", "trapping_id")]






# now new samples
meta_new_run = meta_tot[meta_tot$run == "New_run",]
meta_new_run$extracted <- str_extract(meta_new_run$name, "(?<=16s_)[^_]+(?=_S)")
meta_new_run$extracted %in% poop_link$poop_id_new
meta_new_run <- meta_new_run %>%
  mutate(
    extracted_padded = str_replace(extracted, "^(\\d{2})(?=-)", "0\\1"),   # pad first number (2 digits → 3)
    extracted_padded = str_replace(extracted_padded, "(?<=-)(\\d{1})(?!\\d)", "\\10")  # pad second number (1 digit → add trailing 0)
  )

### Now I need to get the metadata for these samples
meta_new_run$extracted_padded

# Get KRSP poops
#Connect to KRSP database
### Pedigree and litter download -----------------
library(krsp)
con <- krsp_connect(host = "krsp.cepb5cjvqban.us-east-2.rds.amazonaws.com",
                     dbname ="krsp",
                     username = Sys.getenv("krsp_user"),
                     password = Sys.getenv("krsp_password")
)
# Connect to the KRSP tables
krsp_tables(con)
krsp_poop_link = tbl(con, "trapping") %>%
  collect()
#saveRDS(krsp_poop_link, "krsp_poop_link.RDS")
krsp_poop_link$date = as.Date(krsp_poop_link$date)
krsp_poop_link$year  = year(krsp_poop_link$date)
krsp_poop_link = krsp_poop_link[krsp_poop_link$year > 2010,]
krsp_poop_link$poop_id = gsub("\\.", "-",krsp_poop_link$tvariable1)

patterns <- c("K1094", "K2824", "K5314")

for (p in patterns) {
  matched_rows <- grepl(p, krsp_poop_link$comments, ignore.case = FALSE)
  krsp_poop_link$poop_id[matched_rows] <- p
}


krsp_poop_link_1 = krsp_poop_link[!is.na(krsp_poop_link$poop_id),] %>%
  mutate(
    extracted_padded = str_replace(poop_id, "^(\\d{1})(?=-)", "0\\1"),  
    extracted_padded = str_replace(extracted_padded, "^(\\d{2})(?=-)", "0\\1"), # pad first number (2 digits → 3)
    extracted_padded = str_replace(extracted_padded, "(?<=-)(\\d{1})(?!\\d)", "\\10")  # pad second number (1 digit → add trailing 0)
  )
length(which(krsp_poop_link_1$extracted_padded %in% meta_new_run$extracted_padded)) ### 136 we have
'%nin%' = Negate("%in%")
nrow(meta_new_run)
meta_new_run$extracted_padded[which(meta_new_run$extracted_padded %nin% krsp_poop_link_1$extracted_padded )]


### now I have to do it for the new sample:
krsp_poop_link_1_to_match = krsp_poop_link_1[which(krsp_poop_link_1$extracted_padded  %in%  meta_new_run$extracted_padded),]
meta_new_run_to_match  = meta_new_run[which(meta_new_run$extracted_padded %in% krsp_poop_link_1_to_match$extracted_padded ),]
new_meta_combined1 <- merge(meta_new_run_to_match,krsp_poop_link_1_to_match,  by = "extracted_padded", all = TRUE)
colnames(new_meta_combined1)
old_meta_combined <- old_meta_combined[, c("name_old", "name", "run.x", "sample_type", "trapping_id")]
new_meta_combined = data.frame(name_old = new_meta_combined1$extracted_padded,
                               name = new_meta_combined1$name,
                               run.x = new_meta_combined1$run,
                               sample_type = new_meta_combined1$sample_type,
                               trapping_id = new_meta_combined1$id
                               )
### rbind the two tables
metadata_combined = rbind(old_meta_combined, new_meta_combined)
length(unique(metadata_combined$trapping_id))
dupes <- metadata_combined$trapping_id[duplicated(metadata_combined$trapping_id)]
unique(dupes)



metadata_combined$id = metadata_combined$trapping_id
krsp_poop_link = tbl(con, "trapping") %>%
  collect()
krsp_poop_link_match = krsp_poop_link[krsp_poop_link$id %in% metadata_combined$id,]
new_meta_combined_final <- merge(metadata_combined,krsp_poop_link_match,  by = "id", all = TRUE)
saveRDS(new_meta_combined_final,"new_meta_all_samples.RDS")
merged_fin = merged_fin[rownames(merged_fin) %in% new_meta_combined_final$name,]
saveRDS(merged_fin,"new_ASVtable_all_samples.RDS")
merged_fin = merged_fin[m]
colnames(merged_fin)[1:10]
ncol(merged_fin)




#### New taxonomy table = 
TAXtab16sold1 = readRDS("taxa_table_silva_run1.RDS")
TAXtab16sold2 = readRDS("taxa_table_silva_run2.RDS")
TAXtab16snew = readRDS("taxa_table_silva_new.RDS")

'%nin%' = Negate("%in%")
Taxtab_total = rbind(TAXtab16sold1,TAXtab16sold2)
Taxtab_total = rbind(Taxtab_total,TAXtab16snew)
Taxtab_total <- Taxtab_total[!duplicated(rownames(Taxtab_total)), ]
ASVtable = readRDS("new_ASVtable_all_samples.RDS")
ncol(ASVtable)
nrow(Taxtab_total)
unique(colnames(ASVtable)  %nin% rownames(Taxtab_total)) ## all there, just subset
Taxtab_total_final = Taxtab_total[rownames(Taxtab_total) %in% colnames(ASVtable),]
nrow(Taxtab_total_final)
Taxtab_total_final <- Taxtab_total_final[colnames(ASVtable), ]
unique(rownames(Taxtab_total_final) == colnames(ASVtable))

seqs <- colnames(merged_fin)
lens <- nchar(seqs)
table(lens)

saveRDS(Taxtab_total_final,"Taxtab_total_silva.RDS")










