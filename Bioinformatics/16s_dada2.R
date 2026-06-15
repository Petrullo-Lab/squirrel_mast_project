# Set library path
.libPaths("/contrib/laurenpetrullo/dada2/lib/R/library")

# Load required packages
library("dada2")
library("ggpubr")

args <- commandArgs(trailingOnly = TRUE)

# Ensure both arguments are provided
if (length(args) < 2) {
  stop("Please provide both the input directory and the output directory.")
}

# Set input and output directories
input_dir <- args[1]  # First argument: input directory
output_dir <- args[2]  # Second argument: output directory

# Set working directory to the output directory
setwd(output_dir)
path = input_dir

# Print paths for verification
print(paste("Input directory:", input_dir))
print(paste("Output directory set to:", getwd()))

#list the files
fnFs <- sort(list.files(path, pattern="_R1_", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_R1_"), `[`, 1)

#Check for matching lengths
if ( length(fnFs) != length(fnRs)){ stop("Number of R1 and R2 files is different") } 

# Check if the samples are less than six, in case just print all, otherwise print 6 random samples
if ( length(fnFs) < 6)
{ 
  random_numbers <- c(1:length(fnFs))
} else
{
  random_numbers <- sample(1:length(fnFs), 6, replace = FALSE)
}



a = plotQualityProfile(fnFs[1]) # forward
ggsave(filename = "qualityF.pdf", a)
a = plotQualityProfile(fnRs[1]) # reverse
a
ggsave(filename = "qualityR.pdf", a)
print("Quality profiles printed")

# Filter the samples
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(145,145),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=24, trimLeft = 5) # On Windows set multithread=FALSE
out = as.data.frame(out)
out$Retained = out$reads.out/out$reads.in*100


filtFs = filtFs[out$reads.out > 0]
filtRs = filtRs[out$reads.out > 0]
sample.names = sample.names[out$reads.out > 0]

# Learn errors
errF <- learnErrors(filtFs, multithread=24)
saveRDS(errF, "filtFs.RDS")

errR <- learnErrors(filtRs, multithread=24)
saveRDS(errR, "errR.RDS")
# Plot error rates that are based on nucleotide transition probabilities and quality scores; observed scores (black) should match red lines (expected)
a = plotErrors(errF, nominalQ=TRUE)
ggsave("plotErrorsF.pdf", a)
a = plotErrors(errR, nominalQ=TRUE)
ggsave("plotErrorsR.pdf", a)
# Reduce computation time by dereplicating sequences
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)




# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names
# Sample inference
dadaFs <- dada(derepFs, err=errF, multithread=24)
saveRDS(dadaFs, "dadaFs.RDS")
dadaRs <- dada(derepRs, err=errR, multithread=24)
saveRDS(dadaRs, "dadaRs.RDS")

### Merge etc
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=24, verbose=TRUE)
saveRDS(seqtab.nochim,"asv_table.RDS")


getN <- function(x) sum(getUniques(x))
new_df = data.frame( sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))


rownames(out) = gsub("_R1_001.fastq.gz", "",rownames(out))
track <- merge(out, new_df, by=0, all = T)
track[is.na(track)] = 0
track1 = data.frame(track)

# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track1) <- c("sample","input", "after_filt", "quality filt Retained %","denoised", "merged", "nonchim")
#Total percent of dataset remaining
track1$Tot_Retained = track1$nonchim/track1$input*100
head(track1)





#### save all tables
write.table(track1, "sequence_pipeline_stats.txt", sep = "\t", quote = F)
write.table(seqtab.nochim, "asv_table.txt", sep = "\t", quote = F)






### Assign taxonomy
taxaGTDBfull <- assignTaxonomy(seqtab.nochim, "/groups/laurenpetrullo/ref_databases/GTDB_bac120_arc53_ssu_r207_Genus.fa.gz", multithread=TRUE)
taxaGTDBfull <- addSpecies(taxaGTDBfull, "/groups/laurenpetrullo/ref_databases/GTDB_bac120_arc53_ssu_r207_Species.fa.gz")
taxaS138 <- assignTaxonomy(seqtab.nochim, "/groups/laurenpetrullo/ref_databases/silva_nr99_v138.1_train_set.fa", multithread=TRUE)
taxaS138 <- addSpecies(taxaS138, "/groups/laurenpetrullo/ref_databases/silva_species_assignment_v138.1.fa")
#save tables in RDS format 
saveRDS(taxaGTDBfull, "taxa_tableGTDB.RDS")
saveRDS(taxaS138, "taxa_table_silva.RDS")
#### save all tables
write.table(taxaGTDBfull, "taxa_tableGTDB.txt", sep = "\t", quote = F)
write.table(taxaS138, "taxa_table_silva.txt", sep = "\t", quote = F)
write.table(seqtab.nochim, "asv_table.txt", sep = "\t", quote = F)



