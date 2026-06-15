# code to run ordbetareg on my computer
cat("R script started/n")
flush.console()
library(brms)
library(ordbetareg)
library(marginaleffects)
cat("brms loaded/n")
flush.console()

library(cmdstanr)
cat("cmdstanr loaded/n")
flush.console()


Sys.setenv(OMP_NUM_THREADS = 6)
cat("Sys.setenv set to 6 threads/n")
flush.console()

setwd("/xdisk/laurenpetrullo/schiro/dyadic_models/ord_beta_reg")
cat("setwd done/n")
flush.console()

# Get the array task ID
args <- commandArgs(trailingOnly = TRUE)
index <- as.integer(args[1])
cat("args done/n")
flush.console()

# Load ASV list and names
asv_list <- readRDS("BC_aligned_by_genus_filtered_new.RDS")
asv_names <- readRDS("genera_names_new.rds")
cat("readRDSs done")
flush.console()

# Safety check
if (index > length(asv_list)) {
  stop("Index out of bounds: ", index)
}

# Extract data and ASV name
df <- asv_list[[index]]
asv <- asv_names[index]
print("data extracted")

# Make sure output folder exists
if (!dir.exists("output_models_BC")) dir.create("output_models_BC")

print("model started")
# Fit the dyadic model
priors <- c(
  set_prior("normal(1, 0.5)", class="b", coef="mother_offspring"),
  set_prior("normal(1, 0.5)", class="b", coef="mother_offspring:masting"),
  set_prior("normal(1, 0.5)", class="b", coef="mother_offspring:same_season")
)

model <- ordbetareg(
  BC_Similarity ~ 1 + mother_offspring*masting_factor + mother_offspring*same_season +
    year_dist + same_year + dist_m + same_sex + run + age_dist +
    (1 | mm(Sample1, Sample2)) + (1 | mm(squirrel1, squirrel2)),
  data = df,
  extra_prior = priors,
  warmup = 1000, iter = 3000, chains = 4, cores = 4,
  threads = threading(6), init = 0,
  control = list(adapt_delta = 0.95),
  backend = "cmdstanr"
)

# Save as RDS as well
saveRDS(model, file = paste0("output_models_BC/brms_model_BC_", asv, ".rds"))

