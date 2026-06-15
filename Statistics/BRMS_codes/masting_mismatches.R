library(brms)




options(brms.backend = "cmdstanr")
setwd("/xdisk/laurenpetrullo/schiro/dyadic_models/masting_dfferent_code")
melted_diffs_selected_try = readRDS("Full_data_set_BRMS.RDS")


###  model 3 rar jacc ------
model1<-brm(BC_Similarity_rar~1+mother_offspring*masting_1+mother_offspring*same_season+year_dist+same_year+
              dist_m+same_sex+run+age_dist+
              (1|mm(Sample1,Sample2)) + (1|mm(squirrel1,squirrel2)),
            data = melted_diffs_selected_try,
            family= "Beta",
            warmup = 1000, iter = 3000,
            threads = threading(6),
            cores = 24,
            chains = 4,
            init=0)


saveRDS(model1,"model1_BC_rarefied_masting_mismatches.RDS")



###  model 3 rar jacc ------

model1<-brm(Jaccard_Similarity_rar~1+mother_offspring*masting_1+mother_offspring*same_season+year_dist+same_year+
              dist_m+same_sex+run+age_dist+
              (1|mm(Sample1,Sample2)) + (1|mm(squirrel1,squirrel2)),
            data = melted_diffs_selected_try,
            family= "Beta",
            warmup = 1000, iter = 3000,
            threads = threading(6),
            cores = 24,
            chains = 4,
            init=0)


saveRDS(model1,"model1_jacc_rarefied_mismatches.RDS")
