library(brms)




options(brms.backend = "cmdstanr")
setwd("/xdisk/laurenpetrullo/schiro/dyadic_models/new_subselection")
melted_diffs_selected_try = readRDS("melted_diffs_selected_mast.RDS")


### first Same year models (remove different years) ------
### Model 1 jacc 
model1<-brm(Jaccard_Similarity~1+mother_offspring+same_season+year_dist+
              dist_m+same_sex+read_dist+run+age_dist+
              (1|mm(Sample1,Sample2)) + (1|mm(squirrel1,squirrel2)),
            data = melted_diffs_selected_try,
            family= "Beta",
            warmup = 1000, iter = 3000,
            threads = threading(6),
            cores = 24,
            chains = 4,
            init=0)


saveRDS(model1,"model1_jacc_only_mast_non_rar.RDS")




### Model 2 BC sim
model1<-brm(BC_Similarity~1+mother_offspring+same_season+year_dist+
              dist_m+same_sex+read_dist+run+age_dist+
              (1|mm(Sample1,Sample2)) + (1|mm(squirrel1,squirrel2)),
            data = melted_diffs_selected_try,
            family= "Beta",
            warmup = 1000, iter = 3000,
            threads = threading(6),
            cores = 24,
            chains = 4,
            init=0)

saveRDS(model1,"model1_BC_only_mast_non_rar.RDS")



###  model 3 rar jacc ------

model1<-brm(Jaccard_Similarity_rar~1+mother_offspring+same_season+year_dist+
              dist_m+same_sex+read_dist+run+age_dist+
              (1|mm(Sample1,Sample2)) + (1|mm(squirrel1,squirrel2)),
            data = melted_diffs_selected_try,
            family= "Beta",
            warmup = 1000, iter = 3000,
            threads = threading(6),
            cores = 24,
            chains = 4,
            init=0)


saveRDS(model1,"model1_jacc_only_mast_rar.RDS")




### Model 24BC sim rar ------------
model1<-brm(BC_Similarity_rar~1+mother_offspring+same_season+year_dist+
              dist_m+same_sex+read_dist+run+age_dist+
              (1|mm(Sample1,Sample2)) + (1|mm(squirrel1,squirrel2)),
            data = melted_diffs_selected_try,
            family= "Beta",
            warmup = 1000, iter = 3000,
            threads = threading(6),
            cores = 24,
            chains = 4,
            init=0)

saveRDS(model1,"model1_BC_only_mast_rar.RDS")
