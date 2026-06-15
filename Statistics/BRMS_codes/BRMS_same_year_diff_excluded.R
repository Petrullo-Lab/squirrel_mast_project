library(brms)




options(brms.backend = "cmdstanr")
setwd("/xdisk/laurenpetrullo/schiro/dyadic_models/new_subselection")
melted_diffs_selected_try = readRDS("Only_same_year_15th_Aug_BRMS.RDS")


### first Same year models (remove different years) ------
### Model 1 jacc 
model1<-brm(Jaccard_Similarity~1+mother_offspring*masting+mother_offspring*same_season+
            +dist_m+same_sex+read_dist+run+age_dist+
            (1|mm(Sample1,Sample2)) + (1|mm(squirrel1,squirrel2)),
            data = melted_diffs_selected_try,
            family= "Beta",
            warmup = 1000, iter = 3000,
            threads = threading(6),
            cores = 24,
            chains = 4,
            init=0)


saveRDS(model1,"model1_jacc_same_year_ag_15th_diff_excluded.RDS")




### Model 2 BC sim
model1<-brm(BC_Similarity~1+mother_offspring*masting+mother_offspring*same_season+
              dist_m+same_sex+read_dist+run+age_dist+
              (1|mm(Sample1,Sample2)) + (1|mm(squirrel1,squirrel2)),
            data = melted_diffs_selected_try,
            family= "Beta",
            warmup = 1000, iter = 3000,
            threads = threading(6),
            cores = 24,
            chains = 4,
            init=0)

saveRDS(model1,"model1_BC_same_year_ag_15th_diff_excluded.RDS")



### first Same year models (remove different years) ------
### Model 1 jacc 
model1<-brm(Jaccard_Similarity~1+masting+same_season+
              +dist_m+same_sex+read_dist+run+age_dist+
              (1|mm(Sample1,Sample2)) + (1|mm(squirrel1,squirrel2)),
            data = melted_diffs_selected_try,
            family= "Beta",
            warmup = 1000, iter = 3000,
            threads = threading(6),
            cores = 24,
            chains = 4,
            init=0)


saveRDS(model1,"model1_jacc_same_year_ag_15th_diff_excluded_no_mom_off.RDS")




### Model 2 BC sim
model1<-brm(BC_Similarity~1+mother_offspring*masting+mother_offspring*same_season+
              dist_m+same_sex+read_dist+run+age_dist+
              (1|mm(Sample1,Sample2)) + (1|mm(squirrel1,squirrel2)),
            data = melted_diffs_selected_try,
            family= "Beta",
            warmup = 1000, iter = 3000,
            threads = threading(6),
            cores = 24,
            chains = 4,
            init=0)

saveRDS(model1,"model1_BC_same_year_ag_15th_diff_excluded_no_mom_off.RDS")

