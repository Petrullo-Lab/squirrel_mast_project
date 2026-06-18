# Code for mast project paper

This is a guide for the analysis presented in the paper  "Maternal effects shape scale-dependent convergence in gut microbial responses to environmental change"

To start process the downaloded fastq files present at the accession number 

then run the following scripts:

1) Bioinformatics/16s_dada2.R On each separate run file this runs. Dada2 and creates ASVtables
2) Bioinformatics/combination_runs.R to cobine the different runs.  This code is generated through the KRSP database which has limited access. Please contact me for a copy of the data used if interested in repeating the analysis

Neverthless I provide here the data to perform all statistical analysis, in the folder data.

Then you can move to statistics

1) Alpha_beta_diff_abb_analysis.R for alpha beta and differential abudance analysis
2) dietary_analysis.R for the dietary analysis
3) picust_normal_analysis_masting.R for the analysis of the picrust predictions


In the folder "BRMS" there are all codes to run BRMS models

1) BRMS_data_preparation.R to prepare the metadata for the analysis
2) rarefied_brms.R to run the models (in an HPC)
3) model_analysis_BRMS_final.R to generate BRMS plots
4) genera_metadata_preparation.R to prepare the metadata for the genera analysis
5) code_genera_new.R to run the genera analysis
6) BRMS_genera_plot.R to plot them




