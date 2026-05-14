# -----------------------------------------------------------------------
# Title: Estimate time-varying interaction coefficients (script to run function) 
# Author: Ewa Merz (e2merz@ucsd.edu)
# Description: This script is used to run a function to estimate the interactions between microbes at the Scripps Pier using 16S and 18S sequencing data. The output of the function are the following folders: 

## 1. multivariate_embeddings: Contains the created random multivariate embeddings (number defined by MV_embeddings) and their predictive performance using multivariate simplex projection
## 2. distance_matrix: Contains the multiview distance matrix for each target (e.g., taxa, species, etc.) created by averaging the distances between sampling points from the best performing (number defined by top_embeds) multivariate embeddings from 1.
## 3. parametrization: Contains the parametrization results from tuning the regularized regression using the defined values for alpha, lambda and theta (seq_alpha, seq_lambda and seq_theta) for each target
## 4. coefficients: Estimated interaction coefficients from the best performing parameter combinations defined in 3. for each target
## 5. performance: Contains the model performance for each target of the best parameter combinations used in 4. to estimate the interaction coefficients
## 6. interactions: Averaged interaction coefficients from 4. for each target, weighted based on performance in 4.

# R version: R version 4.4.0 (2024-04-24) -- "Puppy Cup"
# -----------------------------------------------------------------------

## 1. Required libraries
library(parallel) # used to parallelize processes

## 2. Load function
# Sources the custom function script to estimate interactions.
source("scripts/functions/interaction_coefficients_function_rEDM_1.15.3_pred_tp.R") # this function is used to estimate the interaction coefficients

## 3. Read data
# CSV containing relative abundance of 16S/18S ASVs and temperature.
data <- read.csv("data/data_sequences_0.1_rel_ab_0.5_occ_binned_4_days_with_temperature.csv") # relative abundance of microbial community and water temperature

## 4. Target variables
# Identify all ASVs (columns) for which we want to estimate interactions.
# We exclude the 'date' column to leave only biological and environmental variables.
target_vars <- names(select(data,-date)) # the target variables for which interactions should be calculated (here the microbial taxa)

## 5. Run function (with the option to parallelize)
# Use mclapply to iterate the inter_coeffs function over the list of target ASVs.
# Note for local testing: To test the script quickly, you can set target_vars to a 
# single ASV name and reduce the seq_parameters to single values.

mclapply(X=target_vars, # target (e.g. species or taxon)
         FUN=inter_coeffs, # interaction function, previously loaded
         input_dat=data, # data (first column is date, the rest species or sequence relative or absolute abundance)
         out_dir="model_out/ASVs_Scripps_Pier_28092024/", # output directory (where should the model output be saved?)
         
         # --- Multiview Embedding Parameters ---
         E=10, # embedding dimension used for simplex projection and distance calculations
         MV_embeddings=1000, # number of random multivariate embeddings created
         top_embeds=100, # number of best multivariate embeddings used to calculate the multiview distance
         
         # --- Hyperparameter Grid for S-map and Elastic Net ---
         seq_theta=c(0,0.1,1,3,8), # sequence of theta (nonlinearity) for parameter tuning
         seq_alpha=c(0.1,0.3,0.5,0.7,0.9), # sequence of alpha (relationship of L2 and L1 penalty) for parameter tuning
         seq_lambda=c(1.000, 0.178, 0.032, 0.006, 0.001), # sequence of lambda (penalization of regression coefficients) for parameter tuning
         
         # --- Prediction Settings ---
         pred_tp=2, # prediction horizon (here 2 time points = 8 days)
         
         # --- Parallelization Settings ---
         mc.cores=1) # number of cores to use (i.e., when running on a server)