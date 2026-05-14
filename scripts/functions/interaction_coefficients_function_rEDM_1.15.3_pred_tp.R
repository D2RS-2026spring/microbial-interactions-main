# -----------------------------------------------------------------------
# Title: Estimate time-varying interaction coefficients with the MDR S-map (function) 
# Author: Ewa Merz (e2merz@ucsd.edu)
# Description: This script is used to create a function that can estimate the interactions between taxa or species within a community using the MDR S-map approach described by Chang et al. (2021): Reconstructing large interaction networks from empirical time series data. 

# The analysis encompasses several steps:
# 1. Scale and lag the data
# 2. Create random multivariate embeddings
# 3. Evaluate forecast performance using simplex of the in 2. created multivariate embeddings
# 4. Calculate the multiview distance between points in state-space based on the best performing multivariate embeddings from 4.
# 5. MVD S-map parametrization: Find the best parameter combination of theta, lambda and alpha for the regularized S-maps using "leave-future-out cross-validation".
# 6. Estimate interaction coefficients based on the best performing MDR S-map models from 6.
# 7. Average interaction coefficients

# R version: R version 4.4.0 (2024-04-24) -- "Puppy Cup"
# -----------------------------------------------------------------------

options(scipen = 6, digits = 4) # controls the digits of numeric variables

# -----------------------------------------------------------------------

## Required libraries
library(tidyverse) # awesome to reshape and wrangle data
library(lubridate) # great for dealing with dates
library(rEDM) # used for multivariate simplex projection
library(vegan) # used to calculate the distance between points
library(glmnet) # run regularized regression (elastic-net)
library(zoo) # used to interpolate values before Simplex (and only for simplex)

# -----------------------------------------------------------------------

## Define helper functions
liner_interp <- function(x){ # linear interpolation to approximate NA's (for simplex only)
  # EDM methods are sensitive to gaps; thus we approximate them using linear interpolation. However, if the gap is larger than 4, we don't interpolate.
  x <- na.approx(x,na.rm=F,maxgap=4)
}

# -----------------------------------------------------------------------

## Define function to estimate interaction coefficients
inter_coeffs <- function(target_i, # target (e.g. species or taxon)
                         input_dat, # data (first column is date, the rest species or sequence relative or absolute abundance)
                         out_dir, # output directory (where should the model output be saved?)
                         E=10,
                         MV_embeddings=1000, # number of random multivariate embeddings created
                         top_embeds=100, # number of best multivariate embeddings used to calculate the multiview distance
                         seq_theta=c(0,0.1,1,3,8), # sequence of theta (nonlinearity) for parameter tuning
                         seq_alpha=c(0.1,0.3,0.5,0.7,0.9), # sequence of alpha (relationship of L2 and L1 penalty) for parameter tuning
                         seq_lambda=c(1.000, 0.178, 0.032, 0.006, 0.001), # sequence of lambda (penalization of regression coefficients) for parameter tuning
                         pred_tp=1){ # prediction horizon (here 2 time points = 8 days)
  
  # -----------------------------------------------------------------------
  
  ## Create output directories
  # This section ensures that the standard folder structure exists to save the output for each model run.
  dir.create(paste0(out_dir,"multivariate_embeddings"), showWarnings = FALSE)
  dir.create(paste0(out_dir,"distance_matrix"), showWarnings = FALSE)
  dir.create(paste0(out_dir,"parametrization"), showWarnings = FALSE)
  dir.create(paste0(out_dir,"coefficients"), showWarnings = FALSE)
  dir.create(paste0(out_dir,"performance"), showWarnings = FALSE)
  dir.create(paste0(out_dir,"interactions"), showWarnings = FALSE)
  
  # -----------------------------------------------------------------------
  
  ## List of all target variables in the data (e.g., taxa, sequences, etc.)
  target_vars <- names(select(input_dat,-1)) # the -1 removes the date or time column
  
  
  # -----------------------------------------------------------------------
  # 1. Scale and lag the data
  # -----------------------------------------------------------------------
  # Data must be standardized (e.g., zero mean and standard deviation of 1) so that Euclidean distances in state-space 
  # are not dominated by the ASVs with the highest mean abundance.
  
  for (i in 1:length(target_vars)) {
    for (tp in c(0,pred_tp)) {
      
      t <- pull(input_dat,1) # get time
      var <- target_vars[i] # get target variable name
      val <- pull(input_dat,var) # get target variable from the data frame
      val_scaled <- scale(val) # scale to zero mean and standard deviation of 1 by subtracting the mean and then dividing by the standard deviation
      
      ## Save attributes from scale function (this is used later to back-transform the rmse into a more meaningful number)
      if(var==target_i&tp==0){
        attributes_i <- attributes(val_scaled) # used to back-transform data for output (from the scale function) from this post https://stackoverflow.com/questions/10287545/backtransform-scale-for-plotting
        
      }
      
      ## Assemble the scaled data
      # Create vectors for the variables at current time (tp = 0) and at the lagged time (tp = -1)
      data_vector <- data.frame(time=t,
                                target=var,
                                val_scaled_lagged=lag(val_scaled,tp), # scale and lag the variable by tp
                                tp=tp)
      
      ## Save scaled output data
      if (i==1&tp==0) {out <- data_vector} else {out <- rbind(out, data_vector)}
    }
  }
  
  ## Some data wrangling and formatting on the output from the previous function
  data <- out %>%
    mutate(target_tp=paste0(target,"_tp_",tp)) %>% # add the lag information to the target
    select(-target,-tp) %>%
    spread(target_tp,val_scaled_lagged) %>%
    select(-time)
  
  # -----------------------------------------------------------------------
  ## 2. Multivariate embeddings and their predictive performance
  # -----------------------------------------------------------------------
  
  ## Create a list of embedding variables for the target
  # Define all possible taxa that could contribute to the embedding of the focal taxon target_i.
  embed_vars_target_i <- expand_grid(target_1=target_i,target_2=target_vars,tp=pred_tp) %>%
    filter(!target_1==target_2) %>%
    mutate(target_2_tp=paste0(target_2,"_tp_",tp)) %>%
    pull("target_2_tp")
  
  # Save the output from before into a data frame
  embed_vars <- data.frame(var_embed=embed_vars_target_i) # embedding variables for target i
  
  
  # -----------------------------------------------------------------------
  ## 2.1 Create random multivariate embeddings
  
  ## Calculate the maximum number of possible embeddings (this will determine which randomization procedure we use)
  possible_MV_embeddings <- choose(nrow(embed_vars),E-1) # number of maximum possible embedding combinations. We use E_best-1 because we include an autoregressive term in the MDR s-map models
  
  ## Check if we have enough possible embeddings to create the number of MV embeddings specified, if not, we set it to the highest possible number and give a waring
  if (MV_embeddings > possible_MV_embeddings) {
    MV_embeddings <- possible_MV_embeddings
    
    print(paste0("Warning: The choosen number of MV embeddings is larger than the possible number of embeddings, thus we lower it to the number of possible embeddings which is ",possible_MV_embeddings,"."))  # give user a warning
    
  }
  
  ## If MV_embeddings is small, we use the combn function to generate random multivariate embeddings
  if(MV_embeddings<1000){
    
    embeddings = combn(embed_vars$var_embed, E-1, simplify = FALSE)
    
  } else {
    
    ## If MV_embeddings is large, we use this loop to sample n (=MV_embeddings) embeddings from all possible embeddings, if we would use combn it would take way too long since it has to compute all possible embeddings first and then subsample
    j <- 1
    
    while (!j==MV_embeddings+1){
      embed <- NULL # create empty list to store results
      embed=list(sort(sample(unique(embed_vars$var_embed),E-1,replace=FALSE))) # sort orders the variables inside the list alphabetically --> keep always the same structure, we don't care of the order of variables (unique probably would)
      if(j==1){embeddings <- embed
      j = j+1 } else {
        if(!embed%in% embeddings) { # prevent duplicates
          embeddings <- rbind(embeddings,embed)
          j=j+1
          print(j)
        }
      }
      if(i==MV_embeddings+1){rownames(embeddings) <- NULL} # remove the rownames (we don't need those)
    }
  }
  
  
  # -----------------------------------------------------------------------
  ## 2.2 Evaluate predictive performance of random mulivariate embeddings
  
  ## Create output data frame to store results form multivariate simplex projection
  output <- expand.grid(target=rep(target_i,length(embeddings)), # data frame to save the output
                        embedding=NA, # previously generate random embeddings
                        rho = NA, # predictive skill
                        rmse = NA) # predictive error
  
  ## Choose data points as library and used for prediction
  # Use the first 50% of the time series to train the model and the second 50% to test it.
  lib <- paste0("1 ",as.integer(nrow(data)/2)) # use 50% of data as library
  pred <- paste(as.integer(nrow(data)/2)+1,nrow(data)) # and 50% for prediction
  
  ## Data interpolation
  data_interp <- data %>%
    mutate_all(.funs=liner_interp) # interpolate data, since Simplex function is very sensitive to missing values
  
  ## Estimate predictive performance with simplex
  for (j in 1:length(embeddings)) { #
    
    autoreg <- paste0(target_i,"_tp_",pred_tp) # autoregressive term
    
    columns <- paste(c(autoreg,embeddings[[j]]),collapse=" ") # variables used for embedding in simplex
    
    # Run simplex to evaluate the forecast skill of this specific 'view'.
    simplex <- Simplex(dataFrame=data_interp,pred=pred,lib=lib,target=autoreg,columns=columns,noTime=TRUE,embedded=TRUE,Tp=pred_tp) # run simplex
    
    predicted <- simplex$Predictions* as.numeric(attributes_i[["scaled:scale"]]) + as.numeric(attributes_i[["scaled:center"]]) # back-transform predictions
    observed <- simplex$Observations* as.numeric(attributes_i[["scaled:scale"]]) + as.numeric(attributes_i[["scaled:center"]]) # back-transform observations
    
    output$rho[j] <- cor(observed,predicted,use="pairwise.complete.obs") # Pearson's correlation between predicted and observed values
    output$rmse[j] <- sqrt(mean((observed-predicted)^2, na.rm=T)) # rooted mean squared error
    output$embedding[j] <- list(embeddings[[j]]) # embedding used
    
    print(paste("Multivariate Embeddings: Estimated ",j," of ",MV_embeddings," forecast skills")) # let user know of progress
    
  }
  
  ## Change the structure of the list in the embedding variable
  output$embedding <- vapply(output$embedding, paste0, collapse = ", ", character(1L))
  
  ## Save output in folder "multivariate_embeddings"
  write.csv(output,paste0(out_dir,"multivariate_embeddings/multivariate_embeddings_target_",target_i,".csv"),row.names = F)
  
  
  # -----------------------------------------------------------------------
  ## 3. Estimate the multiview distance
  # -----------------------------------------------------------------------
  
  ## Find the top performing multiview embeddings
  # Select the top 100 views that yielded the highest forecasting performance.
  embeddings <- read.csv(paste0(out_dir,"multivariate_embeddings/multivariate_embeddings_target_",target_i,".csv")) %>% # load data (saved from previous section)
    mutate(rho=ifelse(rho<0,0,rho)) %>% # set negative rho's to zero
    mutate(inv_rmse=max(rmse)-rmse+min(rmse)) %>% # we have to inverse the value for rmse
    mutate(performance=sqrt(inv_rmse^2+rho^2)) %>% # calculate performance, where high rho and low rmse get the highest performance
    arrange(-performance) %>% # arrange according to performance, putting highest performing embeddings on top
    slice(1:top_embeds) %>% # keep only best embeddings
    mutate(embedding=strsplit(embedding, ", ", fixed = TRUE)) # some string formatting
  
  # Weight each of the top views proportional to its forecast skill.
  weights <- embeddings[,'performance']/sum(embeddings[,'performance']) # weights each multivariate embedding based on it's performance in forecasting
  
  dist.matrix <- matrix(0,nrow(data),nrow(data)) # create empty distance matrix for the target
  
  ## Computation of multiview distance as the weighted average of multivariate distances
  for(j in 1:nrow(embeddings)){
    # Calculate Euclidean distance between points for each view.
    dist <- as.matrix(select(data,paste0(target_i,"_tp_",pred_tp),embeddings$embedding[j][[1]])) # estimate the distance among points in each multivariate embedding
    # Combine individual distances into a weighted averaged distance matrix.
    dist.matrix <- dist.matrix+as.matrix(dist(dist))*weights[j] # calculates the weighted average among all multivariate distances (= multiview distance)
  }
  
  print("Multiview distance (MVD): Estimated distance among sampling points in SSR") # let user know of progress
  
  # Save multiview distance matrix in folder "distance_matrix"
  write.csv(dist.matrix,paste0(out_dir,"distance_matrix/distance_matrix_",target_i,".csv"),row.names=F)
  
  # -----------------------------------------------------------------------
  ## 4. Parametrize the MDR S-map
  # -----------------------------------------------------------------------
  
  ## Create data set with response (target) variable and lagged explanatory variables (according to pred_tp)
  ds_0 <- select(data,paste0(target_i,"_tp_0")) # take the response/target variable at t = 0
  ds_tp <- select(data,paste0(names(select(select(input_dat,-1),target_i,everything())),"_tp_",pred_tp)) # take the explanatory variables at t = -pred_tp (e.g. 1 or 2 time steps before)  
  ds <- bind_cols(ds_0,ds_tp) # assemble the data for the MVD S-map
  
  ## Read in the multiview distance matrix
  dist.matrix <- read.csv(paste0(out_dir,"distance_matrix/distance_matrix_",target_i,".csv")) # this is from the previous part
  
  ## Define penalty factors
  # We do not want to remove the autoregressive effect.
  penalty.factor <- rep(1,ncol(ds)-1)
  penalty.factor[1] <- 0 # zero-penalty for diagonal coefficients (always the first column, we made sure of that above), effect of one node on itself, those coefficients are not going to be penalized
  
  ## Grid of parameters to tune (here lambda, theta and alpha)
  param.grid <- data.frame(expand.grid(theta=seq_theta,alpha=seq_alpha,lambda=seq_lambda)) %>%
    mutate(target=target_i,
           cv_min_dat_perc=NA,
           rmse=NA,
           rho=NA,
           const_pred_rmse=NA,
           const_pred_rho=NA)
  
  ## Evaluate model performance for each parameter combination using leave-future-out cross-validation
  for(p in 1:nrow(param.grid)){
    
    ## extract parameters for p
    lambda <- param.grid$lambda[p]
    theta <- param.grid$theta[p]
    alpha <- param.grid$alpha[p]
    
    ## Create empty vector for the predicted values
    predicted <- rep(NA,as.numeric(nrow(ds)))
    
    ## Leave-future-out cross-validation throughout the time series
    # "Day-forward chaining" — increase the training set incrementally to predict the next point.
    for(j in round(0.1*nrow(ds),0):nrow(ds)){ # do this for each row in the data, j=data point
      
      ds_j <- ds[1:j,] # increase the data for each point (= day forward chaining)
      
      distance_j <-  dist.matrix[1:j,j] # distance between j and all other points
      
      if (any(is.na(ds_j[j,]))==TRUE|(all(is.na(distance_j[-j])))==TRUE) {
        predicted[j] <- NA # the result becomes NA, if data for point j contains NA
      }else{ # run the analysis
        
        # Prevent data leakage: exclude the current point and its immediate past from the training library.
        distance_j[(j-(pred_tp-1)):j] <- NA # set data point j and points before it equal to pred_tp to NA (we don't want to include j in model fitting and by excluding points before j equal to pred_tp we prevent data leakage)
        
        mean.distance_j <- mean(distance_j,na.rm=T) # mean distance between j and all other points
        
        # Points closer in the manifold have higher weight in the regression.
        weight_j <- exp(-theta*distance_j/mean.distance_j) # creates weights using an exponentially decay function
        
        ds_j.weighted_j <- sqrt(weight_j)*ds_j # weight data for point j
        ds_j.weighted_j <- ds_j.weighted_j[!apply(is.na(ds_j.weighted_j),1,any),] # excludes NA from the data frame
        ds_j.weighted_j  <- apply(ds_j.weighted_j,2,as.numeric) # convert everything to numeric
        
        if(length(unique(ds_j.weighted_j[,1]))>1){ # make sure that y is not constant
          
          ## Solving sparse regression with elastic-net regularization
          # alpha controls the mix between Lasso (L1) and Ridge (L2) penalties.
          fit0 <- glmnet(x=ds_j.weighted_j[,-1], y=ds_j.weighted_j[,1], alpha=alpha,lambda=lambda,family="gaussian",penalty.factor=penalty.factor) #da.j[,-1] this removes the variable y from dat.j
          
          ## Make the prediction
          pred.ds_j <- as.numeric(ds_j[j,-1]) # extract library (x values) to predict y at time j
          
          pred_j <- NULL # create empty value for prediction at point j
          
          pred_j <- predict(fit0,newx=matrix(pred.ds_j,nrow=1)) # predict point j based on previously generated model
          
          ## save prediction
          predicted[j] <- pred_j
        } else {
          predicted[j] <- NA # set prediction to NA if y is constant   
        }
      } # end of else
    } # end of j
    
    ## Assemble predictions and back-transform data
    cv.predictions <- data.frame(predicted=predicted, # predictions, first 10% points should be NA because of the day-forward-chaining cross validation
                                 observed=ds[,1], # observations
                                 const_pred=ds[,2]) %>% # persistence model (assumes no change from one time point to the next)
      mutate(predicted=predicted* as.numeric(attributes_i[["scaled:scale"]]) + as.numeric(attributes_i[["scaled:center"]]), # back-transform predicted values
             observed=observed* as.numeric(attributes_i[["scaled:scale"]]) + as.numeric(attributes_i[["scaled:center"]]), # back-transform observed values
             const_pred=const_pred* as.numeric(attributes_i[["scaled:scale"]]) + as.numeric(attributes_i[["scaled:center"]])) #  # back-transform persistence model values
    
    cv.predictions_no_NA <- drop_na(cv.predictions) # drop NAs for prediction performance evaluation
    
    param.grid$cv_min_dat_perc[p]=min(which(!is.na(cv.predictions$predicted)))/nrow(ds)*100 # minimum amount of data used during cross-validation
    param.grid$rmse[p]=sqrt(mean((cv.predictions_no_NA$observed-cv.predictions_no_NA$predicted)^2, na.rm=T)) # calculate prediction error as rooted mean square error, where error is the difference between observed and predicted values
    param.grid$rho[p]=cor(x=cv.predictions_no_NA$predicted,y=cv.predictions_no_NA$observed,use='pairwise.complete.obs') # calculate prediction skill as Pearson's correlation between observations and predictions
    param.grid$const_pred_rmse[p]=sqrt(mean((cv.predictions_no_NA$observed-cv.predictions_no_NA$const_pred)^2, na.rm=T)) # calculate the prediction error from the persistence model
    param.grid$const_pred_rho[p]=cor(x=cv.predictions_no_NA$const_pred,y=cv.predictions_no_NA$observed,use='pairwise.complete.obs') # calculate the prediction skill from the persistence model
    
    print(paste("MVD S-map parametriztation: Done ",p," of ",nrow(param.grid)," parameter combinations."))  # let user know of progress
    
  } # end of p, parameter tuning
  
  ## Save results from parametrization in folder "parametrization"
  write.csv(param.grid,paste0(out_dir,"parametrization/smap_parametrization_",target_i,".csv"),row.names = F)
  
  
  # -----------------------------------------------------------------------
  ## 5. Estimate the MDR S-map coefficients
  # -----------------------------------------------------------------------
  
  ## Results from cross validation (tuning lambda, theta and alpha), from the previous part
  parameters_target_i <- read.csv(paste0(out_dir,"parametrization/smap_parametrization_",target_i,".csv")) 
  
  ## Keep only the best parameter combinations
  # Define an "Ensemble" of the top 5% performing parameter configurations.
  param.grid <- parameters_target_i %>%
    mutate(rho=ifelse(rho<0,0,rho)) %>% # set negative rho's to zero
    mutate(inv_rmse=max(rmse)-rmse+min(rmse)) %>% # we have to inverse the value for rmse to calculate performance
    mutate(performance_cv=sqrt(inv_rmse^2+rho^2)) %>% # calculate performance as a combination of skill and error
    filter(rho>const_pred_rho&rmse<const_pred_rmse) %>% # make sure that model performs better than the persistence model
    filter(performance_cv>=quantile(performance_cv,.95)) %>% # only keep the best performing models (95% quantile)
    arrange(-performance_cv) %>% # arrange based on performance, where highest performing model comes first
    rename(rho_cv=rho, # rename some of the variables
           rmse_cv=rmse,
           inv_rmse_cv=inv_rmse) %>%
    mutate(rmse_fit=NA,
           rho_fit=NA)
  
  n_models_realized <- as.numeric(nrow(param.grid)) # number of models (parameter combinations) to estimate coefficients
  
  if (n_models_realized >0){
    
    ## Fit MDR S-map for each parameter combination
    # This loop extracts the Jacobian coefficients for every point in the dataset.
    for(p in 1:nrow(param.grid)){ 
      
      ## extract parameters for p
      lambda <- param.grid$lambda[p]
      theta <- param.grid$theta[p]
      alpha <- param.grid$alpha[p]
      
      ## Create empty vector for the predicted values (we want to evaluate how well our final MDR S-map performs)
      predicted <- rep(NA,as.numeric(nrow(ds)))
      
      ## Create empty vector for model coefficients (interaction coefficients)
      jcoefficients <- matrix(NA,nrow(ds),ncol(ds))
      colnames(jcoefficients) <- c('intercept',colnames(ds_tp)) # set explanatory variables as column names
      
      ## Fit MDR S-map using the whole data
      for(j in 1:nrow(ds)){ # do this for each row in the data, j=data point
        
        distance_j <-  dist.matrix[,j] # distance between j and all other points
        
        if (any(is.na(ds[j,]))==TRUE|(all(is.na(distance_j[-j])))==TRUE) {
          predicted[j] <- NA # the predicted value becomes NA, if data for j contains NA
          jcoefficients[j,] <- NA # the coefficient values become NA, if data for j contains NA
        }else{ # run the analysis
          
          distance_j[(j-(pred_tp-1)):j] <- NA # Prevent leakage
          
          mean.distance_j <- mean(distance_j,na.rm=T) # mean distance between j and all other points
          
          weight_j <- exp(-theta*distance_j/mean.distance_j) # creates weights using an exponentially decay function
          
          ds_weighted_j <- sqrt(weight_j)*ds # weight data for point j
          ds_weighted_j <- ds_weighted_j[!apply(is.na(ds_weighted_j),1,any),] # excludes NA from the data frame
          ds_weighted_j  <- apply(ds_weighted_j,2,as.numeric) # convert everything to numeric
          
          ## Solving sparse regression with elastic-net regularization
          fit0 <- glmnet(x=ds_weighted_j[,-1], y=ds_weighted_j[,1], alpha=alpha,lambda=lambda,family="gaussian",penalty.factor=penalty.factor) #da.j[,-1] this removes the variable y from dat.j
          
          ## Coefficients are extracted from the fitting of elastic-net
          # The beta values here are the interaction strengths.
          jcoefficients[j,1] <- as.numeric(fit0$a0) # intercept
          jcoefficients[j,2:(ncol(ds))] <- as.numeric(fit0$beta) # interaction coefficients
          
          ## Make the prediction
          pred.ds_j <- as.numeric(ds[j,-1]) # extract library (x values) to predict y at time j
          
          pred_j <- NULL # create empty value for prediction at point j
          
          pred_j <- predict(fit0,newx=matrix(pred.ds_j,nrow=1)) # predict point j based on previously generated model
          
          ## save prediction
          predicted[j] <- pred_j
          
        } # end of else
      }# end of j
      
      ## Assemble predictions and back-transform data
      fit.predictions <- data.frame(predicted=predicted, # predictions
                                    observed=ds[,1], # observations
                                    const_pred=ds[,2]) %>% # persistence model (assumes no change form one time point to the next)
        mutate(predicted=predicted* as.numeric(attributes_i[["scaled:scale"]]) + as.numeric(attributes_i[["scaled:center"]]), # back-transform predicted values
               observed=observed* as.numeric(attributes_i[["scaled:scale"]]) + as.numeric(attributes_i[["scaled:center"]]), # back transform observed values
               const_pred=const_pred* as.numeric(attributes_i[["scaled:scale"]]) + as.numeric(attributes_i[["scaled:center"]])) # back-transform persistence model values
      
      fit.predictions_no_NA <- drop_na(fit.predictions) # drop NAs for forecast performance skill evaluation
      
      ## Evaluate the skills (rmse & rho) for a given parameter combination when using the whole data for model fitting
      param.grid$rmse_fit[p]=sqrt(mean((fit.predictions_no_NA$observed-fit.predictions_no_NA$predicted)^2, na.rm=T)) # calculate prediction error as rooted mean square error, where error is the difference between observed and predicted values
      param.grid$rho_fit[p]=cor(x=fit.predictions_no_NA$predicted,y=fit.predictions_no_NA$observed,use='pairwise.complete.obs') # calculate prediction skill as Pearson's correlation between observations and predictions
      
      ## Save coefficients in folder "coefficients"
      write.csv(jcoefficients,paste0(out_dir,"coefficients/coefficients_target_",target_i,"_theta_",theta,"_alpha_",alpha,"_lambda_",lambda,".csv"),row.names = F)
      
      print(paste("MVD S-map coefficients: Estimated ",p," of ",nrow(param.grid)," interaction coefficients")) # let user know of progress
      
    } # end of p, best parameters
    
    ## Add performance when using the full data set for MDR S-map (fit) to output data
    param.grid <- param.grid   %>%
      mutate(rho_fit=ifelse(rho_fit<0,0,rho_fit)) %>% # set negative rho's to zero
      mutate(inv_rmse_fit=max(rmse_fit)-rmse_fit+min(rmse_fit)) %>% # we have to inverse the value for rmse
      mutate(performance_fit=sqrt(inv_rmse_fit^2+(rho_fit)^2))
    
    ## Save model performances in folder "performance"
    write.csv(param.grid,paste0(out_dir,"performance/smap_performance_",target_i,".csv"),row.names = F)
    
    
    # -----------------------------------------------------------------------
    # 6. Estimate interactions (average MDR S-map coefficients)
    # -----------------------------------------------------------------------
    # Loop through the results of all top models and calculate their weighted mean.
    
    for (p in 1:nrow(param.grid)) { # do this for the best parameter combinations from before
      
      theta <- param.grid$theta[p]
      alpha <- param.grid$alpha[p]
      lambda <- param.grid$lambda[p]
      
      ## Read in all coefficient files
      file_p <- read.csv(paste0(out_dir,"coefficients/coefficients_target_",target_i,"_theta_",theta,"_alpha_",alpha,"_lambda_",lambda,".csv")) %>% # read in coefficients files
        slice(-c(1:pred_tp)) %>% # remove first rows of data frame according to pred_tp (since we don't have predictions here and we want to match it up with the original time, because we are using lags, we need to move everything up according to pred_tp)
        add_case(intercept=rep(NA,pred_tp)) %>% # we need to add empty rows at the end of the data frame and take advantage here that intercept will be in every file_p, independent of the other variables (maybe there is a better way to do this)
        mutate(time=unique(pull(input_dat,1))) %>% # get the time index from the input data frame and add it to file_p
        gather(key="target_2",value="coefficient",-time) %>% # reshape the data frame
        filter(!target_2=="intercept") %>% # remove the intercept (it has no biological meaning, see discussion in Chang et al. 2021)
        mutate(target_1=paste(target_i), # add the target as a variable
               performance=as.numeric(paste(param.grid$performance_fit[p])), # add performance (we will use this as a weight)
               link=ifelse(!coefficient==0,1,0)) # record if there is a link (coefficient different from zero)
      
      if(p==1){coefficients_i <- file_p}else{coefficients_i <- rbind(coefficients_i,file_p)}
      
    }
    
    ## Calculate weighted average of coefficients based on performance
    coefficients_weighted_i <- coefficients_i  %>%
      filter(!is.na(link)) %>% # remove NA's
      mutate(weight = performance) %>% # this column is going to be used for the weights, giving the highest performing model the highest weight
      group_by(time,target_1,target_2) %>%
      summarise(coefficient_w_mean=sum(coefficient * weight)/sum(weight), # weighted mean
                coefficient_w_sd=sum((weight/sum(weight)) * (coefficient - coefficient_w_mean)^2), # weighted standard deviation
                link_perc=sum(link)/n_models_realized*100, .groups = "drop") # percentage of models where a link between those two taxa was recoreded
    
    ## Save interactions in the folder "interactions"
    write.csv(coefficients_weighted_i,paste0(out_dir,"interactions/interaction_coefficients_target_",target_i,".csv"),row.names=F)
    
    print(paste0("Estimated average interaction coefficients for target ",target_i," and saved them in ",out_dir,"coefficients/")) # let user know of progress
    
  } else {print(paste0("No MVD S-map model for ",target_i," outperformed the constant predictor. No interaction coefficients have been estimated."))} # give a warning if no interaction coefficients could be estimated due to poor model performance
}