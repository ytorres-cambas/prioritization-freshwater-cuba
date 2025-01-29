# Run Boosted Regression Trees (BRT) models for species distribution modeling.
# This script is adapted from:
# Valavi, R., Guillera-Arroita, G., Lahoz-Monfort, J., & Elith, J. (2021). 
# Predictive performance of presence-only species distribution models: 
# a benchmark study with reproducible code. 
# Ecological Monographs, 92, e01486. 
# https://doi.org/10.1002/ecm.1486
#
# The script performs the following steps:
# - Loads required spatial and statistical libraries.
# - Reads species occurrence data, predictor variables, and calibration datasets.
# - Iteratively fits Boosted Regression Trees (BRT) models using cross-validation.
# - Evaluates model performance using AUC and TSS metrics.
# - Generates predictions across multiple model runs.
# - Saves model outputs and predictions as CSV and spatial files.


# Load required libraries
if (!require(dplyr)) install.packages("dplyr")
if (!require(dismo)) install.packages("dismo")
if (!require(sf)) install.packages("sf")

# Import functions
source("scr_functions.R")

# Get file names (used to save results)
files_names <- list.files(path = "./data/sdm/input/occurr", full.names = FALSE)

# Load map with subcatchments
map_subcatch <- st_read("./data/subc_habit_clasif.gpkg") %>%
  st_set_crs(3795)

# Import dataset with prediction sites
predict_sites <- read.csv("./data/sdm/input/pred_sites/pred_sites.csv")  

# Loop through species datasets
for (sp in seq_along(files_names)) {
  time_sp <- system.time({
    tryCatch({
      print(paste0("Modelling ", files_names[sp]))
      
      # Import dataset with presence-absences and predictors
      dataset <- read.csv(paste0("./data/sdm/input/occurr/", files_names[sp], "/occurr.csv"))
      
      # Define output path
      output_path <- paste0("./data/sdm/output/indep_models/", files_names[sp])
      if (!dir.exists(output_path)) dir.create(output_path)
      
      # Initialize matrices for results
      thres <- matrix(nrow = 1, ncol = 10)
      perform <- matrix(nrow = 2, ncol = 10)
      PredProb <- matrix(nrow = nrow(predict_sites), ncol = 11, 
                         dimnames = list(NULL, c(paste0("split", 1:10), "stream")))
      PredProb[, 11] <- predict_sites$stream_ID 
      
      predPresAbs <- matrix(nrow = nrow(predict_sites), ncol = 11, 
                            dimnames = list(NULL, c(paste0("split", 1:10), "stream")))
      predPresAbs[, 11] <- predict_sites$stream_ID 
      
      # Import calibration and evaluation datasets
      CalDf <- read.csv(paste0("./data/sdm/input/occurr/", files_names[sp], "/calib.csv"))
      EvalDf <- read.csv(paste0("./data/sdm/input/occurr/", files_names[sp], "/eval.csv"))
      
      set.seed(1234)
      
      # Run 10 replicates of the model
      for (r in 1:10) {
        CalIds <- CalDf[, r]
        CalibData <- dataset[dataset$stream %in% CalIds, ]
        EvaIds <- EvalDf[, r]
        EvaData <- dataset[dataset$stream %in% EvaIds, ]
        
        # Calculate case weights
        prNum <- sum(CalibData$pres_abs == 1)
        bgNum <- sum(CalibData$pres_abs == 0)
        wt <- ifelse(CalibData$pres_abs == 1, 1, prNum / bgNum)
        
        # Initialize BRT parameters
        ntrees <- 50
        tcomplexity <- ifelse(prNum < 50, 1, 5)
        lrate <- 0.001
        mod_brt <- NULL
        
        # Model fitting loop
        for (m in 1:40) {
          if (is.null(mod_brt)) {
            if (m == 11) lrate <- 0.0001
            if (m == 21) ntrees <- 25
            if (m == 31) tcomplexity <- ifelse(prNum < 50, 1, 3)
            
            mod_brt <- tryCatch({
              gbm.step(
                data = CalibData,
                gbm.x = 2:13,
                gbm.y = 14, 
                family = "bernoulli",
                site.weights = wt,
                tree.complexity = tcomplexity,
                learning.rate = lrate,
                n.trees = ntrees,
                n.folds = 5,
                max.trees = 10000,
                plot.main = FALSE
              )
            }, error = function(e) NULL)
          }
        }
        
        if (is.null(mod_brt)) next
        
        # Bias correction for evaluation data
        EvaData_bias_cor <- EvaData %>%
          mutate(z_d_town = 1, z_d_roads = 1) %>% 
          dplyr::select(PC1:PC10, z_d_town, z_d_roads)
        
        # Predict on evaluation dataset
        PredProbEval <- predict(mod_brt, EvaData_bias_cor, 
                                n.trees = mod_brt$gbm.call$best.trees, 
                                type = "response")
        
        # Compute evaluation metrics
        idsPresabsProb <- data.frame(EvaData[, c("stream", "pres_abs")], PredProbEval)
        TSS <- eval_mod(idsPresabsProb)$performance[5, 2]
        AUC <- eval_mod(idsPresabsProb)$performance[6, 2]
        th <- eval_mod(idsPresabsProb)$threshold[1, 2]
        thres[, r] <- th
        perform[, r] <- c(TSS, AUC)
        
        if (AUC > 0.7 & TSS > 0.7) {
          # Bias correction for prediction sites
          pred_var_bias_cor <- predict_sites %>%
            mutate(z_d_town = 1, z_d_roads = 1) %>% 
            dplyr::select(PC1:PC10, z_d_town, z_d_roads)
          
          # Predict on full dataset
          pred <- predict(mod_brt, pred_var_bias_cor,
                          n.trees = mod_brt$gbm.call$best.trees,
                          type = "response")
          PredProb[, r] <- pred
          predPresAbs[, r] <- ifelse(pred >= th, 1, 0)
          
          # Save results
          write.csv(AUC, paste0(output_path, "/auc_brt_", r, ".csv"), row.names = FALSE)
          write.csv(TSS, paste0(output_path, "/tss_brt_", r, ".csv"), row.names = FALSE)
          write.csv(th, paste0(output_path, "/th_brt_", r, ".csv"), row.names = FALSE)
        }
      }
      
      # Generate final prediction map
      mean_pro <- rowMeans(PredProb[, 1:10], na.rm = TRUE)
      th_mean <- mean(thres)
      mean_pres_abs <- ifelse(mean_pro >= th_mean, 1, 0)
      
      prediction_df <- data.frame(PredProb, mean_pro, mean_pres_abs)
      prediction_map <- map_subcatch %>%
        inner_join(prediction_df, by = c("id" = "stream"))
      
      st_write(prediction_map, paste0(output_path, "/", files_names[sp], "_brt_predictions.gpkg"),
               driver = "GPKG", append = FALSE, quiet = TRUE)
    }, error = function(e) {})
  })
}











    

    
    
    



