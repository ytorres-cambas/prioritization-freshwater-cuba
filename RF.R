# Run Random Forest models for species distribution modeling.  
# This script is adapted from:  
# Valavi, R., Guillera-Arroita, G., Lahoz-Monfort, J., & Elith, J. (2021).  
# Predictive performance of presence-only species distribution models:  
# a benchmark study with reproducible code.  
# Ecological Monographs, 92, e01486.  
# https://doi.org/10.1002/ecm.1486  
#  
# The script performs the following steps:  
# - Loads necessary libraries and functions.  
# - Imports species occurrence data and environmental predictors.  
# - Trains Random Forest models using cross-validation.  
# - Evaluates model performance using AUC, TSS, and other metrics.  
# - Generates and saves spatial predictions of species distribution.  


# Load necessary libraries
if (!require(sf)) install.packages("sf")
if (!require(randomForest)) install.packages("randomForest")
if (!require(scales)) install.packages("scales")
if (!require(raster)) install.packages("raster")
if (!require(dplyr)) install.packages("dplyr")

# Import functions
source("scr_functions.R")

# Get file names for input occurrences
data_path <- "./data/sdm/input/occurr"
files_names <- list.files(path = data_path, full.names = FALSE)

# Load subcatchment map
data_map <- "./data/subc_habit_clasif.gpkg"
map_subcatch <- st_read(data_map) %>% st_set_crs(3795)

# Load prediction sites
predict_sites <- read.csv("./data/sdm/input/pred_sites/pred_sites.csv")

# Loop over species files
for (sp in seq_along(files_names)) {
  time_sp <- system.time({
    tryCatch({
      cat("Modelling", files_names[sp], "\n")
      
      # Import dataset
      occurr_file <- file.path(data_path, files_names[sp], "occurr.csv")
      dataset <- read.csv(occurr_file)
      
      # Define output path
      output_path <- file.path("./data/sdm/output/indep_models", files_names[sp])
      if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
      
      # Initialize result matrices
      thres <- matrix(nrow = 1, ncol = 10)
      perform <- matrix(nrow = 2, ncol = 10)
      
      PredProb <- matrix(nrow = nrow(predict_sites), ncol = 11)
      PredProb[, 11] <- predict_sites$stream_ID
      colnames(PredProb) <- c(paste0("prob_split", 1:10), "stream")
      
      predPresAbs <- matrix(nrow = nrow(predict_sites), ncol = 11)
      predPresAbs[, 11] <- predict_sites$stream_ID
      colnames(predPresAbs) <- c(paste0("pres_abs_split", 1:10), "stream")
      
      rf_var_importance <- matrix(nrow = ncol(dataset[, 2:11]), ncol = 11)
      rownames(rf_var_importance) <- colnames(dataset[, 2:11])
      colnames(rf_var_importance)[11] <- 'mean'
      
      # Load calibration and evaluation datasets
      CalDf <- read.csv(file.path(data_path, files_names[sp], "calib.csv"))
      EvalDf <- read.csv(file.path(data_path, files_names[sp], "eval.csv"))
      
      set.seed(1234)
      
      # Run 10 model replicates
      for (r in 1:10) {
        CalibData <- dataset[dataset$stream %in% CalDf[, r], ]
        EvaData <- dataset[dataset$stream %in% EvalDf[, r], ]
        
        CalibData$pres_abs <- as.factor(CalibData$pres_abs)
        CalibData <- CalibData[, c(2:14)]
        
        prNum <- sum(CalibData$pres_abs == "1")
        bgNum <- sum(CalibData$pres_abs == "0")
        samsize <- c("0" = prNum, "1" = prNum)
        
        mod_rf <- try(randomForest(pres_abs ~ ., data = CalibData,
                                   ntree = 1000, sampsize = samsize,
                                   replace = TRUE, importance = TRUE), silent = TRUE)
        
        if (inherits(mod_rf, "try-error")) {
          cat("Error for species", files_names[sp], "replicate", r, "\n")
          next
        }
        
        # Apply bias correction and predict
        EvaData_bias_cor <- EvaData %>% mutate(z_d_town = 1, z_d_roads = 1) %>% select(PC1:PC10, z_d_town, z_d_roads)
        PredProbEval <- as.numeric(predict(mod_rf, EvaData_bias_cor, type = "prob")[, "1"])
        
        # Evaluate model
        idsPresabsProb <- data.frame(EvaData[, c("stream", "pres_abs")], PredProbEval)
        eval_res <- eval_mod(idsPresabsProb)
        
        TSS <- eval_res$performance[5, 2]
        AUC <- eval_res$performance[6, 2]
        th <- eval_res$threshold[1, 2]
        
        thres[, r] <- th
        perform[, r] <- c(TSS, AUC)
        
        if (AUC > 0.7 & TSS > 0.7) {
          pred_var_bias_cor <- predict_sites %>% mutate(z_d_town = 1, z_d_roads = 1) %>% select(PC1:PC10, z_d_town, z_d_roads)
          pred <- as.numeric(predict(mod_rf, pred_var_bias_cor, type = "prob")[, "1"])
          PredProb[, r] <- pred
          
          presabs <- ifelse(pred >= th, 1, 0)
          predPresAbs[, r] <- presabs
          
          # Save results
          write.csv(AUC, file.path(output_path, paste0("auc_rf_", r, ".csv")), row.names = FALSE)
          write.csv(TSS, file.path(output_path, paste0("tss_rf_", r, ".csv")), row.names = FALSE)
          write.csv(th, file.path(output_path, paste0("th_rf_", r, ".csv")), row.names = FALSE)
          write.csv(pred, file.path(output_path, paste0("pred_rf_", r, ".csv")), row.names = FALSE)
          write.csv(presabs, file.path(output_path, paste0("presabs_rf_", r, ".csv")), row.names = FALSE)
        }
      }
      
      # Compute mean probability and CV
      mean_pro <- rowMeans(PredProb[, 1:10], na.rm = TRUE)
      cv_mean_pro <- apply(PredProb[, 1:10], 1, function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE))
      th_mean <- mean(thres)
      mean_pres_abs <- ifelse(mean_pro >= th_mean, 1, 0)
      
      prediction_df <- data.frame(PredProb, mean_pro, cv_mean_pro, predPresAbs[, 1:10], mean_pres_abs)
      prediction_map <- map_subcatch %>% inner_join(prediction_df, by = c("id" = "stream"))
      
      st_write(prediction_map, file.path(output_path, paste0(files_names[sp], "_rf_predictions.gpkg")), driver = "GPKG", append = FALSE, quiet = TRUE)
    }, error = function(e) {})
  })
}






