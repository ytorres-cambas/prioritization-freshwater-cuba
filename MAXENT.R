# Run MaxEnt models for species distribution modeling.  
# This script is adapted from:  
# Valavi, R., Guillera-Arroita, G., Lahoz-Monfort, J., & Elith, J. (2021).  
# Predictive performance of presence-only species distribution models:  
# a benchmark study with reproducible code.  
# Ecological Monographs, 92, e01486.  
# https://doi.org/10.1002/ecm.1486  
#  
# The script performs the following steps:  
# - Loads necessary libraries and functions.  
# - Reads species presence data and environmental predictors.  
# - Trains MaxEnt models using cross-validation.  
# - Evaluates model performance using AUC and other metrics.  
# - Generates and saves spatial predictions of species distribution.  


# Loading libraries
if (!require(dplyr)) install.packages("dplyr")
if (!require(dismo)) install.packages("dismo")
if (!require(biomod2)) install.packages("biomod2")
if (!require(sf)) install.packages("sf")

Sys.setenv(JAVA_HOME = 'C:/Program Files/Java/jre1.8.0_341')

# Import functions
source("scr_functions.R")

# Get the file names for saving results
files_names <- list.files(path = "./data/sdm/input/occurr", full.names = FALSE)

# Map with subcatchments
map_subcatch <- st_read("./data/subc_habit_clasif.gpkg") %>%
  st_set_crs(3795)

# Import dataset with prediction sites
predict_sites <- read.csv("./data/sdm/input/pred_sites/pred_sites.csv")

for (sp in seq_along(files_names)) {
  tryCatch({
    print(paste0("Modelling ", files_names[sp]))
    
    # Import dataset with presence-absences and predictors
    dataset <- read.csv(paste0("./data/sdm/input/occurr/", files_names[sp], "/occurr.csv"))
    
    # Output path
    output_path <- paste0("./data/sdm/output/indep_models/", files_names[sp])
    if (!dir.exists(output_path)) dir.create(output_path)
    
    # Initialize matrices to save results
    thres <- matrix(nrow = 1, ncol = 10)
    perform <- matrix(nrow = 2, ncol = 10)
    PredProb <- matrix(nrow = nrow(predict_sites), ncol = 11)
    PredProb[, 11] <- predict_sites$stream_ID 
    colnames(PredProb) <- c(paste0("split", 1:10), "stream")
    
    predPresAbs <- matrix(nrow = nrow(predict_sites), ncol = 11)
    predPresAbs[, 11] <- predict_sites$stream_ID 
    colnames(predPresAbs) <- c(paste0("split", 1:10), "stream")
    
    # Import data frames with calibration and evaluation datasets
    CalDf <- read.csv(paste0("./data/sdm/input/occurr/", files_names[sp], "/calib.csv"))
    EvalDf <- read.csv(paste0("./data/sdm/input/occurr/", files_names[sp], "/eval.csv"))
    
    for (r in 1:10) {
      tryCatch({
        # Split calibration and evaluation data
        CalibData <- dataset[dataset$stream %in% CalDf[, r], ]
        EvaData <- dataset[dataset$stream %in% EvalDf[, r], ]
        
        # Tune Maxent model parameters
        pred <- names(dataset[, 2:13])  # Predictor names
        occ <- "pres_abs"  # Response variable
        params <- maxent_param(data = CalibData, pred = pred, occ = occ, k = 5)
        
        # Run model with tuned parameters
        CalibMaxent <- maxent(x = CalibData[, pred], p = CalibData[, occ], args = params)
        
        # Correct for observer bias in evaluation dataset
        EvaData_bias_cor <- EvaData %>%
          mutate(z_d_town = 1, z_d_roads = 1) %>% 
          select(PC1:PC10, z_d_town, z_d_roads)
        
        # Predict on evaluation dataset
        PredProbEval <- predict(CalibMaxent, EvaData_bias_cor[, pred], args = "outputformat=cloglog")
        
        # Evaluation metrics
        idsPresabsProb <- data.frame(EvaData[, c("stream", "pres_abs")], PredProbEval)
        TSS <- eval_mod(idsPresabsProb)$performance[5, 2]
        AUC <- eval_mod(idsPresabsProb)$performance[6, 2]
        th <- eval_mod(idsPresabsProb)$threshold[1, 2]
        thres[, r] <- th
        perform[, r] <- c(TSS, AUC)
        
        if (AUC > 0.7 & TSS > 0.7) {
          # Correct for observer bias in prediction sites
          pred_var_bias_cor <- predict_sites %>%
            mutate(z_d_town = 1, z_d_roads = 1) %>% 
            select(PC1:PC10, z_d_town, z_d_roads)
          
          # Predict on complete dataset
          pred <- predict(CalibMaxent, pred_var_bias_cor[, pred], args = "outputformat=cloglog")
          PredProb[, r] <- pred
          
          # Binary transformation
          presabs <- ifelse(pred >= th, 1, 0)
          predPresAbs[, r] <- presabs
          
          # Save results
          write.csv(EvalDf[, r], paste0(output_path, "/eva_max_", r, ".csv"), row.names = FALSE)
          write.csv(AUC, paste0(output_path, "/auc_max_", r, ".csv"), row.names = FALSE)
          write.csv(TSS, paste0(output_path, "/tss_max_", r, ".csv"), row.names = FALSE)
          write.csv(th, paste0(output_path, "/th_max_", r, ".csv"), row.names = FALSE)
          write.csv(pred, paste0(output_path, "/pred_max_", r, ".csv"), row.names = FALSE)
          write.csv(presabs, paste0(output_path, "/presabs_max_", r, ".csv"), row.names = FALSE)
        }
      }, error = function(e) {})  
    }
    
    # Compute mean probability and coefficient of variation
    mean_pro <- rowMeans(PredProb[, 1:10], na.rm = TRUE)
    cv <- function(x) sqrt(var(x, na.rm = TRUE) / length(x)) / mean(x, na.rm = TRUE)
    cv_mean_pro <- apply(PredProb[, 1:10], 1, cv)
    th_mean <- mean(thres, na.rm = TRUE)
    mean_pres_abs <- ifelse(mean_pro >= th_mean, 1, 0)
    
    # Create prediction dataframe
    prediction_df <- data.frame(PredProb, mean_pro, cv_mean_pro, predPresAbs[, 1:10], mean_pres_abs)
    prediction_map <- map_subcatch %>% inner_join(prediction_df, by = c("id" = "stream"))
    
    # Save prediction map
    st_write(prediction_map, paste0(output_path, "/", files_names[sp], "_max_predictions.gpkg"),
             driver = "GPKG", append = FALSE, quiet = TRUE)
    
  }, error = function(e) {})
}














