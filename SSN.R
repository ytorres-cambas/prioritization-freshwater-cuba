### Spatial Stream Network Models

# This script implements Spatial Stream Network (SSN) models for species distribution modeling.
# It imports necessary packages, loads and processes spatial data, calibrates models, evaluates performance,
# and generates prediction maps. The approach considers different correlation models for different taxa
# (fish, plants, amphibians, and invertebrates) and incorporates observer bias correction.
# Ver Hoef, J. M., Peterson, E. E., Hooten, M. B., Hanks, E. M., & Fortin, M. J. (2018). Spatial autoregressive models for statistical inference from ecological data. *Ecological Monographs, 88*(1), 36-59.


# Load necessary packages
if (!require(dplyr)) install.packages("dplyr")
if (!require(SSN)) install.packages("SSN")
if (!require(sf)) install.packages("sf")

# Import custom functions
source("scr_functions.R")

# Load spatial data (subcatchment map)
map_subcatch <- st_read("./data/subc_habit_clasif.gpkg") %>%
  st_set_crs(3795)

# Define path to SSN objects
list_ssn_dirs <- list.dirs(path = "./data/sdm/input/ssn/ssn_obj", recursive = FALSE)

# Get file names
file_names <- list.files(path = "./data/sdm/input/ssn/ssn_obj", full.names = FALSE)

# Loop over SSN directories
for (sp in 1:length(list_ssn_dirs)) {
  time_sp <- system.time({
    tryCatch({
      
      # 1. Import SSN object
      print(paste("Start", files_names[sp], sep = " "))
      print("Importing SSN object...")
      
      ssn_dir <- list_ssn_dirs[sp]
      ssn_sp_i <- importSSN(ssn_dir, predpts = "pred_sites")
      
      # 2. Generating an additive function value (necessary for Tail-up models)
      print("Generating an additive function value...")
      ssn_sp_i <- additive.function(ssn_sp_i, "H2OArea", "computed.afv")
      
      # 3. Calculate distance matrix
      print("Calculating distance matrix...")
      createDistMat(ssn_sp_i, predpts = "pred_sites", o.write = TRUE)
      
      # Model parameter calibration
      # fishes: tail-down models only; plants: tail-up and Euclidean models only;
      # reptiles, amphibians, and invertebrates: tail-down, Euclidean;
      
      obs_sites <- getSSNdata.frame(ssn_sp_i, Name = "Obs")
      load("./data/sdm/input/ssn/CorMdls.RData")
      
      if (obs_sites$clas[1] == "Actinopterygii") {
        # only tail-down
        CorMdls <- CorMdls[(171:175)]
      } else {
        # Euclidean and tail-up
        if (obs_sites$king[1] == "Plantae") {
          CorMdls <- CorMdls[(c(126:145, 166:170, 176:179))]
        } else {
          # tail-down, Euclidean
          CorMdls <- CorMdls[(c(146:165, 171:179))]
        }
      }
      
      glmssnParam <- SSNparam(ssn_sp_i, CorMdls)
      
      # Output path
      output_path <- paste0("./data/sdm/output/indep_models/", files_names[sp])
      ifelse(!dir.exists(output_path), dir.create(output_path), FALSE)
      
      # Matrix to save results
      thres <- matrix(nrow = 1, ncol = 10)
      perform <- matrix(nrow = 2, ncol = 10)
      predict_sites <- getSSNdata.frame(ssn_sp_i, Name = "pred_sites")
      PredProb <- matrix(nrow = nrow(predict_sites), ncol = 11)
      PredProb[, 11] <- predict_sites$stream
      colnames(PredProb) <- c(paste0("split", 1:10), "stream")
      
      predPresAbs <- matrix(nrow = nrow(predict_sites), ncol = 11)
      predPresAbs[, 11] <- predict_sites$stream
      colnames(predPresAbs) <- c(paste0("split", 1:10), "stream")
      
      # Import data frames with sub-catchment IDs for calibration and evaluation
      CalDf <- read.csv(paste0("./data/sdm/input/occurr/", files_names[sp], "/calib.csv"))
      EvalDf <- read.csv(paste0("./data/sdm/input/occurr/", files_names[sp], "/eval.csv"))
      
      # Loop over 10 repetitions for model fitting and evaluation
      for (r in 1:10) {
        tryCatch({
          # Data frame with observation points
          ssn_test_dataDF <- getSSNdata.frame(ssn_sp_i)
          
          # IDs of observations selected for evaluation
          EvaIds <- EvalDf[, r]
          
          # Insert NAs in the response variable column for the evaluation observations
          ssn_test_dataDF[ssn_test_dataDF$stream %in% EvaIds, "pres_abs"] <- NA
          
          # Put the data frame with NAs into the SSN object
          ssn_calib_data <- putSSNdata.frame(ssn_test_dataDF, ssn_sp_i)
          
          # Fit the model with the calibration data and parameters
          options(na.action = "na.omit")
          CalibSsn <- try(glmssn(formula = glmssnParam$formula,
                                 family = "binomial",
                                 ssn_calib_data,
                                 CorModels = glmssnParam$model,
                                 addfunccol = "computed.afv",
                                 control = list(trunc.pseudo = 100)))
          
          # Correct bias covariate to a constant value across all sub-catchments
          eval_temp <- getSSNdata.frame(CalibSsn$ssn.object, "_MissingObs_") %>%
            mutate(z_d_town = 1, z_d_roads = 1)
          
          CalibSsn_cor_eval <- putSSNdata.frame(eval_temp, CalibSsn, "_MissingObs_")
          
          # Predict points used for evaluation with the model trained on the calibration dataset
          PredEval <- predict.glmssn(CalibSsn_cor_eval, "_MissingObs_")
          
          # Extract predictions and transform from logit to probabilities
          SSNProb <- function(p) {
            logit <- getPreds(p, pred.type = "pred")
            prob <- 1 / (1 + exp(-logit[, 2]))
            return(prob)
          }
          
          PredProbEval <- SSNProb(PredEval)
          
          # Extract observed presence-absences for evaluation
          obsPresAbs <- getSSNdata.frame(ssn_sp_i) %>%
            filter(stream %in% EvaIds) %>%
            dplyr::select(stream, pres_abs)
          
          # Evaluation metrics
          idsPresabsProb <- data.frame(obsPresAbs, PredProbEval)
          TSS <- eval_mod(idsPresabsProb)$performance[5, 2]
          AUC <- eval_mod(idsPresabsProb)$performance[6, 2]
          th <- eval_mod(idsPresabsProb)$threshold[1, 2]
          thres[, r] <- th
          perform[, r] <- c(TSS, AUC)
          
          if (AUC > 0.7 & TSS > 0.7) {
            # Correct bias covariate to a constant value across all sub-catchments
            temp_pred_df <- getSSNdata.frame(CalibSsn, "pred_sites") %>%
              mutate(z_d_town = 1, z_d_roads = 1)
            
            CalibSsn_cor_eval <- putSSNdata.frame(temp_pred_df, CalibSsn, "pred_sites")
            
            # Predict on the complete dataset
            PredAll <- predict.glmssn(CalibSsn_cor_eval, "pred_sites")
            
            # Extract probabilities
            PredProb[, r] <- SSNProb(PredAll)
            
            # Binary transformation
            presabs <- ifelse(PredProb[, r] >= th, 1, 0)
            predPresAbs[, r] <- presabs
            
            # Save results for models with AUC and TSS > 0.7
            write.csv(EvaIds, row.names = FALSE, paste0(output_path, "/eva_ssn_", r, ".csv"))
            write.csv(AUC, row.names = FALSE, paste0(output_path, "/auc_ssn_", r, ".csv"))
            write.csv(TSS, row.names = FALSE, paste0(output_path, "/tss_ssn_", r, ".csv"))
            write.csv(th, row.names = FALSE, paste0(output_path, "/th_ssn_", r, ".csv"))
            write.csv(PredProb[, r], row.names = FALSE, paste0(output_path, "/pred_ssn_", r, ".csv"))
            write.csv(presabs, row.names = FALSE, paste0(output_path, "/presabs_ssn_", r, ".csv"))
          }
          
        }, error = function(e) {})
      }
      
      # Calculate mean probability and CV
      mean_pro <- apply(PredProb[, 1:10], 1, mean, na.rm = TRUE)
      cv <- function(x) (sqrt(var(x, na.rm = TRUE) / length(x))) / mean(x, na.rm = TRUE)
      cv_mean_pro <- apply(PredProb[, 1:10], 1, cv)
      
      # Calculate mean threshold
      th_mean <- mean(thres, na.rm = TRUE)
      
      # Presence/absence based on mean probability
      mean_pres_abs <- ifelse(mean_pro >= th_mean, 1, 0)
      
      # Create prediction dataframe
      prediction_df <- data.frame(PredProb, mean_pro, cv_mean_pro, predPresAbs[, 1:10], mean_pres_abs)
      
      # Join with subcatchment map
      prediction_map <- map_subcatch %>%
        inner_join(prediction_df, by = c("id" = "stream"))
      
      # Save prediction map
      st_write(prediction_map, paste0(output_path, "/", files_names[sp], "_ssn_predictions.gpkg"),
               driver = "GPKG", append = FALSE, quiet = TRUE)
      
    }, error = function(e) {})
    
  })
  
  # Clean up environment
  rm(list = ls())
  
  # Re-import custom functions
  source("scr_functions.R")
  
  # Re-load spatial data (subcatchment map)
  map_subcatch <- st_read("./data/subc_habit_clasif.gpkg") %>%
    st_set_crs(3795)
  
  # Redefine path to SSN objects
  list_ssn_dirs <- list.dirs(path = "./data/sdm/input/ssn/ssn_obj", recursive = FALSE)
  
  # Redefine file names
  file_names <- list.files(path = "./data/sdm/input/ssn/ssn_obj", full.names = FALSE)
}

