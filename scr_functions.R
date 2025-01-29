# Set of functions called by scripts that execute species distribution modelling
# and spatial conservation planning

library(prioritizr)
library(gurobi)
library(dplyr)
library(sf)
library(withr)
library(parallel)
library(foreach)
library(ggplot2)
library(tmap)
library(scales)
library(pathviewr)
library(units)
library(stringr)

#' Function to export a solution as a GeoPackage
#'
#' This function exports a solution to a GeoPackage (GPKG) format, optionally 
#' adding an importance column to the solution data. The solution is joined 
#' with subcatchments based on a common identifier and written to the specified 
#' output path.
#'
#' @param solution A data frame containing the solution, with a column for 
#'   the identifier (`id`) to match with the `subcatchments` data.
#' @param subcatchments A spatial data frame or Simple Features (SF) object 
#'   representing subcatchments, containing the common identifier (`id`).
#' @param importance An optional numeric vector representing the importance values 
#'   to be added to the solution data.
#' @param output_path A character string specifying the path to save the GeoPackage 
#'   (GPKG) file, including the filename and extension.
#'
#' @return A GeoPackage file containing the joined solution and subcatchments data.
#' @export
#' @import dplyr
#' @import sf
#'
export_sol <- function(solution, subcatchments, importance = NULL, output_path) {
  
  # Join solution data with subcatchments
  map <- subcatchments %>%
    inner_join(solution, by = c("id" = "id"))
  
  # If importance is provided, add it to the map
  if (!is.null(importance)) {
    map$importance <- importance
  } 
  
  # Write the result to a GeoPackage (GPKG)
  st_write(map, output_path, append = FALSE)
}


#' Function to assess feature representation
#'
#' This function calculates and compares the representation of features 
#' based on existing protected areas and prioritization solutions. 
#' It computes the number of features adequately represented by the 
#' existing protected areas and the prioritization solution, as well as 
#' the relative coverage for each.
#'
#' @param problem A prioritization problem object containing the data for the features.
#' @param solution A data frame representing the prioritization solution, with a column indicating the protected areas.
#' @param pa_colname A character string representing the column name in the solution that indicates existing protected areas.
#'
#' @return A list containing:
#' \item{n_features_protected}{A data frame with the number of features protected by the existing protected areas and by the prioritization solution.}
#' \item{representation}{A data frame showing the species, their relative protection in existing protected areas, and in the prioritization solution, along with whether the target is met.}
#' @export
#' @import prioritizr
#'
target_coverage <- function(problem, solution, pa_colname) {
  
  # Calculate feature representation based on existing protected areas
  pa <- as.data.frame(ifelse(solution[, pa_colname] == TRUE, 1, 0))
  tc_pa <- eval_target_coverage_summary(problem, pa)
  
  # Number of features adequately represented by existing protected areas
  n_feature_protected_pa <- sum(tc_pa$met)
  
  # Summarize representation (percent coverage)
  relative_protection_pa <- tc_pa$relative_held * 100
  
  # Calculate feature representation based on the prioritization solution
  tc_sol <- eval_target_coverage_summary(problem, 
                                         as.data.frame(solution[, "solution_1"]))
  
  # Number of features adequately represented by the solution
  n_feature_protected_solution <- sum(tc_sol$met)
  
  # Summarize representation (percent coverage)
  relative_protection_solution <- tc_sol$relative_held * 100
  
  # Prepare the output
  output <- list(
    "n_features_protected" = data.frame(n_feature_protected_pa, 
                                        n_feature_protected_solution),
    "representation" = data.frame(
      "species" = tc_pa$feature,
      relative_protection_pa,
      "met_target_pa" = tc_pa$met,
      relative_protection_solution,
      "met_target_sol" = tc_sol$met)
  )
  
  return(output)
}



#' Function to visualize representation of features
#'
#' This function creates a histogram to visualize the comparison between existing protected areas (PA)
#' and prioritization solutions based on the coverage of features. A vertical line indicating the target 
#' coverage is also plotted.
#'
#' @param represent A data frame containing the representation data for both protected areas and prioritization solutions.
#' @param col_prot_pa A character string representing the column name for the protected areas coverage in the data.
#' @param col_prot_sol A character string representing the column name for the prioritization solution coverage in the data.
#' @param target A numeric value representing the target percentage coverage of features to be plotted as a vertical line.
#'
#' @return A ggplot object displaying the histogram of coverage comparisons.
#' @export
#' @import ggplot2
#' @importFrom RColorBrewer scale_color_brewer
#'
hist_representation <- function(represent, col_prot_pa, col_prot_sol, target) {
  
  # Combine the existing PA and prioritization data into a representation vector
  representation <- c(rep("Existing PA", times = nrow(represent)), 
                      rep("Prioritization", times = nrow(represent)))
  
  # Combine the coverage values for plotting
  coverage <- c(represent[, col_prot_pa], represent[, col_prot_sol])
  
  # Create a data frame for plotting
  df_rep <- data.frame(coverage, representation)
  
  # Create a data frame for the target line
  ln <- data.frame("ln" = target) 
  
  # Create the plot
  plot_rep <- ggplot(df_rep, aes(x = coverage, fill = representation, color = representation)) +
    geom_histogram(position = "dodge", fill = "white") +
    xlab("Percent coverage of features (%)") +
    ylab("Frequency") +
    scale_color_brewer(palette = "Dark2") +
    theme_minimal() +
    geom_vline(data = ln, aes(xintercept = ln), linetype = "dashed")
  
  return(plot_rep)
}



#' Function to solve problems with different connectivity penalties
#'
#' This function solves optimization problems by applying different connectivity penalties and
#' using the Gurobi solver. The results are then used to calibrate the connectivity penalties.
#'
#' @param penalty A numeric value representing the connectivity penalty to be applied.
#' @param problem An object representing the optimization problem (e.g., from the \code{prioritizr} package).
#' @param conn_data A data frame containing connectivity data used to apply the penalty.
#' @param gap A numeric value for the optimality gap. This is the stopping criterion for the solver.
#' @param time_limit A numeric value specifying the maximum time (in seconds) for solving the problem.
#'   Defaults to the largest integer possible.
#'
#' @return A data frame with the solution to the problem, with a column for each penalty value.
#' @export
#' @import prioritizr
#' @importFrom dplyr %>%
#'
test_penal <- function(penalty, problem, conn_data, gap, time_limit = .Machine$integer.max) {
  
  # Solve the problem by adding connectivity penalties and applying the solver
  s <- problem %>%
    add_connectivity_penalties(penalty = penalty, data = conn_data) %>%
    add_gurobi_solver(gap = gap, time_limit = time_limit) %>%
    solve()
  
  # Convert the solution to a data frame and format the column name
  s <- data.frame(s = s$solution_1)
  names(s) <- with_options(list(scipen = 30), paste0("penalty_", penalty))
  
  return(s)
}



#' Function to calculate performance of solutions with different connectivity penalties
#'
#' This function calculates the performance of solutions for a given problem by evaluating their total cost
#' and connectivity using parallel processing. It calculates the cost and connectivity for each column of
#' the solution matrix, where each column represents a potential solution.
#'
#' @param problem An object representing the optimization problem (e.g., from the \code{prioritizr} package).
#' @param solution A matrix or data frame where each column is a proposed solution to the problem.
#' @param connectivity_data A data frame containing connectivity information needed for evaluating the solutions.
#' 
#' @return A matrix with two columns: total cost and connectivity for each solution.
#' @export
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' @import prioritizr
#'
conn_metrics <- function(problem, solution, connectivity_data) {
  
  # Detect the number of available cores
  n.cores <- parallel::detectCores() - 69
  
  # Create a cluster for parallel processing
  my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
  
  # Register the cluster for parallel processing
  doParallel::registerDoParallel(cl = my.cluster)
  
  # Run parallel processing for each solution
  conn_pen_perform <- foreach (m = 1:ncol(solution),
                               .packages = "prioritizr",
                               .combine = rbind) %dopar% {
                                 
                                 # Calculate total cost for the current solution
                                 total_cost <- eval_cost_summary(problem, data.frame(solution[, m]))$cost
                                 
                                 # Calculate total connectivity for the current solution
                                 total_conn <- eval_connectivity_summary(problem, 
                                                                         data.frame(solution[, m]), 
                                                                         data = connectivity_data)$connectivity
                                 
                                 # Combine results into a single vector
                                 result <- c(total_cost, total_conn)
                               }
  
  # Stop the cluster
  parallel::stopCluster(cl = my.cluster)
  
  return(conn_pen_perform)
}




#' Function modified from https://github.com/Pakillo/rSDM/blob/master/R/point_in_cell.R
#'
#' Do points (occurrences) fall out of raster cells?
#'
#' This function examines which points fall on a raster cell without data (NA). 
#' Returns TRUE for points falling in a raster cell without data, and FALSE otherwise.
#' @param locs Points represented by a two-column matrix or data.frame
#' @param ras Raster* object
#'
#' @return A logical vector.
#' @export
#'
point_out_cell <- function(locs, ras) {
  
  # Get NA cells
  rasvals <- raster::extract(ras, locs)
  missing <- is.na(rasvals)
  return(missing)
}

#' Function modified from https://github.com/Pakillo/rSDM/blob/master/R/point_in_cell.R
#' Move occurrences to the nearest raster cell with data
#'
#' Move point occurrences falling in raster cells without data (i.e. NA) to the nearest raster cell with data.
#'
#' @export
#' @import raster
#' @import sp
#' @importFrom class knn1
#' @param locs A matrix or data.frame with coordinates.
#' @param ras Raster* object.
#' @param distance Numeric (optional). Maximum distance to move points. Point coordinates are only changed if the distance to the nearest raster cell is below \code{distance}.
#' @param showchanges Logical. Print table with old and new coordinates.
#' @param showmap Logical. Show map with original and new coordinates?
#' @param leaflet Logical. If TRUE, show leaflet map instead of static map.
#' @return A SpatialPointsDataFrame (with corrected coordinates if move is TRUE).
points2nearestcell <- function(locs = NULL, ras = NULL, 
                               distance = NULL, crs) {
  
  miss <- point_out_cell(locs, ras)
  
  # If there are NA cells...
  if (sum(miss) > 0) {
    
    coord.miss <- locs |> 
      dplyr::slice(which(miss == TRUE))  # points without data
    
    cells.notNA <- raster::rasterToPoints(ras, spatial = TRUE)  # get coordinates of cells with data
    coord.ras <- sp::coordinates(cells.notNA)
    cell.id <- factor(seq_len(nrow(coord.ras)))
    
    # Find the nearest raster cell for each point with missing data
    nearest.cell <- class::knn1(coord.ras, coord.miss, cell.id)
    
    new.coords <- matrix(coord.ras[nearest.cell, ], ncol = 2)
    colnames(new.coords) <- c("long_new", "lat_new")
    
    if (isTRUE(move)) {   # assign new coordinates to those 
      locs_old <- locs
      names(locs_old) <- c("long_old", "lat_old")
      locs[miss, ] <- new.coords
    }
    
    if (!is.null(distance)) {
      
      # Calculate distances between old and new coordinates
      distances <- raster::pointDistance(locs_old, locs, lonlat = FALSE)
      
      # If distance below threshold, accept, otherwise keep old coordinates
      locs <- cbind(locs_old, locs, distances)
    }
    
    locs_points <- locs |> 
      dplyr::select(long, lat) |> 
      sf::st_as_sf(coords = c("long", "lat")) |>
      st_set_crs(crs)
    
    locs_old_points <- locs_old |> 
      dplyr::select(long_old, lat_old) |> 
      sf::st_as_sf(coords = c("long_old", "lat_old")) |>
      st_set_crs(crs)
    
  } else {
    message("All points fall within a raster cell")
  }
  
  return(list("new_points_df" = locs, "new_points_sf" = locs_points, 
              "old_points_sf" = locs_old_points))
}



#' Maxent Parameter Tuning
#'
#' This function performs simultaneous tuning of Maxent's regularization multiplier and feature classes using cross-validation.
#'
#' @param data A data frame containing occurrences (presences coded as '1' and pseudoabsences coded as '0') in the first column, with the remaining columns containing predictor values at each occurrence point.
#' @param occ A character string specifying the name of the column in `data` that contains presence/absence data.
#' @param pred A character vector specifying the names of the columns in `data` that contain predictor variables.
#' @param k An integer specifying the number of folds for cross-validation. Default is 5.
#'
#' @details 
#' The function tunes Maxent's regularization multiplier and feature classes using a grid search approach. It evaluates performance through cross-validated AUC scores and returns the best parameter combination.
#'
#' @return A character vector specifying the best combination of regularization multiplier and feature class arguments for Maxent.
#'
#' @examples
#' \dontrun{
#' library(dismo)
#' library(caret)
#' library(precrec)
#' library(rJava)
#'
#' # Example data
#' data <- data.frame(
#'   occurrence = c(1, 0, 1, 0, 1, 0),
#'   predictor1 = runif(6),
#'   predictor2 = runif(6)
#' )
#' 
#' best_params <- maxent_param(data, occ = "occurrence", pred = c("predictor1", "predictor2"))
#' print(best_params)
#' }
#'
#' @import dismo caret precrec rJava
#' @export
maxent_param <- function(data, occ, pred, k = 5) {
  
  # Check and load required packages
  if (!requireNamespace("dismo", quietly = TRUE)) install.packages("dismo")
  if (!requireNamespace("caret", quietly = TRUE)) install.packages("caret")
  if (!requireNamespace("precrec", quietly = TRUE)) install.packages("precrec")
  if (!requireNamespace("rJava", quietly = TRUE)) install.packages("rJava")
  
  library(dismo)
  library(caret)
  library(precrec)
  library(rJava)
  
  # Generate balanced cross-validation folds
  folds <- caret::createFolds(y = as.factor(data[, occ]), k = k)
  
  # Define regularization multipliers and feature classes
  ms <- c(0.5, 1, 2, 3, 4)
  grid <- expand.grid(
    regmult = paste0("betamultiplier=", ms),
    features = list(
      c("noautofeature", "nothreshold"), # LQHP
      c("noautofeature", "nothreshold", "noproduct"), # LQH
      c("noautofeature", "nothreshold", "nohinge", "noproduct"), # LQ
      c("noautofeature", "nothreshold", "nolinear", "noquadratic", "noproduct"), # H
      c("noautofeature", "nothreshold", "noquadratic", "nohinge", "noproduct")  # L
    ),
    stringsAsFactors = FALSE
  )
  
  AUCs <- c()
  
  # Grid search for tuning parameters
  for (n in seq_along(grid[, 1])) {
    full_pred <- data.frame()
    
    for (i in seq_len(length(folds))) {
      trainSet <- unlist(folds[-i])
      testSet <- unlist(folds[i])
      
      # Train Maxent model
      if (inherits(try(
        maxmod <- dismo::maxent(
          x = data[trainSet, pred],
          p = data[trainSet, occ],
          removeDuplicates = FALSE,
          args = as.character(unlist(grid[n, ]))
        )
      ), "try-error")) {
        next
      }
      
      # Predict using the Maxent model
      modpred <- predict(maxmod, data[testSet, pred], args = "outputformat=cloglog")
      pred_df <- data.frame(score = modpred, label = data[testSet, occ])
      full_pred <- rbind(full_pred, pred_df)
    }
    
    # Evaluate AUC
    AUCs[n] <- precrec::auc(precrec::evalmod(scores = full_pred$score, labels = full_pred$label))[1, 4]
  }
  
  # Find the best parameter combination
  best_param <- as.character(unlist(grid[which.max(AUCs), ]))
  return(best_param)
}


#' Evaluate Model Performance
#'
#' This function calculates metrics of model performance, including AUC of ROC and threshold-based metrics such as sensitivity, specificity, omission, commission, and TSS. The function also determines the optimal threshold that minimizes the difference between sensitivity and specificity.
#'
#' @param idsPresabsProb A data frame with three columns:
#'   \itemize{
#'     \item Column 1: IDs (e.g., unique identifiers for observations).
#'     \item Column 2: Observed occurrences (binary: 1 for presence, 0 for absence).
#'     \item Column 3: Predicted probabilities of occurrence.
#'   }
#'
#' @return A list containing the following components:
#'   \item{threshold}{A data frame with optimal thresholds calculated using the `optimal.thresholds` function.}
#'   \item{performance}{A data frame with performance metrics (e.g., sensitivity, specificity, omission, commission, TSS, and AUC).}
#'
#' @examples
#' \dontrun{
#' library(PresenceAbsence)
#' library(modEvA)
#' library(dplyr)
#' library(tibble)
#'
#' # Example data
#' idsPresabsProb <- data.frame(
#'   IDs = c("ID1", "ID2", "ID3"),
#'   Observed = c(1, 0, 1),
#'   Predicted = c(0.8, 0.2, 0.9)
#' )
#'
#' result <- eval_mod(idsPresabsProb)
#' print(result$threshold)
#' print(result$performance)
#' }
#'
#' @import PresenceAbsence modEvA dplyr tibble
#' @export
eval_mod <- function(idsPresabsProb) {
  
  # Check and load required packages
  if (!requireNamespace("PresenceAbsence", quietly = TRUE)) install.packages("PresenceAbsence")
  if (!requireNamespace("modEvA", quietly = TRUE)) install.packages("modEvA")
  if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
  if (!requireNamespace("tibble", quietly = TRUE)) install.packages("tibble")
  
  library(PresenceAbsence)
  library(modEvA)
  library(dplyr)
  library(tibble)
  
  # Calculate AUC of ROC
  AUC_ROC <- AUC(
    obs = idsPresabsProb[, 2], 
    pred = idsPresabsProb[, 3], 
    plot = FALSE,
    simplif = TRUE
  )
  
  # Calculate optimal thresholds
  thresh <- idsPresabsProb %>%
    optimal.thresholds(opt.methods = c("Sens=Spec"))
  
  # Calculate threshold-based metrics
  perform <- threshMeasures(
    obs = idsPresabsProb[, 2], 
    pred = idsPresabsProb[, 3], 
    simplif = TRUE, 
    thresh = thresh[, 2], 
    standardize = TRUE,
    measures = c("Sensitivity", "Specificity", "Omission", "Commission", "TSS")
  ) %>%
    as.data.frame() %>%
    rownames_to_column(var = "Measures") %>%
    add_row(Measures = "AUC_ROC", Value = AUC_ROC)  # Add AUC as a performance metric
  
  # Output results as a list
  func_output <- list(
    "threshold" = thresh,
    "performance" = perform
  )
  
  return(func_output)
}




#' Ensemble Model Evaluation and Thresholding
#'
#' This function calculates the average presence-absence values, evaluates model performance, and thresholds ensemble predictions based on observed occurrences. 
#'
#' @param ids A vector of unique IDs for the spatial units (e.g., stream IDs).
#' @param presabs A data frame or matrix containing presence-absence predictions, where rows correspond to spatial units and columns to models.
#' @param occurr A data frame containing observed occurrences. The first column must match the `ids` vector, and another column should include observed presence-absence values named `pres_abs`.
#'
#' @return A list with the following components:
#'   \item{avg_pres_abs}{A data frame with the average presence-absence values for each spatial unit.}
#'   \item{ensem_presabs}{A data frame with binary presence-absence values based on the threshold.}
#'   \item{eval_ensem}{A data frame with performance metrics for the ensemble predictions.}
#'   \item{thresh_ensem}{The threshold value used to transform average presence-absence values into binary predictions.}
#'
#' @examples
#' \dontrun{
#' ids <- c(1, 2, 3)
#' presabs <- matrix(c(0.2, 0.7, 0.9, 0.3, 0.8, 0.6), nrow = 3, ncol = 2)
#' occurr <- data.frame(stream = c(1, 2), pres_abs = c(1, 0))
#' result <- ensem_mod(ids, presabs, occurr)
#' print(result)
#' }
#'
#' @import dplyr
#' @export
ensem_mod <- function(ids, presabs, occurr) {
  # Check for required package
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    install.packages("dplyr")
  }
  library(dplyr)
  
  ## Average presence-absence
  avg_presabs <- apply(presabs, 1, mean, na.rm = TRUE) %>%
    as.data.frame()
  avg_presabs$ids <- ids
  names(avg_presabs)[1] <- "avg_pres_ab"
  
  ### Evaluation 
  # Select predictions to compare with observed occurrences
  pred_subset <- avg_presabs %>%
    filter(ids %in% occurr[, 1])
  
  # Inputs for the eval_mod() function
  eval_input <- pred_subset %>%
    inner_join(occurr, by = c("ids" = "stream")) %>%
    dplyr::select(ids, pres_abs, avg_pres_ab)
  
  # Performance metrics
  eval_ensem <- eval_mod(eval_input)$performance
  
  # Threshold
  thresh_ensem <- eval_mod(eval_input)$threshold[1, 2]
  
  ### Transform average presence-absence to binary
  ensem_presabs <- ifelse(dplyr::select(avg_presabs, -ids) >= thresh_ensem, 1, 0)
  colnames(ensem_presabs) <- "ensem_presabs"
  
  ensem_presabs <- data.frame(ids, ensem_presabs)
  
  ## Function output
  output <- list(
    "avg_pres_abs" = avg_presabs,
    "ensem_presabs" = ensem_presabs,
    "eval_ensem" = eval_ensem,
    "thresh_ensem" = thresh_ensem
  )
  
  return(output)
}


#' Model Selection for Spatial Stream Network
#'
#' This function performs model selection for Spatial Stream Network (SSN) analysis by fitting non-spatial and spatial models. 
#' It selects the best non-spatial model based on the lowest AUC and evaluates spatial autocorrelation models using AIC.
#'
#' @param SSNobj An SSN object created using the `SSN` package, containing observed points and network topology.
#' @param CorMdls A list of correlation models to test for spatial autocorrelation. Each correlation model should be an SSN-supported model (e.g., `Exponential.tailup`, `Exponential.taildown`).
#'
#' @return A list with the following components:
#'   \item{formula}{The formula of the best-fitting non-spatial model, enhanced with additional predictors.}
#'   \item{model}{The spatial correlation model with the lowest AIC.}
#'   \item{modelsAIC}{A data frame with AIC values for all tested spatial correlation models.}
#'
#' @examples
#' \dontrun{
#' library(SSN)
#' SSNobj <- importSSN("path/to/ssn_directory", predpts = "pred_points")
#' CorMdls <- list("Exponential.tailup", "Exponential.taildown")
#' result <- SSNparam(SSNobj, CorMdls)
#' print(result$formula)
#' print(result$model)
#' print(result$modelsAIC)
#' }
#'
#' @import SSN dplyr foreach doParallel parallel purrr MuMIn
#' @export
SSNparam <- function(SSNobj, CorMdls) {
  # Check for required packages
  if (!requireNamespace("foreach", quietly = TRUE)) install.packages("foreach")
  if (!requireNamespace("doParallel", quietly = TRUE)) install.packages("doParallel")
  if (!requireNamespace("parallel", quietly = TRUE)) install.packages("parallel")
  if (!requireNamespace("MuMIn", quietly = TRUE)) install.packages("MuMIn")
  if (!requireNamespace("purrr", quietly = TRUE)) install.packages("purrr")
  if (!requireNamespace("SSN", quietly = TRUE)) install.packages("SSN")
  if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
  
  library(SSN)
  library(dplyr)
  library(foreach)
  library(doParallel)
  library(parallel)
  library(purrr)
  library(MuMIn)
  
  # Selecting a non-spatial model
  DataNsp <- getSSNdata.frame(SSNobj, Name = "Obs") %>%
    dplyr::select(PC1:PC10, pres_abs)
  
  ModelNsp <- try(glm(formula = pres_abs ~ ., 
                      data = DataNsp,
                      family = binomial))
  
  options(na.action = "na.fail")
  ListModelNsp <- dredge(ModelNsp)
  BestModelNsp <- eval(attributes(ListModelNsp)$model.calls[[1]])
  
  # Get residuals
  SSNobj@obspoints@SSNPoints[[1]]@point.data$RES <- resid(BestModelNsp)
  
  # Spatial models
  # Setup a parallel backend
  n.cores <- parallel::detectCores() - 1
  my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
  doParallel::registerDoParallel(cl = my.cluster)
  
  ssn_cor_test <- foreach(m = 1:length(CorMdls), .packages = "SSN", .errorhandling = "remove") %dopar% {
    glmssn(RES ~ 1, SSNobj, CorModels = eval(CorMdls[[m]]),
           addfunccol = "computed.afv")
  }
  parallel::stopCluster(cl = my.cluster)
  
  # Models AIC
  cor_modl_AIC <- InfoCritCompare(keep(ssn_cor_test, is.list))
  
  # Results
  formula <- BestModelNsp$formula
  formula <- as.character(formula[[3]])
  formula <- formula[2]
  formula <- as.formula(paste("pres_abs ~ ", formula, " + z_d_town + z_d_roads"))
  model <- eval(CorMdls[[which.min(cor_modl_AIC$AIC)]])
  
  return(list("formula" = formula, "model" = model, "modelsAIC" = cor_modl_AIC))
}


