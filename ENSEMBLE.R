# Ensemble Model for Species Distribution Prediction

# This script performs an ensemble modeling approach to predict species presence-absence
# using multiple independent model outputs. The final ensemble model aggregates results
# from individual models and evaluates performance. The output includes averaged
# presence-absence predictions, model performance
# metrics, threshold values, and spatial predictions mapped to sub-catchments.

# Load required libraries
library(dplyr)
library(sf)
library(raster)
source("scr_functions.R")

# Define input and output paths
pathInput <- "./data/sdm/output/indep_models/"   # Path to independent model predictions
pathOutput <- "./data/sdm/output/ensemble/"      # Path to save ensemble results

# Get the names of species model files
files_names <- list.files(path = pathInput, full.names = FALSE)

# Load sub-catchment map with spatial information
map_subcatch <- st_read("./subc_habit_clasif.gpkg", quiet = TRUE) %>%
  st_set_crs(3795)  # Set coordinate reference system

# Loop through each species model and create ensemble predictions
for (sp in seq_along(files_names)) {
  tryCatch({
    
    # Load stream IDs for prediction sites
    ids <- read.csv("./data/sdm/input/pred_sites/pred_sites.csv")$stream_ID
    
    # Load presence-absence predictions from independent models
    filesPresAbs <- list.files(
      path = file.path(pathInput, files_names[sp]), 
      pattern = ".*presabs.*\\.csv$", full.names = TRUE
    )
    presabs <- lapply(filesPresAbs, read.csv) %>% data.frame()
    
    # Load occurrence records for the species
    occurr <- read.csv(file.path("./data/sdm/input/occurr/", files_names[sp], "occurr.csv"))
    occurr2 <- occurr %>% dplyr::select(stream, pres_abs)
    
    # Apply ensemble modeling function
    ens <- ensem_mod(ids, presabs, occurr2)
    
    # Define output path for the species
    outputpath2 <- file.path(pathOutput, files_names[sp])
    if (!dir.exists(outputpath2)) dir.create(outputpath2, recursive = TRUE)
    
    # Export results
    write.csv(ens$avg_pres_abs, file.path(outputpath2, "avg_pres_abs.csv"), row.names = FALSE)
    write.csv(ens$ensem_presabs, file.path(outputpath2, "ensem_presabs.csv"), row.names = FALSE)
    write.csv(ens$eval_ensem, file.path(outputpath2, "ensem_performance.csv"), row.names = FALSE)
    write.csv(ens$thresh_ensem, file.path(outputpath2, "ensem_threshold.csv"), row.names = FALSE)
    
    # Generate spatial output
    sf_presabs <- map_subcatch %>% inner_join(ens$ensem_presabs, by = c("id" = "ids"))
    st_write(sf_presabs, file.path(outputpath2, paste0(files_names[sp], "_ensem.gpkg")),
             driver = "GPKG", append = FALSE, quiet = TRUE)
    
    print(paste0("Exported results for ", files_names[sp]))
    
  }, error = function(e) {})  # Handle errors silently
}
