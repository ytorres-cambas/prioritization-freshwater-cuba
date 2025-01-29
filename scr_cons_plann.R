# This script is part of a conservation prioritization analysis. 
# It uses the 'prioritizr' and 'gurobi' packages to solve optimization problems 
# related to spatial conservation planning. The goal is to identify the best 
# set of planning units to conserve biodiversity while considering connectivity 
# between units. This script tests various connectivity penalties and computes 
# solution performance. It also exports results, including maps of irreplaceability 
# and feature representation in protected areas.

# Install the Gurobi optimizer (if needed)
install.packages("C:/gurobi952/win64/R/gurobi_9.5-2.zip", repos = NULL)

# Load required libraries
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

# Input and output paths for data files
inputpath <- "./data/conservation_prioritization/input/"
outputpath <- "./data/conservation_prioritization/output/" 

# Geopackage with spatial units (sub-catchments)
subcatchments <- st_read("./data/subc_habit_clasif.gpkg") %>%
  st_set_crs(3795)

# The analysis starts here

#### Import data

# 1. Import datasets
# Planning units:
pu <- read.csv(paste0(inputpath, "pu.csv"))

# Features: ids and names
features <- read.csv(paste0(inputpath, "feature.csv"))
names(features)[1] <- "id"

# Data.frame with planning unit identifiers, feature identifiers, and
# the amount of the feature in the planning unit
rij <- read.csv(paste0(inputpath, "rij.csv"))

# Data.frame with distances between planning units, used for connectivity
conn_data <- read.csv(paste0(inputpath, "longitudinal_distance.csv"))

# Check if planning units in "connectivity" are present in the planning unit data.frame
conn_data <- conn_data %>%
  filter(id1 %in% pu$id & id2 %in% pu$id)

# Change column names
conn_data <- conn_data %>%
  dplyr::select(-c(dist_km))
names(conn_data)[3] <- "boundary"

# Rescale distance between planning units to avoid extremely long running times
conn_data$boundary <- rescale(conn_data$boundary, to = c(0.1, 1))
conn_data$boundary <- round(conn_data$boundary, digits = 2)

###### Calibrate connectivity penalties (30% target)

# 1. Find an upper limit for connectivity penalties. Run this with high gap and time limit

# Set of penalties to test
prelim_lower <- 0.01
prelim_upper <- 5
penalty_preliminary <- round(seq(prelim_lower, prelim_upper, length.out = 10), 3)

# Set the problem
p_30_prelim_penal <- problem(pu, features, cost_colum = "cost", rij = rij) %>% 
  add_min_set_objective() %>% 
  add_relative_targets(0.3) %>% 
  add_binary_decisions()

# Find solutions
sol_30_prelim_penal <- mclapply(penalty_preliminary,
                                test_penal,
                                p_30_prelim_penal,
                                conn_data,
                                0.2,
                                60*10)

# Format results as a single data.frame
sol_30_prelim_penal <- do.call(bind_cols, sol_30_prelim_penal)
sol_30_prelim_penal$id <- pu$id

# Export results
export_sol(sol_30_prelim_penal, subcatchments, importance = NULL, paste0(outputpath, "sol_30_prelim_penal.gpkg"))

# 2. New prioritizations with candidate penalties. Setting 0.5 as upper limit
penalty <- round(seq(0.01, penalty_preliminary[2], length.out = 20), 3)

# Set the problem
p_30_candidat_penal <- problem(pu, features, cost_colum = "cost", rij = rij) %>% 
  add_min_set_objective() %>% 
  add_relative_targets(0.3) %>% 
  add_binary_decisions()

# Find solutions. No time limit and 10% gap
sol_30_candidat_penal <- mclapply(penalty,
                                  test_penal,
                                  p_30_candidat_penal,
                                  conn_data,
                                  0.1)

# Format results as a single data.frame
sol_30_candidat_penal <- do.call(bind_cols, sol_30_candidat_penal)
sol_30_candidat_penal$id <- pu$id

# Export results
export_sol(sol_30_candidat_penal, subcatchments, importance = NULL, paste0(outputpath,"sol_30_candidat_penal.gpkg"))

# 3. New candidate test. This time using 0.068 as upper limit 
penalty_2 <- round(seq(0.01, penalty[3], length.out = 20), 3)

p_30_candidat_penal_2 <- problem(pu, features, cost_colum = "cost", rij = rij) %>% 
  add_min_set_objective() %>% 
  add_relative_targets(0.3) %>% 
  add_binary_decisions()

# Find solutions. No time limit and 10% gap
sol_30_candidat_penal_2 <- mclapply(penalty_2,
                                    test_penal,
                                    p_30_candidat_penal_2,
                                    conn_data,
                                    0.1)

# Format results as a single data.frame
sol_30_candidat_penal_2 <- do.call(bind_cols, sol_30_candidat_penal_2)
sol_30_candidat_penal$id <- pu$id

# Export results
export_sol(sol_30_candidat_penal_2, subcatchments, importance = NULL, paste0(outputpath,"sol_30_candidat_penal_2.gpkg"))

# Performance of solutions with different connectivity penalties
perform_sol_30_candidat_penal_2 <- conn_metrics(problem = p_30_candidat_penal_2,
                                                solution = select(sol_30_candidat_penal_2, -id),
                                                connectivity_data = conn_data)

# Save performance
colnames(perform_sol_30_candidat_penal_2) <- c("total_cost", "total_connectivity")
perform_sol_30_candidat_penal_2$penalty <- colnames(select(sol_30_candidat_penal_2, -id))
write.csv(perform_sol_30_candidat_penal_2, paste0(outputpath,"perform_sol_30_candidat_penal_2.csv"))

# 4. Visual method to select penalties among candidate penalties
# Find the "elbow" of the curve
connect_cost <- rbind(perform_sol_30_penalty_test, perform_sol_30_candidat_penal_2)

elbow <- find_curve_elbow(data_frame = cbind(connect_cost$total_connectivity,
                                             connect_cost$total_cost),
                          export_type = "row_num",
                          plot_curve = FALSE)

con_penal_30 <- connect_cost %>%
  slice(elbow) %>%
  pull(penalty) %>%
  word(2, sep = "_") %>%
  as.numeric()

# 5. Run prioritization with connectivity penalty = con_penal_30 (selected with the visual method). Use a gap of 5%
# Set the problem
p_30_penaly_elbow <- problem(pu, features, cost_colum = "cost", rij = rij) %>% 
  add_min_set_objective() %>% 
  add_relative_targets(0.3) %>% 
  add_binary_decisions()

# Solve the problem
s_30_penaly_elbow <- p_30_penaly_elbow %>%
  add_default_solver(gap = 0.05) %>%
  add_connectivity_penalties(penalty = con_penal_30,
                             data = conn_data) %>%
  solve()

# Export results
export_sol(s_30_penaly_elbow, subcatchments, importance = NULL, paste0(outputpath,"s_30_penaly_elbow.gpkg"))

# Calculate total cost and connectivity of the solution
total_cost_30 <- eval_cost_summary(p_30_penaly_elbow,
                                   s_30_penaly_elbow["solution_1"]) %>%
  pull(cost)

total_conn_30 <- eval_connectivity_summary(p_30_penaly_elbow,
                                           s_30_penaly_elbow["solution_1"],
                                           data = conn_data) %>%
  pull(connectivity)

# Save cost and connectivity
write.csv(data.frame(total_cost_30, total_conn_30),
          row.names = FALSE, paste0(outputpath,
                                    "cost_conn_30_elbow.csv"))

# Planning units importance
# Calculate irreplaceability of planning units
irr_p_30_penaly_elbow <- eval_ferrier_importance(p_30_penaly_elbow, s_30_penaly_elbow["solution_1"])
irr_p_30_penaly_elbow$id <- pu$id

# Save irreplaceability
write.csv(irr_p_30_penaly_elbow,
          row.names = F, paste0(outputpath,
                                "irr_p_30_penaly_elbow.csv"))

# Map with irreplaceability
sf_irr_p_30_penaly_elbow <- subcatchments %>%
  inner_join(irr_p_30_penaly_elbow, by = c("id" = "id"))

st_write(sf_irr_p_30_penaly_elbow, paste0(outputpath,
                                          "irr_p_30_penaly_elbow.gpkg"))


# feature representation in protected areas and the solution
repr_30_penaly_elbow <- target_coverage(problem = p_30_penaly_elbow,
                                        solution = s_30_penaly_elbow, 
                                        pa_colname = "locked_in_pa")
# Save feature representation
write.csv(repr_30_penaly_elbow$representation, 
          paste0(outputpath, "repr_sol_30_elbow.csv"))

# 6. Problem with protected areas locked-in, 30% target, 0.05 gap, penalty = con_penal_30 (selected with the visual method)

# Set the problem
p_30_penaly_elbow_lock_PA <- problem(pu, features, cost_colum = "cost", rij = rij) %>% 
  add_min_set_objective() %>% 
  add_locked_in_constraints(locked_in = "locked_in_pa") %>%
  add_relative_targets(0.3) %>% 
  add_binary_decisions()

# Solve the problem
s_30_penaly_elbow_lock_PA <- p_30_penaly_elbow_lock_PA %>%
  add_default_solver(gap = 0.05) %>%
  add_connectivity_penalties(penalty = con_penal_30, data = conn_data) %>%
  solve()

# Export results
export_sol(s_30_penaly_elbow_lock_PA, subcatchments, importance = NULL,
           paste0(outputpath,"s_30_penaly_elbow_lock_PA2.gpkg"))

# Calculate total cost and connectivity of the solution
total_cost_lock_PA <- eval_cost_summary(p_30_penaly_elbow_lock_PA,
                                        s_30_penaly_elbow_lock_PA["solution_1"]) %>%
  pull(cost)

total_conn_lock_PA <- eval_connectivity_summary(p_30_penaly_elbow_lock_PA,
                                                s_30_penaly_elbow_lock_PA["solution_1"],
                                                data = conn_data) %>%
  pull(connectivity)

# Save cost and connectivity
write.csv(data.frame(total_cost_lock_PA, total_conn_lock_PA),
          row.names = FALSE, 
          paste0(outputpath, "cost_conn_30_elbow_lock_PA.csv"))

# Feature representation in protected areas and the solution
repr_30_penaly_elbow_lock_PA <- target_coverage(problem = p_30_penaly_elbow_lock_PA,
                                                solution = s_30_penaly_elbow_lock_PA, 
                                                pa_colname = "locked_in_pa")

# Save feature representation
write.csv(repr_30_penaly_elbow_lock_PA$representation, 
          row.names = FALSE, 
          paste0(outputpath, "repr_sol_30_elbow_lock_PA.csv"))

# Planning units importance
# Calculate irreplaceability of planning units
irr_30_penaly_elbow_lock_PA <- eval_ferrier_importance(p_30_penaly_elbow_lock_PA, 
                                                       s_30_penaly_elbow_lock_PA["solution_1"])

irr_30_penaly_elbow_lock_PA$id <- pu$id

# Save irreplaceability
write.csv(irr_30_penaly_elbow_lock_PA,
          row.names = F, 
          paste0(outputpath, "irr_30_penaly_elbow_lock_PA.csv"))

# Map with irreplaceability
sf_irr_30_penaly_elbow_lock_PA <- subcatchments %>%
  inner_join(irr_30_penaly_elbow_lock_PA, by = c("id" = "id"))

st_write(sf_irr_30_penaly_elbow_lock_PA, 
         paste0(outputpath, "irr_30_penaly_elbow_lock_PA.gpkg"))

