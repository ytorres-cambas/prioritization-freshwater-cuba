# Spatial Conservation Prioritization for Freshwater Biodiversity in Cuba

## Overview
This repository contains scripts for spatial conservation prioritization of freshwater biodiversity in Cuba. The analysis uses species distribution models (SDMs) to predict the distribution of freshwater species. Predictions from these models are assessed for overlap with existing protected areas to identify conservation gaps, and are then used as features in the conservation prioritization analysis.

The main objectives of this work are:

- Assess the Effectiveness of Current Protected Areas: Evaluate how well current protected areas provide protection for various freshwater species groups.

- Conservation Planning Under Post-2020 Global Biodiversity Framework: Develop solutions that align with the 30% conservation target for freshwater species, ensuring connectivity.

## Code
The repository includes the following scripts:

- `BRT.R` - Boosted Regression Trees (BRT) for species distribution modeling.
- `MAXENT.R` - Maximum Entropy (MaxEnt) modeling for SDMs.
- `RF.R` - Random Forest (RF) modeling for SDMs.
- `SSN.R` - Spatial Stream Network analysis for freshwater ecosystems.
- `ENSEMBLE.R` - Ensemble modeling approach for SDMs.
- `scr_cons_plann.R` - Main script for conservation planning analysis.
- `scr_functions.R` - Utility functions for SDMs and conservation planning.
- `conserv_prior_cuba.Rproj` - R project file for managing the analysis environment.

## Instructions to run the analysis 
### Download input data

1. Download dataset data.zip.
2. Unzip data.zip and ensure that unziped folder is in the same folder as the R project file `conserv_prior_cuba.Rproj`.
3. Unzip ssn_obj.zip .

### Run the SDMs

1. Open `conserv_prior_cuba.Rproj` in RStudio.
2. Open the scripts to run the SDMs (`BRT.R`, `MAXENT.R`, `RF.R`, `SSN.R`). 
Please note that running all models for all species in the dataset can take a significant amount of 
time (e.g., approximately one week on a computer with 32 GB RAM and 16 cores), especially the SSN model.
3. Run `ENSEMBLE.R` after the last modeling algorith has ended.

### Run the spatial conservation prioritization analysis

1. Run `scr_mk_input_cons_plann.R`
2. Run `scr_cons_plann_30.R`
2. Run `scr_cons_plann_17.R`

## Data 

Input data to run the models and spatial conservation prioritization analyses are orginized as follows:

```
├── data.zip/
│   ├── subc_habit_clasif.gpkg
│   ├── sdm/
│      ├── input/
│         ├── occurr/
│         ├── pred_sites/pred_sites.csv
│         ├── ssn/
│            ├── ssn_obj.zip 
│            ├── CorMdls.RData
│      ├── output/
│         ├── ensemble/
│         ├── indep_models/
│   ├── conservation_prioritization/
│      ├── input/
│         ├── feature.csv
│         ├── longitudinal_distance.csv
│         ├── pu.csv
│      ├── output

Output folders store SDM and spatial conservation prioritization results.

```

### Data Description

- subc_habit_clasif.gpkg: GeoPackage file with spatial units used for the analysis.

- occurr: Contains 230 folders, one for each species. Folder names are in the format genus_species. Files in each folder:

- calib.csv: Table with 10 columns, each containing IDs of spatial units used to train models following a ten-fold repeated split-sampling strategy.

- eval.csv: Same format as calib.csv, but used for model evaluation.

- occurr.csv: Table with 20 columns, including:

   Column 1 (stream): Spatial unit IDs.

   Columns 2-11 (PC1-PC10): Principal components from Principal Component Analysis (PCA).

   Columns 12-13 (z_d_roads, z_d_town): Distance to road and distance to town from the midpoint of each stream segment.

   Column 14 (pres_abs): Presence-absence data.

   Columns 15-20: Taxonomic information (species, family, order, class, phylum, kingdom).

- pred_sites.csv: Table with predictor values (10 components from PCA) for each spatial unit in the study area used for predictions. Columns:

   Columns 1-10: PC1 to PC10 – Principal components from Principal Component Analysis (PCA).

   Columns 12-13: z_d_roads, z_d_town – Distance to road and distance to town from the midpoint of each spatial unit.

- ssn_obj.zip: Files to create an SSN (Spatial Stream Network) object for each species. This object is an R object required as input to run a Spatial Stream Network model.

- CorMdls.RData: R object containing a list of spatial autocorrelation models.

- feature.csv: Table with two columns:

- id_feature: A unique identification number for each feature (species) in the spatial conservation prioritization analysis.

- name: Species names in the format genus_species.

- longitudinal_distance.csv: Table with longitudinal distances between stream segments within the same river basin. Columns:

- id1 and id2: Unique IDs for each stream segment.

- dist_km and dist_km_inv: Distance and inverse distance between id1 and id2, expressed in kilometers.

- pu.csv: Table with information about each planning unit. Columns:

- id: Unique ID for each planning unit.

- cost: Value of the Human Footprint Index.

- locked_in_pa: TRUE if a planning unit is located in a protected area, otherwise FALSE.



## Figures




## Usage
To reproduce the analyses:
1. Clone this repository:
   ```sh
   git clone https://github.com/ytorres-cambas/prioritization-freshwater-cuba
   ```
2. Open the R project file (`conserv_prior_cuba.Rproj`).
3. Run the scripts sequentially according to the workflow.

## Acknowledgments
This work was funded by the Georg Foster Postdoctoral Fellowship of the Alexander von Humboldt Foundation (Ref 3.2 - CUB - 1212347 - GF-P). We also aknowledge support by the German Federal Ministry of Education and Research (BMBF grant agreement number no. 033W034A)

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

