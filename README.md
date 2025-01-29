# Spatial Conservation Prioritization for Freshwater Biodiversity in Cuba

## Overview
This repository contains scripts used for spatial conservation prioritization of freshwater biodiversity in Cuba. The analysis focuses on species distribution modeling (SDMs) and systematic conservation planning to identify key areas for conservation efforts.

## Abstract
The study aims to support freshwater biodiversity conservation in Cuba by integrating species distribution models with systematic conservation planning approaches. Using ecological and environmental data, we applied multiple modeling techniques, including Boosted Regression Trees (BRT), Random Forest (RF), Maximum Entropy (MaxEnt), and ensemble modeling, to predict the distribution of freshwater species. The results informed a conservation prioritization analysis to identify high-priority areas for protection. This workflow provides a robust methodological framework for biodiversity conservation planning in freshwater ecosystems.

## Repository Structure
The repository includes the following scripts:

- `BRT.R` - Boosted Regression Trees (BRT) for species distribution modeling.
- `ENSEMBLE.R` - Ensemble modeling approach for SDMs.
- `MAXENT.R` - Maximum Entropy (MaxEnt) modeling for SDMs.
- `RF.R` - Random Forest (RF) modeling for SDMs.
- `SSN.R` - Spatial Stream Network analysis for freshwater ecosystems.
- `scr_cons_plann.R` - Main script for conservation planning analysis.
- `scr_functions.R` - Utility functions for SDMs and conservation planning.
- `conserv_prior_cuba.Rproj` - R project file for managing the analysis environment.

## Usage
To reproduce the analyses:
1. Clone this repository:
   ```sh
   git clone https://github.com/your-username/your-repo-name.git
   ```
2. Open the R project file (`conserv_prior_cuba.Rproj`).
3. Run the scripts sequentially according to the workflow.
4. Modify parameters as needed for specific analyses.

## Dependencies
Ensure the following R packages are installed:
```r
install.packages(c("dismo", "raster", "sf", "ggplot2", "randomForest", "gbm"))
```

## Acknowledgments
This work was supported by the Georg Foster Postdoctoral Fellowship of the Alexander von Humboldt Foundation (Ref 3.2 - CUB - 1212347 - GF-P).

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

