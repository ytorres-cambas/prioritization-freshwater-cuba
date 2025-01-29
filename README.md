# Spatial Conservation Prioritization for Freshwater Biodiversity in Cuba

## Overview
This repository contains scripts for spatial conservation prioritization of freshwater biodiversity in Cuba. The analysis uses species distribution models (SDMs) to predict the distribution of freshwater species. Predictions from these models are assessed for overlap with existing protected areas to identify conservation gaps, and are then used as features in the conservation prioritization analysis.

The main objectives of this work are:

    Assess the Effectiveness of Current Protected Areas: Evaluate how well current protected areas provide protection for various freshwater species groups.
    Conservation Planning Under Post-2020 Global Biodiversity Framework: Develop solutions that align with the 30% conservation target for freshwater species, ensuring connectivity.

## Repository Structure
The repository includes the following scripts:

- `BRT.R` - Boosted Regression Trees (BRT) for species distribution modeling.
- `MAXENT.R` - Maximum Entropy (MaxEnt) modeling for SDMs.
- `RF.R` - Random Forest (RF) modeling for SDMs.
- `SSN.R` - Spatial Stream Network analysis for freshwater ecosystems.
- `ENSEMBLE.R` - Ensemble modeling approach for SDMs.
- `scr_cons_plann.R` - Main script for conservation planning analysis.
- `scr_functions.R` - Utility functions for SDMs and conservation planning.
- `conserv_prior_cuba.Rproj` - R project file for managing the analysis environment.

## Input data
Input data for scripts can be downloaded from 

## Usage
To reproduce the analyses:
1. Clone this repository:
   ```sh
   git clone https://github.com/your-username/your-repo-name.git
   ```
2. Open the R project file (`conserv_prior_cuba.Rproj`).
3. Run the scripts sequentially according to the workflow.

## Acknowledgments
This work was funded by the Georg Foster Postdoctoral Fellowship of the Alexander von Humboldt Foundation (Ref 3.2 - CUB - 1212347 - GF-P). We also aknowledge support by the German Federal Ministry of Education and Research (BMBF grant agreement number no. 033W034A)

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

