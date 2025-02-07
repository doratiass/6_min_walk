# Impact of Repeated 6-Minute Walk Test Measurements on Mortality Prediction

**Status:** *Work in Progress*

## Overview

This project explores how repeated measurements from the 6-minute walk test (6MWT) can be used to predict mortality risk by employing joint modeling techniques. The analysis integrates longitudinal (repeated measures) and time-to-event (survival) data to generate dynamic, patient-specific risk predictions over time.

**Important:**\
The raw and processed data for this project are confidential and cannot be shared due to privacy restrictions. This repository contains only the code for data cleaning, model development, evaluation, and dynamic predictions.

## Project Contents

The repository is organized into several key components:

-   **Data Preparation and Cleaning**\
    *Scripts:* `scripts/01_data_preparation.R`\
    This section of the code handles:
    -   Importing and cleaning raw clinical and mortality data.
    -   Merging datasets and handling missing values.
    -   Creating derived variables (e.g., follow-up times, 6MWT distance adjustments).
    -   Generating a baseline summary table for patient characteristics.
-   **Model Development and Evaluation**\
    *Scripts:* `scripts/02_model_development.R`\
    This part includes:
    -   Fitting linear mixed effects models (LME) to capture the trajectory of 6MWT measurements.
    -   Developing Cox proportional hazards models to assess survival outcomes.
    -   Constructing joint models that link the longitudinal process with survival outcomes.
    -   Generating time-dependent ROC curves and calibration plots to evaluate model discrimination and calibration.
-   **Dynamic Predictions and Bootstrap Evaluation**\
    *Scripts:* Included within the evaluation code and helper functions in `scripts/00 - funcs.R`\
    Features include:
    -   Generating dynamic prediction plots that provide individualized risk estimates over time.
    -   Performing bootstrap evaluations using parallel processing (with `future` and `furrr`) to assess the stability of model performance metrics (e.g., AUC).

## How to Use the Code

Since the data are confidential, this repository is intended for internal use and further development. The code is modular and can be adapted if you obtain similar data in a secure environment. In summary:

1.  **Data Preparation:** Run the data preparation script to process and merge raw data (if you have access to similar datasets). This will generate an RData file (`CHF_data.RData`) with cleaned data ready for modeling.

2.  **Modeling:** Use the model development scripts to fit the LME, Cox, and joint models. The code includes evaluation components that generate ROC and calibration curves.

3.  **Dynamic Predictions:** The code in the helper functions and dynamic prediction section creates plots that display patient-specific risk trajectories over time. These can be customized or extended to suit your analysis needs.

4.  **Bootstrap Evaluation:** The bootstrap code demonstrates how to assess the stability of the joint model’s discrimination metrics using resampling techniques and parallel processing.

## Requirements

The project relies on the following R packages:

-   **tidyverse:** Data manipulation and visualization.
-   **readxl:** Reading Excel files.
-   **zoo:** Time-series handling.
-   **survival:** Survival analysis.
-   **JMbayes2:** Joint modeling of longitudinal and survival data.
-   **lmtest, MuMIn:** Model testing and evaluation.
-   **future, furrr:** Parallel processing and bootstrap evaluation.
-   **ggpubr, gridExtra, flextable, gtsummary:** Visualization and reporting.

Install the required packages using:

\`\`\`r install.packages(c("tidyverse", "readxl", "zoo", "survival", "future", "furrr", "ggpubr", "gridExtra", "flextable", "gtsummary")) \# Note: JMbayes2 and MuMIn may require installation from CRAN or GitHub.

**Dynamic Predictions**

Dynamic prediction plots are generated using a custom function (plot_dyn_pred()) defined in scripts/00 - funcs.R. These plots provide patient-specific risk predictions over time based on the outputs of the joint models.

**Future Work**

-   **Model Refinement:** Further refine the joint models and explore additional functional forms to improve prediction accuracy.

-   **External Validation:** Validate the models using external datasets or additional cross-validation folds.

-   **Enhanced Reporting:** Develop interactive dashboards or Shiny apps for real-time dynamic predictions.

-   **Manuscript Preparation:** Prepare documentation and manuscripts for potential publication.

------------------------------------------------------------------------

Note: Due to privacy restrictions, the raw and processed data used in this project are not publicly available.
