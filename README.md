<br />
<div align="center">
  <a href="https://github.com/cfjimenezv07/HDFTS_CDF">
    <img src="MQ.png" alt="MQ Logo" height="150">
    <img src="UoY.png" alt="York Logo" height="150">
  </a>

  <h3 align="center">Modeling and Forecasting Subnational Age Distribution of Death Counts</h3>
</div>

## Abstract
<p align="justify">
  Existing mortality forecasting methods focus on age-specific mortality rates, which lie in an unconstrained space, thereby overlooking the distributional nature of life-table death counts. Few studies have developed and compared forecasting methods that model the shape and dynamics of the age distribution of deaths, especially at the subnational level, where data quality varies greatly. This paper presents several forecasting methods to model and forecast the subnational age distribution of death counts. The age distribution of death counts has many similarities to probability density functions, which are non-negative and have a constrained integral, and thus live in a constrained nonlinear space. To address the nonlinear nature of these objects, we implement a cumulative distribution function transformation that is scale-free and additionally imposes monotonicity. Using subnational Japanese life-table death counts from the Japanese Mortality Database (2025), we evaluate the forecast accuracy of the transformation and forecasting methods. The improved forecast accuracy of life-table death counts achieved here will be of great interest to demographers in estimating regional age-specific survival probabilities and life expectancy, and to actuaries as a foundation for exploring potential applications in determining annuity prices for various ages and maturities.
</p>

---

## Data Provenance & Source
The empirical analysis uses subnational Japanese life-table death counts (d_x) for 47 prefectures covering the period 1975–2023 across ages 0 to 110+.

- Source: Japanese Mortality Database (JMD), National Institute of Population and Social Security Research.
- URL: https://www.ipss.go.jp/p-toukei/JMD/index-en.asp
- Data Access Date: September 18, 2025.

---

## Software Dependencies & Environment

This replication package was tested and verified under R version 4.6.0.

### Environment & Package Verification
To reproduce the exact R session environment and inspect pinned package versions:
source("load_packages.R")
sessionInfo()

---

## Code Workflow & Script Execution Guide

To reproduce the analysis, tables, and figures from the paper, execute the R scripts in the following step-by-step order.

### Step 1: Environment Setup & Package Dependencies
Before running any analysis or data scripts, load all necessary R libraries in the current session using:
- load_packages.R: Installs and loads all required package dependencies in the current session.

### Step 2: Data Preprocessing
- read_data.R: Reads raw subnational death probability data, converts q_x to life-table death counts (d_x), computes national benchmark population averages, and exports the .rds list objects (containing 47 prefecture matrices) directly into the data/ directory.

### Step 3: Auxiliary Functions (auxiliary_source/ folder)
The auxiliary_source/ directory contains helper functions used by the point, interval, and evaluation routines:
- CDF_transformation.R: Core function for performing the Cumulative Distribution Function transformation.
- auxiliary_point.R: Helper functions for point forecast computations.
- auxiliary_interval.R: Helper functions for prediction interval estimation.
- forecast_Arima.R: Univariate forecasting routines using ARIMA.
- forecast_errors.R: Calculation routines for evaluating point and interval forecast accuracy metrics.

### Step 4: Point & Interval Forecasting Methods

Once data preparation (read_data.R) is finished, execute the forecasting routines. Output files are automatically saved into designated subdirectories within results/.

> Note: Inside the point/ and interval/ folders, each R script contains a header description detailing which specific paper results, tables, or figures it corresponds to.

#### 4.1 Point Forecasts (point/ folder)
Executes point forecasting models across five methods:
- CDF_UFTS.R: Univariate Functional Time Series method.
- CDF_MFTS.R: Multivariate Functional Time Series method.
- CDF_MLFTS.R: Multilevel Functional Time Series method.
- CDF_HDFPCA.R: High-Dimensional Functional Principal Component Analysis method.
- CDF_FANOVA.R: Functional ANOVA combined with MFTS.

Output destination: Output files are stored in results/Results_Point/.

#### 4.2 Interval Forecasts (interval/ folder)
Contains scripts for uncertainty quantification split into two methodology subfolders:

##### A. Conformal Prediction (interval/conformal/)
- conformal_UFTS.R
- conformal_MFTS.R
- conformal_MLFTS.R
- conformal_HDFPCA.R
- conformal_FANOVA.R

##### B. Standard Deviation Tuning (interval/sd/)
- CDF_UFTS_interval.R
- CDF_MFTS_interval.R
- CDF_MLFTS_interval.R
- CDF_HDFPCA_interval.R
- CDF_FANOVA_interval.R

Output destination: Output files are stored in results/Results_Interval/.

---

### Step 5: Tables, Diagnostics, and Visualization

After generating and storing the point and interval results in results/, run the root scripts to reproduce paper tables, diagnostic tests, and figures.

#### Tables & Statistical Analysis
- Table1_PointAnalysis_EVR.R: Main text point forecast accuracy evaluation using EVR criterion.
- Table2_SelectedK.R: Evaluates the eigenvalue ratio criterion (EVR) component selection across methods.
- Table3_IntervalAnalysis_EVR.R: Main text prediction interval coverage and interval score evaluation using EVR.
- Table4_StationaryAnalysis.R: Stationarity and functional conditional heteroskedasticity diagnostics.
- TableA1_PointAnalysis_K6.R: Appendix A point forecast sensitivity analysis (K=6).
- TableB1_IntervalAnalysis_K6.R: Appendix B prediction interval sensitivity analysis (K=6).
- TableC1_ConformalAnalysis.R: Appendix C conformal prediction interval evaluation.

#### Visualizations & Figures
- Figure_1_and_3.R: Rainbow plots of subnational mortality distributions and Kullback–Leibler divergence heatmaps across prefectures.
- Figure_5_Heatmaps_point.R: Heatmaps comparing point forecast error metrics across prefectures and horizons.
- Figure_6_Heatmaps_Interval.R: Heatmaps evaluating prediction interval performance (Coverage and Interval Score).
- Figure_7_Example_Prediction_Intervals.R: Out-of-sample prediction interval visualizations for representative cases.
- Figure_8_Equiv_MortalityInstruments.R: Visual comparison of forecast d_x, q_x, and e_x for Hokkaido (2007–2023).
- MCS_FigureD1.R: Model Confidence Set (MCS) statistical evaluation plot for Appendix D (Figure D1).

---


## Contact

- Han Lin Shang — hanlin.shang@mq.edu.au
- Cristian Felipe Jiménez-Varón — cristian.jimenezvaron@york.ac.uk

<br />
<div align="center">
  <a href="https://github.com/cfjimenezv07/HDFTS_CDF"><strong>Explore R Code Repository »</strong></a>
</div>