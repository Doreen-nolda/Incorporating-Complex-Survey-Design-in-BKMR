A Practical Framework for Incorporating Complex Survey Design in BKMR

This repository contains code for a design-aware Bayesian Kernel Machine Regression (BKMR) framework that integrates survey-weighted 
resampling to account for complex survey designs when modeling nonlinear exposure–response relationships in environmental mixtures.

Contents
	•	Simulation_Code.R
Simulation framework for evaluating the performance of the proposed method under informative sampling.
The current implementation uses 3 exposure variables for computational efficiency, but the code is fully generalizable and can be 
extended to higher dimensions (e.g., 10 exposures) as used in the full study.
	•	Real_Data_Code.R
Application of the method to NHANES data, including preprocessing, exposure standardization, survey design handling, resampling-based 
BKMR fitting, and visualization of exposure–response functions.

Methods Overview

The workflow:
	•	Generates or processes data with complex survey structure
	•	Applies PSU-preserving, weight-proportional resampling
	•	Fits BKMR models under both naïve and design-aware approaches
Produces: These outputs describe individual, joint, and overall effects of exposures, along with their relative importance in the mixture model.
	•	Univariate exposure–response functions
	•	Bivariate interaction surfaces
	•	Overall mixture effect estimates
	•	Posterior inclusion probabilities (PIPs)

Purpose
The goal of this repository is to provide a practical and reproducible implementation of BKMR that accounts for informative sampling, 
addressing a key limitation of standard BKMR in survey-based studies.

Requirements
	•	R (≥ 4.0)
	•	Packages:
bkmr, dplyr, tidyr, ggplot2, survey, readr

Notes
	•	Simulation settings can be modified to increase the number of exposures.
	•	 Due to computational demands, BKMR may require substantial runtime, especially under resampling, bootstrap replication, or higher-dimensional mixture settings.
