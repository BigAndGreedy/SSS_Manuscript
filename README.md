# Data Files

## westMN.csv

Three fields, each describing a feature of an observed bigmouth buffalo: 1) Age (in years), 2) Length (total length in mm), and 3) Sex (M for male, F for female, blank cell for unsexed)

## jamestown.csv

Three fields, each describing a feature of an observed bigmouth buffalo: 1) Age (in years), 2) Length (total length in mm), and 3) Sex (M for male, F for female, blank cell for unsexed)

# R Files

## setup.R

Loading libraries and importing data files

## helper_functions.R

-   Defines functions to generate samples via a certain sampling design from a dataset, those being generate_random_sample(), generate_stratified_sample(), and generate_two_phase_sample().
-   Defines VBGF(age, L_inf, k, t0), which gives length at age for given values of parameters. Also defines the function fit_VBGF(), which fits a sex-specific VBGF (with t0 restricted to 0) to a dataset and indicates if it successfully converged (return field Converged = 1 if it successfully converged). If it didn't, returns Converged = 0 and
-   Defines miscellaneous helper functions, get_statistics(), stats_rand_2phase(), and stats_stratified().

## simulated_population.R

Code used to create simulated populations of bigmouth buffalo and compare random/two-phase sampling of those populations. Defines functions for simulating populations of bigmouth buffalo: successful_years(), VBGF(), length_age_norm(), and generate_pop(). Functions were defined to assist simulating samples, gen_pop_conv(), vbgf_funcs2(), and rand_twoPhase_simPops(). Sets of parameters used to simulate populations were defined as jt_pars2, mn_pars2, system3_pars, and system4_pars. Sets of parameters used to simulate samples of different sizes were defined as pars_small, pars_med, and pars_large.

## PA_FA_comparison.R

Code used to compare random, proportional allocation, and fixed allocation samples of empirical datasets and to report results as formatted tables.

## RS_TP_comparison.R

Code used to compare random and two-phase samples of empirical datasets and to report results as formatted tables.

## catch_curve_analysis.R

Code for conducting catch curve analysis on simulated samples from empirical datasets. Analyzes Type 1 and Type 2 samples on empirical datasets.

# Result Files

Files created to store results from long simulations

## jt_MSEs.csv (or mn/s3/s4_MSEs.csv)

Stores results for average MSE of parameter estimates of random/two phase samples in simulated systems using jamestown/minnesota/system 3/system 4 parameters

