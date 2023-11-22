# Racing Diffusion Model of the Implicit Association Test

## Overview

This repo contains all processing, modeling, visualization, analysis scripts in support of the paper: LaFollette K.J., Rubez D., Demaree H.A., & Goldenberg A. (2023) "Challenging the Mechanism for the Implicit Assoication Test".

Implicit Association Test data collected from the Ideology 2.0 Project (https://osf.io/2483h/; Schmidt et al., 2022). Directory included for pilot / exploratory analyses conducted at the time of Stage 1 submission. Separate directory for confirmatory analyses to be added for Stage 2 submission.

## Contents

### processing.ipynb

Processing notebook to calculate D-scores from raw trial-level data and to extract explicit preference scores from subject-level data. Note that the full csvs containing all raw data from Ideology 2.0 are not uploaded at the time of Stage 1 submission, but csvs containing all data necessary to replicate results are available on OSF https://osf.io/2r3zd/

### raceddm_new.stan and submit.py

Stan modeling code for standard IATs (raceddm_new). submit.py is set to compile model as C++, and fit the model to data with 4 MCMC chains, parallelized across 16 cores (within-chain parallelization accomplished with Stan's reduce_sum function). You will need a C++ compiler such as GCC or Clang, and cmdstanpy installed. Submit separate jobs for each IAT by changing the filename variable in submit.py. For review purposes, already fitted models are available on OSF as netcdf files https://osf.io/2r3zd/

### analysis.ipynb

Visualization and analysis notebook. Reproduces all figures and analyses included in Stage 1 submission.