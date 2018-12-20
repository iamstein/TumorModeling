This repository contains the code for the paper that has been uploaded to biorxiv here:
https://www.biorxiv.org/content/early/2018/10/17/432500

## Parameter Files
The parameter files are stored in the data directory and are called "ModelF_DRUG.xlsx" where DRUG can be Atezolizumab, Bevacizumab, Pembrolizumab, or Trastuzumab.

## Key R Scripts and Rmarkdown files
These are short Tasks in Rmarkdown or R for generating key figures for the manuscript or other helpful codes.
They are stored in the ModelF folder.
-	Task10  - illustartes some simple graphs for how code can be executed
- Task09c - Varying dose and checking AFTIR for each drug
- Task05  - Varying other parameters and checking AFTIR accuracy for each drug
- Task09b - comparison of simple formulation of AFTIR to the full formulation of AFTIR
- Task04  - look at how changing target accumulation affects AFTIR

## Key helper files and functions
These contain the key functions that are called by multiple scripts above
- AFIRT_calculation.R
  -	read.param.file() = read parameter file from Excel
  -	lumped.parameters.theory() = calculates AFIRT from theory
  -	lumped parameters.simulation() = calculates 
  -	simulation() - simulates model for one set of parmeters
-	ams_graphics_v2.R - Contains some useful graphics functions
  -	scale.x.log10(), scale.y.log10() - for properly labeling log scaling
-	ams_initialize_scrpt.R - To be called at the top of every Rmd file.  Initialization code and some useful constants.
-	ams_tmdd_helper.R - Could potentially be helpful for doing sensitivity analysis.
-	ivsc_2cmtc_shedct.R - ODE for Model F 

