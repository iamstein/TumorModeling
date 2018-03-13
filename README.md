# TumorModeling
Tumor Model IMA 2017

These contain the key functions that are called my multiple tasks below
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

## Tasks
These are short Tasks in Rmarkdown for exploring results and testing out new functionality
-	Task00 - Hongshan’s first CodeDemo of the AFIRT_calculation.R
-	Task01 - Andy getting rid of enough bugs taht the AFIRT theory matches the simulation.  Unrealistic dosing scenarios were explored
- Task03 - Implement AFIRTM and AFIRTS for membrane-bound and soluble target
-	Task02 - Andy trying out realistic doses to see how AFIRT holds up
-	Task04 - Sameed working on some sensitivity analyses of other parameters. 
- Task05 - Create the 3x6 plot for Sensitivity Analysis for Atezo, Pembro, and Herceptin.

