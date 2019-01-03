#run all key scripts for AFTIR publication

knit("Task04_Sensitivity_Analysis_Individual.Rmd")
knit("Task05e_Sensitivity_Analysis_Aggregate_wThiele.Rmd")
source("Task09d_AFIRT_vs_AFTIRsimple_wThiele.R")
source("Task09e_AFIRT_Kssd_Kss_Kd_wThiele.R")