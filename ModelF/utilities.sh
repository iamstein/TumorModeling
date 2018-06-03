#!/usr/bin/bash

# Render a Rmarkdown file into an html file
R -e "rmarkdown::render('Task11_Jumps_In_Sensitivity_Analysis.Rmd', output_file='Task11_Jumps_In_Sensitivity_Analysis.html')"
