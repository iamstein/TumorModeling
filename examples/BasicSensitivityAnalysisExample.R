# In this example, we perform sensitivity analysis on Dtot3 against the variable kon3,
# The model is Model F
# The range of the variable is in [0.01, 1] with 10 folds.




source("SensitivityAnalysis.R")

sim = basic.sensitivity.analysis(parameter.file="~/IMAProjs/TumorModeling/ModelF/ivsc_4cmtct_shedct_param.csv",
                                 model.file="~/IMAProjs/TumorModeling/ModelF/ivsc_4cmtct_shedct.R",
                                 model=ivsc_4cmtct_shedct(),
                                 quantity="Dtot3", 
                                 variable="kon3",
                                 variable.from=0.01,
                                 variable.to=1,
                                 fold.number=10,
                                 dose.nmol=80,
                                 tau=14,
                                 tmax=3*28,
                                 compartment=2)

print(dim(sim))
print(head(sim))

