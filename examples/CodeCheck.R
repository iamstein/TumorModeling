# Suppose the code the model F does what it supposes to to,
# then the following should happen
# 1. Set kon3 = 0, then DM3 and DS3 will be constantly 0 throughout the simulation
# 2. Set kon1 = 0, then DM1 and DS1 will be constantly 0

source("SensitivityAnalysis.R")

 
sim = basic.sensitivity.analysis(parameter.file="~/IMAProjs/TumorModeling/ModelF/ivsc_4cmtct_shedct_param.csv",
                                 model.file="~/IMAProjs/TumorModeling/ModelF/ivsc_4cmtct_shedct.R",
                                 model=ivsc_4cmtct_shedct(),
                                 quantity="DM3", 
                                 variable="kon3",
                                 variable.from=0,
                                 variable.to=0,
                                 fold.number=1,
                                 dose.nmol=80,
                                 tau=14,
                                 tmax=3*28,
                                 compartment=2)

print(dim(sim))
print(head(sim))
