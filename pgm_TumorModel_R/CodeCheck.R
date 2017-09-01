<<<<<<< HEAD
# For model C, if everything does what is suppose to do, then the following should happen
# 1. Set k13D = 0, then D3 should be constantly 0 through out the simulation
# 2. Set kon3 = 0, then DM3 should be constantly 0
# 3. Set kon1 = 0, then DM1 should be constantly 0
# These can be verified by sensitivity analysis

source("SensitivityAnalysis.R")
parameter_file = "../data/ivsc_3cmtct_shed3_param.csv"
model_file = "./ivsc_3cmtct_full.R"

# 1
sim = basic.sensitivity.analysis(parameter.file="/homes/li108/IMAProjs/TumorModeling/data/ivsc_3cmtct_shed3_param.csv",
                                 model = ivsc_3cmtct_full(),
                                 model.file="ivsc_3cmtct_full.R",
                                 quantity = "D3",
                                 variable= "k13D",
                                 variable.from = 0, 
                                 variable.to = 0,
                                 fold.number = 1,
                                 dose.nmol=80,
                                 tau=14,
                                 tmax = 3*28,
                                 compartment=2)
=======
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
>>>>>>> 8ad03243a2daaf1a89a375048a473adaac33ac44
print(head(sim))
