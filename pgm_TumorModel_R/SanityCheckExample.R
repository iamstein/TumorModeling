
source("SensitivityAnalysis.R")

sim = sanity.check(parameter.file="~/IMAProjs/TumorModeling/ModelF/ivsc_4cmtct_shedct_param.csv",
                                 model.file="~/IMAProjs/TumorModeling/ModelF/ivsc_4cmtct_shedct.R",
                                 model=ivsc_4cmtct_shedct(),
                                 quantities=list("DM1","DS1","DM3","DS3","D3", "D2"), 
                                 variables=list("kon3","kon1","k12D", "k12D"),
                                 dose.nmol=80,
                                 tau=14,
                                 tmax=3*28,
                                 compartment=2)

quantities = list("DM1","DS1", "DM3", "DS3", "D3", "D2")
for (quan in quantities){
    value = sim$quan
    for (j in value){
        if(j !=0){
            print("At time =", sim['time',j],quan,"is",value[j])
        } 
      }
}
    
    

