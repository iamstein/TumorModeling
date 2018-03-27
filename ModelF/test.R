suppressMessages(source("ams_initialize_script.R"))

model = ivsc_4cmtct_shedct()

tmax = 26*7 #days
tau  = 21   #days
compartment = 2


param.as.double =  read.param.file("../data/ModelF_Herceptin_Params.xlsx")
dose.nmol.range = lseq(.1,100,7)*scale.mpk2nmol
dose.nmol = dose.nmol.range[1]

a = lumped.parameters.theory(param.as.double = param.as.double, 
                            dose.nmol = dose.nmol, tau=tau, soluble=FALSE)

print(a)

