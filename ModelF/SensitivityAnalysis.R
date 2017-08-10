# design a function that systematically performs sensitivity analysis

source("ams_initialize_script.R")
source("ivsc_4cmtct_shedct.R")

# Global Variables
model = ivsc_4cmtct_shedct()
tmax = 3*28 # End of the observation time, unit=day
tau = 14 # dosing interval, unit=day
compartment = 2 # compartment to which dosing is applied

# Import parameters
d = read.csv("ivsc_4cmtct_shedct_param.csv")
param.as.double = d$Value
names(param.as.double) = d$Parameter

# Theoretical AFIRT
AFIRT_theory.Kss = function(dose.nmol){
  p = as.data.frame(t(param.as.double))
  Kss = with(p, (koff3 + keDM3 + kshedM3)/kon3)
  
  # numerator and denomenator for Mtot3.ss(Mtot3 at steady state)
  numerator = with(p, k13DM*(VD1/VD3)*ksynM1+(keDM1+kshedM1+k13DM)*kshedM3)
  denomenator = with(p, (keDM1+kshedM1+k13DM)*(keDM3+kshedM3+k31DM)-k31DM*k13DM)
  Mtot3.ss = numerator / denomenator
  
  # numerator and denomenator for M3.0 (M3 at initial state)
  numerator = with(p, k13M*(VD1/VD3)*ksynM1+(keM1+kshedM1+k13D)*kshedM3)
  denomenator = with(p, (keM1+kshedM1+k13D)*(keD3+kshedM3+k31D)-k31D*k13D)
  M3.0 = numerator / denomenator
  
  # Target accumulation in the tumor compartment
  Tacc.tum = Mtot3.ss / M3.0
  
  CL = with(p, keD1 / VD1)
  B = with(p, (k13D*VD1/VD3)/(keD3 + k31D))
  
  AFIRT.theory.Kss = Kss*Tacc.tum*(CL*tau)/(dose.nmol*B)
  return(AFIRT.theory.Kss)
}

AFIRT_theory.Kss = function(dose.nmol){
  p = as.data.frame(t(param.as.double))
  Kd = with(p, koff3 / kon3)
  
  # numerator and denomenator for Mtot3.ss(Mtot3 at steady state)
  numerator = with(p, k13DM*(VD1/VD3)*ksynM1+(keDM1+kshedM1+k13DM)*kshedM3)
  denomenator = with(p, (keDM1+kshedM1+k13DM)*(keDM3+kshedM3+k31DM)-k31DM*k13DM)
  Mtot3.ss = numerator / denomenator
  
  # numerator and denomenator for M3.0 (M3 at initial state)
  numerator = with(p, k13M*(VD1/VD3)*ksynM1+(keM1+kshedM1+k13D)*kshedM3)
  denomenator = with(p, (keM1+kshedM1+k13D)*(keD3+kshedM3+k31D)-k31D*k13D)
  M3.0 = numerator / denomenator
  
  # Target accumulation in the tumor compartment
  Tacc.tum = Mtot3.ss / M3.0
  
  CL = with(p, keD1 / VD1)
  B = with(p, (k13D*VD1/VD3)/(keD3 + k31D))
  
  AFIRT.theory.Kd = Kd*Tacc.tum*(CL*tau)/(dose.nmol*B)
  return(AFIRT.theory.Kd)
}


AFIRT_sim = function(dose.nmol){
  ev = eventTable(amount.units="nmol", time.units="days")
  sample.points = c(seq(-7, tmax, 0.1), 10^(-3:0)) # sample time, increment by 0.1
  sample.points = sort(sample.points)
  sample.points = unique(sample.points)
  ev$add.sampling(sample.points)
  ev$add.dosing(dose=dose.nmol, nbr.doses=floor(tmax/tau)+1, dosing.interval=tau,
                dosing.to=compartment)
  
  init = model$init(param.as.double)
  out = model$rxode$solve(param.as.double, ev, init)
  out = model$rxout(out)
  out = out %>% 
    mutate(Sfree.pct = S1/init["S1"],
           Mfree.pct = M3/init["M3"],
           dose.nmol = dose.nmol)
  
  last.two.doses = out %>%
    filter(time > 2*28 & time < 3*28)
  
  AFIRT.sim = mean(last.two.doses$Mfree.pct)
  
  return(AFIRT.sim)
}

lseq = function(from, to, length.out){
  sequence = seq(log(from), log(to), length.out=length.out)
  sequence = exp(sequence)
  return(sequence)
}

sensitivity_analysis = function(target, variable, variable.range){
  # Arguments:
  # target: a list of quantities the sensitivity analysis is perfomed on
  # variable: a list of variable the sensitivity is performed against
  # variable.range: a list of range for each variable 
  
  # Return: 
  # a data frame with two columns(column1 = dose.nmol, column2 = AFIRT.sim)
  
  LIST = list()
  for (variable.index in 1:length(variable)){
    df = data.frame()
    for (variable.value in variable.range[variable.index]){
      row = c(variable.value)
      for (tar in target){
        tar.value = tar(var=variable.value)
        row = append(row, tar.value)
      }
      df = rbind(df, row)
    }
    LIST = append(LIST, df)
  }
  
  return(LIST) 
}  
