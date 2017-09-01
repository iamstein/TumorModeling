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
AFIRT_theory = function(dose.nmol){
    p = as.data.frame(t(param.as.double))
    Kss = with(p, (koff3 + keDM3 + kshedM3)/kon3)
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
      
    AFIRT.theory.Kss = Kss*Tacc.tum*(CL*tau)/(dose.nmol*B)
    AFIRT.theory.Kd = Kd*Tacc.tum*(CL*tau)/(dose.nmol*B)
    return(c(AFIRT.theory.Kss, AFIRT.theory.Kd))
}

lseq = function(from, to, length.out){
    sequence = seq(log(from), log(to), length.out=length.out)
    sequence = exp(sequence)
    return(sequence)
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

sensitivity_analysis_wrt_dose.nmol = function(dose.nmol.range){
  # Return: a data frame with two columns(column1 = dose.nmol, column2 = AFIRT.sim)
  df = data.frame()
  for (dose.nmol in dose.nmol.range){
    AFIRT.sim = AFIRT_sim(dose.nmol=dose.nmol)
    AFIRT.theory = AFIRT_theory(dose.nmol = dose.nmol) 
    row = append(c(dose.nmol, AFIRT.sim), AFIRT.theory)
    df = rbind(df, row)
  }
  colnames(df) = c("dose.nmol", "AFIRT.sim", "AFIRT.theory.Kss", "AFIRT.theory.Kd")
  return(df)
}

dose.nmol.range = 80*scale.mg2nmol*lseq(0.01, 1000, 20)

df = sensitivity_analysis_wrt_dose.nmol(dose.nmol.range = dose.nmol.range)
print(df)


# sensitivity analysis on AFIRT wrt ksynM3
df = data.frame()
ksynM3.range = lseq(1, 100, 20)
ksynM3.init = param.as.double["ksynM3"]
for (i in ksynM3.range){
    ksynM3 = ksynM3.init*i
    param.as.double["ksysM3"] = ksynM3
    AFIRT.sim = AFIRT_sim(dose.nmol=0.8)
    AFIRT.theory = AFIRT_theory(dose.nmol=0.8)
    row = append(c(ksynM3, AFIRT.sim), AFIRT.theory)
    df = rbind(df, row)
}
colnames(df) = c("ksynM3", "AFIRT.sim", "AFIRT.theory.Kss", "AFIRT.theory.Kd")
print(df)



df = data.frame()
ksynM3.range = lseq(1, 100, 20)
ksynM3.init = param.as.double["VD1"]
for (i in ksynM3.range){
    VD1 = ksynM3.init*i
    param.as.double["VD1"] = ksynM3
    AFIRT.sim = AFIRT_sim(dose.nmol=0.8)
    AFIRT.theory = AFIRT_theory(dose.nmol=0.8)
    row = append(c(VD1, AFIRT.sim), AFIRT.theory)
    df = rbind(df, row)
}
colnames(df) = c("VD1", "AFIRT.sim", "AFIRT.theory.Kss", "AFIRT.theory.Kd")
print(df)


