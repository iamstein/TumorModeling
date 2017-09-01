library(RxODE)
library(dplyr)

ivsc_2cmtc_full = function() {
  model           = list()
  model$name      = as.character(sys.calls()[[sys.nframe()]])
  
  #COMPARTMENTS AND INITIAL CONDITIONS
  model$cmtname   =             c('Amt.depot','Amt.central','Amt.periph','Conc.T'                 ,'Conc.DT')            
  model$cmtshort  =             c('Ad'       ,'Ac'         ,'Ap'        ,'T'                      ,'DT'     )
  model$init      = function(p){return(c(0   , 0           , 0          ,p["ksyn"]/p["keT"]       , 0       ))}

  #PARAMEETRS IN MODEL
  model$pin       = c('F','ka','Vc','keD','k12','k21','ksyn','keT','keDT','koff','kon'); #input parameters
  model$pode      = c('F','ka','Vc','keD','k12','k21','ksyn','keT','keDT','koff','kon'); #ode   parameters
  model$repar     = function(p){return(p)} #reparameterization function

  model$rxode.str = '
     D            =  Ac/Vc;
     d/dt(Ad)     =  -ka *Ad;
     d/dt(Ac)     = F*ka *Ad - k12*Ac + k21*Ap    - keD *Ac  + (- kon*D*T + koff*DT)*Vc;
     d/dt(Ap)     =            k12*Ac - k21*Ap;
     d/dt(T)      = ksyn                          - keT * T     - kon*D*T + koff*DT;
     d/dt(DT)     =                            		- keDT*DT     + kon*D*T - koff*DT;
  '
  model$rxode     = RxODE(model = model$rxode.str, modName = model$name)
  
  model$rxout     = function(result,p) {
    result        = as.data.frame(result)
    T0            = p["ksyn"]/p["keT"]
    result$Tpchg  = result$T/T0
    result$Dtot   = result$D + result$DT
    result$Ttot   = result$T + result$DT
    return(result) }
  
  model$check     = "Task01_full_example.R"
  
  #FUNCTIONS OF KEY QUANTITIES FROM MODEL PARAMETERS
  model$fun.Cnonlin = function(p){
    Cnonlin = unname(p["ksyn"]/(p["keD"]))
    stopifnot(!is.na(Cnonlin))
    return(Cnonlin)
  }
  
  model$fun.AFI   = function(p){
    p   = as.list(p)
    AFI = (p$koff/p$kon)*(p$keT/p$keDT)*(p$CL*p$tau/p$dose)
    return(AFI) 
  }
  model$fun.TFIiv = function(p,q){
    p = as.list(p)
    q = as.list(q)
    TFI = (p$koff/p$kon)*(p$keT/p$keDT)/
      (dose*(q$C1*exp(-q$lam1*tau)/(1-exp(-q$lam1*tau))) + 
         dose*(q$C2*exp(-q$lam2*tau)/(1-exp(-q$lam2*tau))))
    return(TFI)
  }
  model$fun.TFIsc = function(p,q){error('not programmed yet')}  
  model$fun.Macro = function(p)  {error('not programmed yet')}   
  
  return(model)
}

ivsc_2cmtc_full_Keq_koff = function() {
  model           = ivsc_2cmtc_full()
  model$name      = as.character(sys.calls()[[sys.nframe()]])
  
  #PARAMEETRS IN MODEL
  model$pin       = c('F','ka','Vc','keD','k12','k21','ksyn','keT','keDT','koff','Keq'); #input parameters
  model$pode      = c('F','ka','Vc','keD','k12','k21','ksyn','keT','keDT','koff','kon'); #ode   parameters
  
  model$repar     = function(p){ #reparameterization function
    p     = as.list(p)
    p$kon = p$koff/p$Keq
    return(unlist(p))
  }  
  
  return(model)
}

ivsc_2cmtc_full_CLV_VmCafirT0TtotssKoff = function() {
  model           = ivsc_2cmtc_full()
  model$name      = as.character(sys.calls()[[sys.nframe()]])
  model$init      = function(p){return(c(0,0,0,p["T0"],0))}
  
  #PARAMEETRS IN MODEL
  model$pin       = c('F','ka','Vc','CL' ,'Q'  ,'Vp' ,'Vm'  ,'T0' ,'Ttotss','koff','Cafir'); #input parameters
  model$pode      = c('F','ka','Vc','keD','k12','k21','ksyn','keT','keDT'  ,'koff','kon'  ); #ode   parameters
  model$repar     = function(p){
    p = p %>%
      t() %>%
      as.data.frame() %>%
      mutate(
        keD = CL/Vc,
        k12 = Q /Vc,
        k21 = Q /Vp,   
        ksyn= Vm/Vc,
        keT = ksyn/T0,
        keDT= ksyn/Ttotss,
        Keq = Cafir/(Ttotss/T0),
        kon = (koff+keDT)/Keq
      ) %>%
      unlist()
    return(unlist(p))
  }
  
  model$rxout     = function(result,p) {
    T0 = p["T0"]
    result = result %>%
      as.data.frame() %>%
      mutate(Tpchg = T/T0,
             Dtot  = D + DT,
             Ttot  = T + DT)
    return(result) }
  
  return(model)
}