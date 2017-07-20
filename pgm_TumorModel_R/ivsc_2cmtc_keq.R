library(RxODE)

ivsc_2cmtc_keq = function() {
  model           = list()
  model$name      = as.character(sys.calls()[[sys.nframe()]])
  
  #COMPARTMENTS AND INITIAL CONDITIONS
  model$cmtname   =                    c('Amt.depot','Amt.central','Amt.periph','Conc.Ttot'         )            
  model$cmtshort  =                    c('Ad'       ,'Actot'      ,'Ap'        ,'Ttot'              )
  model$init      = function(p){return(c( 0         , 0           , 0          , unname(p["ksyn"]/p["keT"])))}

  #PARAMEETRS IN MODEL
  model$pin       = c('F','ka','Vc','keD','k12','k21','ksyn','keT','keDT','Keq'); #input parameters
  model$pode      = c('F','ka','Vc','keD','k12','k21','ksyn','keT','keDT','Keq'); #ode   parameters
  model$repar     = function(p){return(p)} #reparameterization function
  
  model$rxode.str = '
   Dtot         =  Actot/Vc;
   B            =  Dtot  - Ttot - Keq;
   D            = 0.5*( B + (B^2 + 4*Keq*Dtot)^0.5 ); 
   Ac           = D*Vc;

   DT           = Ttot*D/(Keq+D); 
   T            = Ttot - DT;     

   d/dt(Ad)     =  -ka *Ad;
   d/dt(Actot)  = F*ka *Ad - k12*Ac + k21*Ap    - keD *Ac - keDT*DT*Vc;
   d/dt(Ap)     =            k12*Ac - k21*Ap;
   d/dt(Ttot)   = ksyn                          - keT *T  - keDT*DT;
'
  model$rxode     = RxODE(model = model$rxode.str, modName = model$name)
  
  model$rxout     = function(result,p) {
    result        = as.data.frame(result)
    T0            = p["ksyn"]/p["keT"]
    result$Tpchg  = result$T/T0
    result$Dtot   = result$D + result$DT
    result$Ttot   = result$T + result$DT
    return(result)
  }
  model$check     = "not-yet"
    
  #FUNCTIONS OF KEY QUANTITIES FROM MODEL PARAMETERS
    model$fun.Cnonlin = function(p){return(unname(p["ksyn"]/p["keD"]))}
  
    model$fun.AFI   = function(p){
      p   = as.list(p)
      AFI = (p$Keq)*(p$keT/p$keDT)*(p$CL*p$tau/p$dose)
      return(AFI) 
    }
    model$fun.TFIiv = function(p,q){
      p = as.list(p)
      q = as.list(q)
      TFI = (p$Keq)*(p$keT/p$keDT)/
        (dose*(q$C1*exp(-q$lam1*tau)/(1-exp(-q$lam1*tau))) + 
           dose*(q$C2*exp(-q$lam2*tau)/(1-exp(-q$lam2*tau))))
      return(TFI) 
    }
    model$fun.TFIsc = function(p,q){error('not programmed yet')}  
    model$fun.Macro = function(p)  {error('not programmed yet')}  
  
    return(model)
}
ivsc_2cmtc_keq_CLV_Vm = function() {
  model           = ivsc_2cmtc_keq()
  model$name      = as.character(sys.calls()[[sys.nframe()]])
  model$init      = function(p){return(c(0,0,0,unname(p["Vm"]/p["Vc"]/p["keT"])))} #note here, keDT = keT
  
  model$pin       = c('F','ka','Vc','CL' ,'Q'  ,'Vp' ,'Vm'  ,'keT','keDT','Keq'); #input parameters
  model$pode      = c('F','ka','Vc','keD','k12','k21','ksyn','keT','keDT','Keq'); #ode   parameters
  model$repar     = function(p){
    p     = as.list(p)
    p$keD = p$CL/p$Vc
    p$k12 = p$Q /p$Vc
    p$k21 = p$Q /p$Vp    
    p$ksyn= p$Vm/p$Vc
    return(unlist(p))
  } 
  
  model$rxout     = function(result,p) {
    result        = as.data.frame(result)
    ksyn          = p["Vm"]/p["Vc"]
    T0            = ksyn/p["Vc"]
    result$Tpchg  = result$T/T0
    result$Dtot   = result$D + result$DT
    result$Ttot   = result$T + result$DT
    return(result)
  }  
  
  model$check     = "none"
  return(model)
} 


ivsc_2cmtc_keq_CLV_VmCafirT0Acc = function() {
  model           = ivsc_2cmtc_keq()
  model$name      = as.character(sys.calls()[[sys.nframe()]])
  model$init      = function(p){return(c(0,0,0,unname(p["T0"])))} #note here, keDT = keT
  
  model$pin       = c('F','ka','Vc','CL' ,'Q'  ,'Vp' ,'Vm'  ,'T0' ,'Acc'   ,'Cafir'); #input parameters
  model$pode      = c('F','ka','Vc','keD','k12','k21','ksyn','keT','keDT'  ,'Keq'); #ode   parameters
  model$repar     = function(p){
    p     = as.list(p)
    p$keD = p$CL/p$Vc
    p$k12 = p$Q /p$Vc
    p$k21 = p$Q /p$Vp    
    p$ksyn= p$Vm/p$Vc
    p$keT = p$ksyn/p$T0
    p$keDT= p$ksyn/(p$Acc*p$T0)
    p$Keq = p$Cafir/p$Acc
    return(unlist(p))
  } 
  
  model$rxout     = function(result,p) {
    result        = as.data.frame(result)
    T0            = p["T0"]
    result$Tpchg  = result$T/T0
    result$Dtot   = result$D + result$DT
    result$Ttot   = result$T + result$DT
    return(result)
  }  
  
  model$check     = "none"
  return(model)
}  

ivsc_2cmtc_keq_CLV_VmCafirT0Keq = function() {
  model           = ivsc_2cmtc_keq_CLV_VmCafirT0Acc()
  model$name      = as.character(sys.calls()[[sys.nframe()]])

  model$pin       = c('F','ka','Vc','CL' ,'Q'  ,'Vp' ,'Vm'  ,'T0' ,'Cafir' ,'Keq'); #input parameters
  model$pode      = c('F','ka','Vc','keD','k12','k21','ksyn','keT','keDT'  ,'Keq'); #ode   parameters
  model$repar     = function(p){
    p     = as.list(p)
    p$keD = p$CL/p$Vc
    p$k12 = p$Q /p$Vc
    p$k21 = p$Q /p$Vp    
    p$ksyn= p$Vm/p$Vc
    p$keT = p$ksyn/p$T0
    Acc   = p$Cafir/p$Keq
    p$keDT= p$ksyn/(Acc*p$T0)
    p$Keq = p$Cafir/Acc
    return(unlist(p))
  } 
  return(model)
}  

ivsc_2cmtc_keq_CLV_VmCafirT0Ttotss = function() {
  model           = ivsc_2cmtc_keq_CLV_VmCafirT0Acc()
  model$name      = as.character(sys.calls()[[sys.nframe()]])
  
  model$pin       = c('F','ka','Vc','CL' ,'Q'  ,'Vp' ,'Vm'  ,'T0' ,'Ttotss','Cafir'); #input parameters
  model$pode      = c('F','ka','Vc','keD','k12','k21','ksyn','keT','keDT'  ,'Keq'); #ode   parameters
  model$repar     = function(p){
    p     = as.list(p)
    p$keD = p$CL/p$Vc
    p$k12 = p$Q /p$Vc
    p$k21 = p$Q /p$Vp    
    p$ksyn= p$Vm/p$Vc
    p$keT = p$ksyn/p$T0
    p$keDT= p$ksyn/p$Ttotss
    Acc   = p$Ttotss/p$T0
    p$Keq = p$Cafir/Acc
    return(unlist(p))
  } 
  return(model)
} 

ivsc_2cmtc_keq_CLV_VmT0Ttotss = function() {
  model           = ivsc_2cmtc_keq()
  model$name      = as.character(sys.calls()[[sys.nframe()]])
  model$init      = function(p){return(c(0,0,0,unname(p["T0"])))} #note here, keDT = keT
  
  model$pin       = c('F','ka','Vc','CL' ,'Q'  ,'Vp' ,'Vm'  ,'T0' ,'Ttotss','Keq'); #input parameters
  model$pode      = c('F','ka','Vc','keD','k12','k21','ksyn','keT','keDT'  ,'Keq'); #ode   parameters
  model$repar     = function(p){
    p     = as.list(p)
    p$keD = p$CL/p$Vc
    p$k12 = p$Q /p$Vc
    p$k21 = p$Q /p$Vp    
    p$ksyn= p$Vm/p$Vc
    p$keT = p$ksyn/p$T0
    p$keDT= p$ksyn/p$Ttotss
    return(unlist(p))
  } 
  
  model$rxout     = function(result,p) {
    result        = as.data.frame(result)
    T0            = p["T0"]
    result$Tpchg  = result$T/T0
    result$Dtot   = result$D + result$DT
    result$Ttot   = result$T + result$DT
    return(result)
  }  
  
  model$check     = "none"
  return(model)
}  

ivsc_2cmtc_keq_CLV_VmT0Tmult = function() {
  model           = ivsc_2cmtc_keq_CLV_VmT0Ttotss()
  model$name      = as.character(sys.calls()[[sys.nframe()]])

  model$pin       = c('F','ka','Vc','CL' ,'Q'  ,'Vp' ,'Vm'  ,'T0' ,'Tmult' ,'Keq'); #input parameters
  model$pode      = c('F','ka','Vc','keD','k12','k21','ksyn','keT','keDT'  ,'Keq'); #ode   parameters
  
  model$repar     = function(p){
    p     = as.list(p)
    p$keD = p$CL/p$Vc
    p$k12 = p$Q /p$Vc
    p$k21 = p$Q /p$Vp    
    p$ksyn= p$Vm/p$Vc
    p$keT = p$ksyn/p$T0
    Ttotss= p$T0*p$Tmult
    p$keDT= p$ksyn/Ttotss
    return(unlist(p))
  } 
  
  return(model)
}

ivsc_2cmtc_keq_CLV_VmTtotssTmult = function() {
  #T0 = Ttotss*Tmult
  model           = ivsc_2cmtc_keq_CLV_VmT0Ttotss()
  model$name      = as.character(sys.calls()[[sys.nframe()]])
  model$init      = function(p){return(c(0,0,0,unname(p["Ttotss"]*p["Tmult"])))} #note here, keDT = keT
  
  model$pin       = c('F','ka','Vc','CL' ,'Q'  ,'Vp' ,'Vm'  ,'Tmult','Ttotss' ,'Keq'); #input parameters
  model$pode      = c('F','ka','Vc','keD','k12','k21','ksyn','keT'  ,'keDT'   ,'Keq'); #ode   parameters
  
  model$repar     = function(p){
    p     = as.list(p)
    p$keD = p$CL/p$Vc
    p$k12 = p$Q /p$Vc
    p$k21 = p$Q /p$Vp    
    p$ksyn= p$Vm/p$Vc
    p$keDT= p$ksyn/p$Ttotss
    T0    = p$Ttotss*p$Tmult
    p$keT = p$ksyn/T0
    return(unlist(p))
  } 
  
  model$rxout     = function(result,p) {
    result        = as.data.frame(result)
    T0            = p["Ttotss"]*p["Tmult"]
    result$Tpchg  = result$T/T0
    result$Dtot   = result$D + result$DT
    result$Ttot   = result$T + result$DT
    return(result)
  }    
  
  return(model)
}

ivsc_2cmtc_keq_CLV_VmTtotssKeqtmu = function() {
  #T0 = Ttotss*Keq/Keqtmu;  in other words, Keqtmu = Keq/Tmult and T0=Ttotss*Tmult
  model           = ivsc_2cmtc_keq_CLV_VmT0Ttotss()
  model$name      = as.character(sys.calls()[[sys.nframe()]])
  model$init      = function(p){return(c(0,0,0,unname(p["Ttotss"]*p["Keq"]/p["Keqtmu"])))} #note here, keDT = keT
  
  model$pin       = c('F','ka','Vc','CL' ,'Q'  ,'Vp' ,'Vm'  ,'Keqtmu','Ttotss' ,'Keq'); #input parameters
  model$pode      = c('F','ka','Vc','keD','k12','k21','ksyn','keT'   ,'keDT'   ,'Keq'); #ode   parameters
  
  model$repar     = function(p){
    p     = as.list(p)
    p$keD = p$CL/p$Vc
    p$k12 = p$Q /p$Vc
    p$k21 = p$Q /p$Vp    
    p$ksyn= p$Vm/p$Vc
    p$keDT= p$ksyn/p$Ttotss
    T0    = p$Ttotss*p$Keq/p$Keqtmu
    p$keT = p$ksyn/T0
    return(unlist(p))
  } 
  
  model$rxout     = function(result,p) {
    result        = as.data.frame(result)
    T0            = p["Ttotss"]*p["Keq"]/p["Keqtmu"]
    result$Tpchg  = result$T/T0
    result$Dtot   = result$D + result$DT
    result$Ttot   = result$T + result$DT
    return(result)
  }    
  
  return(model)
}


ivsc_2cmtc_keqct = function() {
  model           = ivsc_2cmtc_keq()
  model$name      = as.character(sys.calls()[[sys.nframe()]])
  model$init      = function(p){return(c(0,0,0,unname(p["ksyn"]/p["keDT"])))} #note here, keDT = keT
  
  model$pin       = c('F','ka','Vc','CL' ,'Q'  ,'Vp' ,'ksyn'      ,'keDT','Keq'); #input parameters
  model$pode      = c('F','ka','Vc','keD','k12','k21','ksyn','keT','keDT','Keq'); #ode   parameters
  model$repar     = function(p){
    p["keT"] = p["keDT"]
    return(p)
  }
  
  model$rxout     = function(result,p) {
    result        = as.data.frame(result)
    T0            = p["ksyn"]/p["keDT"] #note here, keDT = keT
    result$Tpchg  = result$T/T0
    result$Dtot   = result$D + result$DT
    result$Ttot   = result$T + result$DT
    return(result)
  }  
  
  model$check     = "Task04_keqct_explore.R"
  return(model)
}  

ivsc_2cmtc_keqct_CLV = function() {
  model           = ivsc_2cmtc_keqct()
  model$name      = as.character(sys.calls()[[sys.nframe()]])
  model$pin       = c('F','ka','Vc','CL' ,'Q'  ,'Vp','ksyn'      ,'keDT','Keq'); #input parameters
  model$pode      = c('F','ka','Vc','keD','k12','k21','ksyn','keT','keDT','Keq'); #ode   parameters
  model$repar     = function(p){
    p     = as.list(p)
    p$keD = p$CL/p$Vc
    p$k12 = p$Q /p$Vc
    p$k21 = p$Q /p$Vp    
    p$keT = p$keDT
    return(unlist(p))
  }  
  return(model)
}

ivsc_2cmtc_keqct_CLV_T0 = function() {
  model           = ivsc_2cmtc_keqct()
  model$name      = as.character(sys.calls()[[sys.nframe()]])
  model$init      = function(p){return(c(0,0,0,p["T0"]))} #note here, keDT = keT
  model$pin       = c('F','ka','Vc','CL' ,'Q'  ,'Vp','ksyn'      ,'T0'  ,'Keq'); #input parameters
  model$pode      = c('F','ka','Vc','keD','k12','k21','ksyn','keT','keDT','Keq'); #ode   parameters
  model$repar     = function(p){
    p     = as.list(p)
    p$keD = p$CL/p$Vc
    p$k12 = p$Q /p$Vc
    p$k21 = p$Q /p$Vp    
    p$keT = p$ksyn/p$T0
    p$keDT= p$keT
    return(unlist(p))
  }  
  model$rxout     = function(result,p) {
    result        = as.data.frame(result)
    result$Tpchg  = result$T/p["T0"]
    result$Dtot   = result$D + result$DT
    result$Ttot   = result$T + result$DT
    return(result)
  }    
  return(model)
}

ivsc_2cmtc_keqct_CLV_VmT0 = function() {
  model           = ivsc_2cmtc_keqct()
  model$name      = as.character(sys.calls()[[sys.nframe()]])
  model$init      = function(p){return(c(0,0,0,p["T0"]))} #note here, keDT = keT
  model$pin       = c('F','ka','Vc','CL' ,'Q'  ,'Vp' ,'Vm'        ,'T0'  ,'Keq'); #input parameters
  model$pode      = c('F','ka','Vc','keD','k12','k21','ksyn','keT','keDT','Keq'); #ode   parameters
  model$repar     = function(p){
    p     = as.list(p)
    p$keD = p$CL/p$Vc
    p$k12 = p$Q /p$Vc
    p$k21 = p$Q /p$Vp    
    p$ksyn= p$Vm/p$Vc
    p$keT = p$ksyn/p$T0
    p$keDT= p$keT
    return(unlist(p))
  }  
  model$rxout     = function(result,p) {
    result        = as.data.frame(result)
    result$Tpchg  = result$T/p["T0"]
    result$Dtot   = result$D + result$DT
    result$Ttot   = result$T + result$DT
    return(result)
  }    
  return(model)
}

ivsc_2cmtc_keqct_CLV_VmCafirT0 = function() {
  model           = ivsc_2cmtc_keqct_CLV_VmT0()
  model$name      = as.character(sys.calls()[[sys.nframe()]])
  
  model$pin       = c('F','ka','Vc','CL' ,'Q'  ,'Vp' ,'Vm'        ,'T0'  ,'Cafir'); #input parameters
  model$pode      = c('F','ka','Vc','keD','k12','k21','ksyn','keT','keDT','Keq'  ); #ode   parameters
  model$repar     = function(p){
    p     = as.list(p)
    p$keD = p$CL/p$Vc
    p$k12 = p$Q /p$Vc
    p$k21 = p$Q /p$Vp    
    p$ksyn= p$Vm/p$Vc
    p$keT = p$ksyn/p$T0
    p$keDT= p$keT
    p$Keq = p$Cafir
    return(unlist(p))
  }  
  
  return(model)
}

ivsc_2cmtc_keqct_CLV_VmCafirTtotss = function() {
  model           = ivsc_2cmtc_keqct_CLV_VmCafirT0()
  model$name      = as.character(sys.calls()[[sys.nframe()]])
  model$init      = function(p){return(c(0,0,0,p["Ttotss"]))} #note here, keDT = keT
  
  model$pin       = c('F','ka','Vc','CL' ,'Q'  ,'Vp' ,'Vm'        ,'Ttotss','Cafir'); #input parameters
  model$pode      = c('F','ka','Vc','keD','k12','k21','ksyn','keT','keDT'  ,'Keq'  ); #ode   parameters
  model$repar     = function(p){
    p     = as.list(p)
    p$keD = p$CL/p$Vc
    p$k12 = p$Q /p$Vc
    p$k21 = p$Q /p$Vp    
    p$ksyn= p$Vm/p$Vc
    p$keDT= p$ksyn/p$Ttotss
    p$keT = p$keDT
    p$Keq = p$Cafir
    return(unlist(p))
  }  
  
  model$rxout     = function(result,p) {
    result        = as.data.frame(result)
    T0 = p["Ttotss"]
    result$Tpchg  = result$T/T0
    result$Dtot   = result$D + result$DT
    result$Ttot   = result$T + result$DT
    return(result)
  } 
  return(model)
}

ivsc_2cmtc_keqct_CLV_VmTtotss = function() {
  model           = ivsc_2cmtc_keqct_CLV_VmT0()
  model$name      = as.character(sys.calls()[[sys.nframe()]])
  model$init      = function(p){return(c(0,0,0,p["Ttotss"]))} #note here, keDT = keT
  model$pin       = c('F','ka','Vc','CL' ,'Q'  ,'Vp' ,'Vm'        ,'Ttotss','Keq'); #input parameters
  model$pode      = c('F','ka','Vc','keD','k12','k21','ksyn','keT','keDT'  ,'Keq'); #ode   parameters
  model$repar     = function(p){
    p     = as.list(p)
    p$keD = p$CL/p$Vc
    p$k12 = p$Q /p$Vc
    p$k21 = p$Q /p$Vp    
    p$ksyn= p$Vm/p$Vc
    T0    = p$Ttotss
    p$keT = p$ksyn/T0
    p$keDT= p$keT
    return(unlist(p))
  }  
  model$rxout     = function(result,p) {
    result        = as.data.frame(result)
    T0            = p["Ttotss"]
    result$Tpchg  = result$T/T0
    result$Dtot   = result$D + result$DT
    result$Ttot   = result$T + result$DT
    return(result)
  }    
  return(model)
}

#wagner model, described in Yan12twodrugs
ivsc_2cmtc_keqctwag_CLV_T0 = function() {
  model           = ivsc_2cmtc_keqct_CLV_T0()
  model$name      = as.character(sys.calls()[[sys.nframe()]])
  
  model$cmtname   =                    c('Amt.depot','Amt.central','Amt.periph')            
  model$cmtshort  =                    c('Ad'       ,'Ac'         ,'Ap'        )
  model$init      = function(p){return(c( 0         , 0           , 0          ))}
  
  model$rxode.str = '
    D            = Ac/Vc;
    denom        = 1 + T0*Keq/(Keq+D)^2

    DT           = T0*D/(Keq+D); 
    T            = T0 - DT;     
    Ttot         = T0;
    
    d/dt(Ad)     =  -ka *Ad;
    d/dt(Ac)     = F*ka *Ad - k12*Ac + k21*Ap    - (keD *D - keDT*DT)/denom*Vc;
    d/dt(Ap)     =            k12*Ac - k21*Ap;
  '
  
  model$rxode     = RxODE(model = model$rxode.str, modName = model$name)

  return(model)
}