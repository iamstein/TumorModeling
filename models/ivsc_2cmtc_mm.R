#parameterizations include
# default (rate constants)
# CLV (clearances and volumes)

library(RxODE)

#------------------------------------------------------------------------------------
ivsc_2cmtc_mm = function() {
  model           = list()
  model$name      = as.character(sys.calls()[[sys.nframe()]])
  
  #COMPARTMENTS AND INITIAL CONDITIONS
  model$cmtname   =                    c('Amt.depot','Amt.central','Amt.periph')            
  model$cmtshort  =                    c('Ad'       ,'Ac'         ,'Ap'        )
  model$init      = function(p){return(c( 0         , 0           , 0   ))}

  #PARAMEETRS IN MODEL
  model$pin       = c('F','ka','Vc','keD','k12','k21','Vm','Km'); #input parameters
  model$pode      = c('F','ka','Vc','keD','k12','k21','Vm','Km'); #ode   parameters
  model$repar     = function(p){return(p)} #reparameterization function
  
  model$fun.Cnonlin = function(p){
    Cnonlin = with(as.list(p),(Vm/Vc)/keD)
    return(Cnonlin)
  }
    
  model$rxode.str = '
     D            =  Ac/Vc;
     d/dt(Ad)     =  -ka *Ad;
     d/dt(Ac)     = F*ka *Ad - k12*Ac + k21*Ap    - keD *Ac  - Vm*D/(D+Km);
     d/dt(Ap)     =            k12*Ac - k21*Ap;
  '
  model$rxode     = RxODE(model = model$rxode.str, modName = model$name)
  
  model$rxout     = function(result,p) {
    result        = as.data.frame(result)
    return(result)
  }
  
  model$check     = "Task02_mm_example.R"
  
  return(model)
}

#------------------------------------------------------------------------------------
ivsc_2cmtc_mm_CLV = function() {
  model       = ivsc_2cmtc_mm()
  
  model$name  = as.character(sys.calls()[[sys.nframe()]])
  model$pin   = c('F','ka','Vc','CL','Q','Vp','Vm','Km'); #input parameters 
  model$repar = function(p) {
    p     = as.list(p)
    p$keD = p$CL/p$Vc
    p$k12 = p$Q /p$Vc
    p$k21 = p$Q /p$Vp    
    return(unlist(p))
  }
  
  model$fun.Cnonlin = function(p){
    Cnonlin = with(as.list(p),Vm/CL)
    return(Cnonlin)
  }
  
  return(model)
}

#------------------------------------------------------------------------------------
ivsc_2cmtc_mm_CLV_Keq = function() {
  model       = ivsc_2cmtc_mm_CLV()
  
  model$name  = as.character(sys.calls()[[sys.nframe()]])
  model$pin   = c('F','ka','Vc','CL','Q','Vp','Vm','Keq'); #input parameters 
  model$repar = function(p) {
    p     = as.data.frame(as.list(p))
    p     = mutate(p,
                   keD = CL/Vc,
                   k12 = Q/Vc,
                   k21 = Q/Vp,
                   Km  = Keq)
    return(unlist(p))
  }

  return(model)
}

#------------------------------------------------------------------------------------
ivsc_2cmtc_mm_CLV_Cafir = function() {
  model       = ivsc_2cmtc_mm_CLV()
  
  model$name  = as.character(sys.calls()[[sys.nframe()]])
  model$pin   = c('F','ka','Vc','CL','Q','Vp','Vm','Cafir'); #input parameters 
  model$repar = function(p) {
    p     = as.data.frame(as.list(p))
    p     = mutate(p,
                   keD = CL/Vc,
                   k12 = Q/Vc,
                   k21 = Q/Vp,
                   Km  = Cafir)
    return(unlist(p))
  }
  
  return(model)
}

#------------------------------------------------------------------------------------
ivsc_2cmtc_mm_CLV_ksynKeq = function() {
  model       = ivsc_2cmtc_mm_CLV()
  
  model$name  = as.character(sys.calls()[[sys.nframe()]])
  model$pin   = c('F','ka','Vc','CL','Q','Vp','ksyn','Keq'); #input parameters 

  model$repar = function(p) {
    p     = as.data.frame(as.list(p))
    p     = mutate(p,
                   keD = CL/Vc,
                   k12 = Q/Vc,
                   k21 = Q/Vp,
                   Km  = Keq,
                   Vm  = ksyn*Vc)
    return(unlist(p))
  }  
  
  model$fun.Cnonlin = function(p) {
    Cnonlin = with(as.list(p),(ksyn*Vc)/CL)
    return(Cnonlin)
  }

  return(model)
}

#------------------------------------------------------------------------------------
ivsc_1cmt_mm_CLV_Keq = function() {
  model       = ivsc_2cmtc_mm_CLV_ksynKeq()
  
  model$name  = as.character(sys.calls()[[sys.nframe()]])
  model$pin   = c('F','ka','Vc','CL','Vm','Keq'); #input parameters 
  
  model$repar = function(p) {
    p     = as.data.frame(as.list(p))
    p     = mutate(p,
                   keD = CL/Vc,
                   k12 = 0,
                   k21 = 0,
                   Km  = Keq)
    return(unlist(p))
  }  
  
  model$fun.Cnonlin = function(p) {
    Cnonlin = with(as.list(p),Vm/CL)
    return(Cnonlin)
  }  
  
  model$check       = ""
  
  return(model)
}

#------------------------------------------------------------------------------------
ivsc_1cmt_mm_CLV_ksynKeq = function() {
  model       = ivsc_2cmtc_mm_CLV_ksynKeq()
  
  model$name  = as.character(sys.calls()[[sys.nframe()]])
  model$pin   = c('F','ka','Vc','CL','ksyn','Keq'); #input parameters 

  model$repar = function(p) {
    p     = as.data.frame(as.list(p))
    p     = mutate(p,
                   keD = CL/Vc,
                   k12 = 0,
                   k21 = 0,
                   Km  = Keq,
                   Vm  = ksyn*Vc)
    return(unlist(p))
  }  
  
  model$check       = ""
  
  return(model)
}

