# Helper function that returns a range of variable when performing 
# sensitvity analysis

read.param.file = function(filename) {
  d = read_excel(filename, 1)
  param.as.double        = d$Value
  names(param.as.double) = d$Parameter
  param.as.double        = param.as.double[model$pin] #keep only parameters used in ODE
}

lseq = function(from, to, length.out){
    # Arguments:
    #   from : initial value of the variable
    #   to : teminal value of the variable
    #   length.out : fold number of <to - from>
    # Return :
    #   A vector of length <length.out> 

    sequence = seq(log(from), log(to), length.out=length.out)
    sequence = exp(sequence)
    return(sequence)
}


# Function computing lumped parameters from theory

lumped.parameters.theory = function(param.as.double=param.as.double, 
                                    dose.nmol=dose.nmol, 
                                    tau=tau){
    # Arguments:
    #   params_file_path: full path of the parameters file.
    #   dose.nmol: dosing amout in nmol
    #   tau: dosing interval in days
    # Return:
    #   A data frame of lumped parameters calculated from theory


    #d <- xlsx::read.xlsx(params_file_path, 1)
    #param.as.double = d$Value
    #names(param.as.double) = d$Parameter
    p    = as.data.frame(t(param.as.double))

    # transmission coefficients for membrane targets
    Kssd.M = with(p, (koff3 + keDM3 + kshedDM3 + k31DM)/kon3)
    Kss.M  = with(p, (koff3 + keDM3 + kshedDM3)        /kon3)
    Kd.M   = with(p,  koff3                            /kon3)

    # transmission coefficients for soluble targets
    Kssd.S = with(p, (koff3 + keDS3 + kshedDM3 + k31DS)/kon3)
    Kss.S  = with(p, (koff3 + keDS3 + kshedDM3)        /kon3)
    Kd.S   = with(p,  koff3                            /kon3)
    
    # numerators and denominator for membrane target at intial state and steady state
    numerator.M1 = with(p, (kshedM3+k31M + keM3)*ksynM1 + (VD3/VD1)*k31M*ksynM3)
    denomenator.M1  = with(p, (keM1 +kshedM1 +k13M) *(keM3 +kshedM3+k31M) -k31M *k13M)

    numerator.DM1 = with(p, (kshedDM3 + k31DM + keDM3)*ksynM1 + (VD3/VD1)*k31DM*ksynM3)
    denomenator.DM1 = with(p, (keDM1+kshedDM1+k13DM)*(keDM3+kshedM3+k31DM)-k31DM*k13DM)

    numerator.M3    = with(p, k13M *(VD1/VD3)*ksynM1+(keM1 +kshedM1 +k13M) *ksynM3)
    denomenator.M3  = with(p, (keM1 +kshedM1 +k13M) *(keM3 +kshedM3+k31M) -k31M *k13M)

    numerator.DM3   = with(p, k13DM*(VD1/VD3)*ksynM1+(keDM1+kshedDM1+k13DM)*ksynM3)
    denomenator.DM3 = with(p, (keDM1+kshedDM1+k13DM)*(keDM3+kshedM3+k31DM)-k31DM*k13DM)
    
    # compute membrane target at initial and steady state
    M10 = numerator.M1 / denomenator.M1
    Mtot1.ss = numerator.DM1 / denomenator.DM1

    M30  = numerator.M3  / denomenator.M3
    Mtot3.ss = numerator.DM3 / denomenator.DM3


    # numerators and denomenators for soluble target
    numerator.S1 = with(p, (k31S+keS3)*(ksynS1+kshedM1*M10)+(VD3/VD1)*k31S*(ksynS3+kshedM3*M30))
    denomenator.S1 = with(p, (k13S+keS1)*(k31S+keS3)-k13S*k31S)

    numerator.DS1 = with(p, (k31DS+keDS3)*(ksynS1+kshedDM1*Mtot3.ss)+(VD3/VD1)*k31DS*(ksynS3+kshedDM3*Mtot3.ss))
    denomenator.DS1 = with(p, (k13DS+keDS1)*(k31DS+keDS3)-k13DS*k31DS)

    numerator.S3 = with(p, (VD1/VD3)*k13S*(ksynS1+kshedM1*M10)+(k13S+keS1)*(ksynS3+kshedM3*M30))
    denomenator.S3 = with(p, (k13S+keS1)*(k31S+keS3)-k13S*k31S)

    numerator.DS3 = with(p, (VD1/VD3)*k13DS*(ksynS1+kshedDM1*Mtot1.ss)+(k13DS+keDS1)*(ksynS3+kshedDM3*Mtot3.ss))
    denomenator.DS3 = with(p, (k13DS+keDS1)*(k31DS+keDS3)-k13DS*k31DS)

    # compute soluble target at initial and steady state
    S10 = numerator.S1 / denomenator.S1
    Stot1.ss = numerator.DS1 / denomenator.DS3

    S30 = numerator.S3 / denomenator.S3
    Stot3.ss = numerator.DS3 / denomenator.DS3


    # Target accumulation in the tumor compartment for membrane target
    Tacc.tum.M = Mtot3.ss / M30

    # Target accumulation in the tumor compartment for soluble target
    Tacc.tum.S = Stot3.ss / S30

    # Biodistribution coefficient (reference: ModelF_Appendix)
    B = with(p, (k13D/(keD3 + k31D) * (VD1/VD3)))
    
    # Clearance 
    CL = with(p, (keD1*VD1))
    
    # Average drug concentration in the central compartment
    Cavg1 = dose.nmol/(CL*tau)
    
    # Average drug concentratio in the tumor compartment (I have no idea how to compute it)
    
    # Compute various AFIRTs for membrane target
    AFIRT.M.Kssd = Kssd.M*Tacc.tum.M/(B*Cavg1)
    AFIRT.M.Kss  = Kss.M *Tacc.tum.M/(B*Cavg1)
    AFIRT.M.Kd   = Kd.M  *Tacc.tum.M/(B*Cavg1)

    # Compute various AFIRTS for soluble target
    AFIRT.S.Kssd = Kssd.S*Tacc.tum.S/(B*Cavg1)
    AFIRT.S.Kss  = Kss.S *Tacc.tum.S/(B*Cavg1)
    AFIRT.S.Kd   = Kd.S  *Tacc.tum.S/(B*Cavg1)

    
  
    
    lumped_parameters_theory = data.frame(type = "theory",
                                          M30=M30, 
                                          Mtot3.ss=Mtot3.ss, 
                                          Tacc.tum=Tacc.tum,
                                          B = B,
                                          Cavg1 = Cavg1,
                                          Cavg3 = B*Cavg1,
                                          AFIRT.M.Kssd = AFIRT.M.Kssd,
                                          AFIRT.M.Kss  = AFIRT.M.Kss, 
                                          AFIRT.M.Kd   = AFIRT.M.Kd,
                                          AFIRT.S.Kssd = AFIRT.S.Kssd,
                                          AFIRT.S.Kss  = AFIRT.S.Kss,
                                          AFIRT.S.Kd   = AFIRT.S.Kd)
    return(lumped_parameters_theory) 
 }


# Function simulates the lumped parameters

lumped.parameters.simulation = function(model=model, param.as.double=param.as.double, 
                                        dose.nmol=dose.nmol, tmax=tmax, tau=tau, compartment){
    
    # Arguments:
    #   model_name: name of the model
    #   params_file_path: full path of the parameters file.
    #   dose.nmol: dosing amount in nmol
    #   tmax: maximum doing period in days
    #   tau: dosing interval in days
    #   compartment: compartment to which dosing is applied
    #   (in model F case, compartment=2)
    # Return:
    #   A data frame of lumped parameters calculated from simulation

    # Run simulation
    #d <- xlsx::read.xlsx(params_file_path, 1)
    #param.as.double = d$Value
    #names(param.as.double) = d$Parameter
    ev = eventTable(amount.units="nmol", time.units="days")
    sample.points = c(seq(-7, tmax, 0.1), 10^(-3:0)) # sample time, increment by 0.1
    sample.points = sort(sample.points)
    sample.points = unique(sample.points)
    ev$add.sampling(sample.points)
    ev$add.dosing(dose=dose.nmol, nbr.doses=floor(tmax/tau)+1, dosing.interval=tau,
                dosing.to=compartment)
  
    #model = do.call(model_name, list()) # model file can contain only one model
    init = model$init(param.as.double)
    out = model$rxode$solve(param.as.double, ev, init)
    out = model$rxout(out)
    out = out %>% 
        mutate(S1free.pct = S1/init["S1"],
               S3free.pct = S3/init["S3"],
               Mfree.pct = M3/init["M3"],
               dose.nmol = dose.nmol)

    # Calculate lumped parameters
    initial_state = out %>%
        filter(time==0)
    M30 = initial_state$M3
    S30 = initial_state$S3
  
    ## Assume the system reaches steady state during the last dosing period
    steady_state = out %>%
    filter(time > (floor(tmax/tau)-1)*tau & time <tmax)
    
    # Total membrane and soluble target
    Mtot3.ss = mean(steady_state$Mtot3)
    Stot3.ss = mean(steady_state$Stot3)

    # Target accumulation rate for membrane target and soluble target
    Tacc.tum.M = Mtot3.ss / M30
    Tacc.tum.S = Stot3.ss / S30
    
  
    ## Average drug concentration in central compartment
    dose_applied = out %>%
        filter(time > 0)
    Cavg1 = mean(dose_applied$D1)
 
  
    # Average drug concentration in tumor compartment
    Cavg3 = mean(dose_applied$D3)
  
    # AFIRT for the membrane target
    AFIRT.M = mean(steady_state$Mfree.pct)

    # AFIRT for the soluble target
    AFIRT.S = mean(steady_state$S3free.pct)

    lumped_parameters_sim = data.frame(type = "simulation",
                                     M30 = M30, 
                                     Mtot3.ss = Mtot3.ss,
                                     S30 = S30,
                                     Stot3.ss = Stot3.ss, 
                                     Tacc.tum.M = Tacc.tum.M,
                                     Tacc.tum.S = Tacc.tum.S,
                                     Cavg1 = Cavg1,
                                     Cavg3 = Cavg3,
                                     B     = Cavg3/Cavg1,
                                     AFIRT.M = AFIRT.M,
                                     AFIRT.S = AFIRT.S,
                                     AFIRT.M.sim = AFIRT.M, #having one named sim will be helpful later on in Task01, Task02, etc.
                                     AFIRT.S.sim = AFIRT.S)
    return(lumped_parameters_sim)
}    


simulation = function(model=model, param.as.double=param.as.double, 
                      dose.nmol=dose.nmol, tmax=tmax, tau=tau){
  #d <- xlsx::read.xlsx(params_file_path, 1)
  #param.as.double = d$Value
  #names(param.as.double) = d$Parameter
  ev = eventTable(amount.units="nmol", time.units="days")
  sample.points = c(seq(-7, tmax, 0.1), 10^(-3:0)) # sample time, increment by 0.1
  sample.points = sort(sample.points)
  sample.points = unique(sample.points)
  ev$add.sampling(sample.points)
  ev$add.dosing(dose=dose.nmol, nbr.doses=floor(tmax/tau)+1, dosing.interval=tau,
                dosing.to=2)
  
  init = model$init(param.as.double)
  out = model$rxode$solve(param.as.double, ev, init)
  out = model$rxout(out)
  out = out %>% 
    mutate(S1free.pct = S1/init["S1"],
           S3free.pct = S3/init["S3"],
           Mfree.pct = M3/init["M3"],
           dose.nmol = dose.nmol)
  return(out)
}

compare.theory.sim = function(model=model, param.as.double=param.as.double,
                          dose.nmol.range = dose.nmol.range,
                          tmax=tmax,tau=tau,dt=dt,compartment=compartment) {

  df_sim = data.frame() # put all simulations for different dose into one data frame
  for (dose.nmol in dose.nmol.range){
    row = lumped.parameters.simulation(model,
      param.as.double, dose.nmol, tmax, tau, compartment)
    df_sim = rbind(df_sim, row)
  }
  df_sim$Dose = dose.nmol.range
  
  df_thy = data.frame() # put all theoretical calculations of lumped parameters at different dose together
  for (dose.nmol in dose.nmol.range){
    row = lumped.parameters.theory(param.as.double=param.as.double, 
                                   dose.nmol=dose.nmol, 
                                   tau=tau)
    df_thy = rbind(df_thy, row)
  }
  df_thy = df_thy %>%
    mutate(Dose = dose.nmol.range,
           AFIRT = AFIRT.M.Kssd)
  
  df_compare = bind_rows(df_thy,df_sim)
  df_compare = df_compare %>%
    arrange(Dose,type) %>%
    mutate_if(is.numeric,signif,2)
    

  #do a single simulation to make sure model is at steady state
  out = simulation(model=model, param.as.double=param.as.double, 
        dose.nmol=max(dose.nmol.range), tmax=tmax, tau=tau) %>%
    select(time,D1,D3,M3,DM3,Mfree.pct) %>%
    gather(cmt,value,-c(time))
  
  g = ggplot(out,aes(x=time,y=value,group=cmt,color=cmt))
  g = g + geom_line()
  g = g + scale.y.log10()
  print(g)
  
  #plot theory vs simulation
  df_plot = df_compare %>%
    select(-contains("AFIRT"),AFIRT) %>%
    gather(param,value,-c(type,Dose))
  
  
  g = ggplot(df_plot, aes(Dose, value, group=type, color=type,linetype=type))
  g = g + scale.x.log10()
  g = g + scale.y.log10()
  g = g + geom_line()
  g = g + geom_point()
  g = g + facet_wrap(~param,scales = "free_y")
  print(g)

  return(df_compare)
}
