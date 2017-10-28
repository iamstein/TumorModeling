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
    p = as.data.frame(t(param.as.double))
    Kss = with(p, (koff3 + keDM3 + kshedM3)/kon3)
    Kd = with(p, koff3 / kon3)
    
    # numerators for M3.0 and Mtot3.ss(Mtot3 at steady state, M3 at initial state)
    numerator.DM   = with(p, k13DM*(VD1/VD3)*ksynM1+(keDM1+kshedDM1+k13DM)*ksynM3)
    denomenator.DM = with(p, (keDM1+kshedDM1+k13DM)*(keDM3+kshedM3+k31DM)-k31DM*k13DM)
    
    numerator.M    = with(p, k13M *(VD1/VD3)*ksynM1+(keM1 +kshedM1 +k13M) *ksynM3)
    denomenator.M  = with(p, (keM1 +kshedM1 +k13M) *(keM3 +kshedM3+k31M) -k31M *k13M)
    
    # numerator and denomenator for M3.0 ()
    Mtot3.ss = numerator.DM / denomenator.DM
    M30      = numerator.M  / denomenator.M

    # Target accumulation in the tumor compartment
    Tacc.tum = Mtot3.ss / M30

    # Biodistribution coefficient (reference: ModelF_Appendix)
    B = with(p, (k13D/(keD3 + k31D) * (VD1/VD3)))
    
    # Clearance 
    CL = with(p, (keD1*VD1))
    
    # Average drug concentration in the central compartment
    Cavg1 = dose.nmol/(CL*tau)
    
    # Average drug concentratio in the tumor compartment (I have no idea how to compute it)
    
    # AFIRT computed with Kss 
    AFIRT.theory.Kss = Kss*Tacc.tum*(CL*tau)/(dose.nmol*B)
    
    # AFIRT computed with Kd
    AFIRT.theory.Kd  = Kd*Tacc.tum*(CL*tau)/(dose.nmol*B)
  
    
    lumped_parameters_theory = data.frame(type = "theory",
                                          M30=M30, 
                                          Mtot3.ss=Mtot3.ss, 
                                          Tacc.tum=Tacc.tum,
                                          B = B,
                                          Cavg1 = Cavg1,
                                          Cavg3 = B*Cavg1,
                                          AFIRT.Kss= AFIRT.theory.Kss, #use Kss as default
                                          AFIRT.Kd = AFIRT.theory.Kd)
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
        mutate(Sfree.pct = S1/init["S1"],
             Mfree.pct = M3/init["M3"],
             dose.nmol = dose.nmol)

    # Calculate lumped parameters
    initial_state = out %>%
        filter(time==0)
    M30 = initial_state$M3
  
    ## Assume the system reaches steady state during the last dosing period
    steady_state = out %>%
    filter(time > (floor(tmax/tau)-1)*tau & time <tmax)
    Mtot3.ss = mean(steady_state$Mtot3)

    Tacc.tum = Mtot3.ss / M30
  
    ## Average drug concentration in central compartment
    dose_applied = out %>%
        filter(time > 0)
    Cavg1 = mean(dose_applied$D1)
  
  
    # Average drug concentration in tumor compartment
    Cavg3 = mean(dose_applied$D3)
  
    # AFIRT
    AFIRT = mean(steady_state$Mfree.pct)
  


    lumped_parameters_sim = data.frame(type = "simulation",
                                     M30=M30, 
                                     Mtot3.ss=Mtot3.ss, 
                                     Tacc.tum=Tacc.tum,
                                     Cavg1 = Cavg1,
                                     Cavg3 = Cavg3,
                                     B     = Cavg3/Cavg1,
                                     AFIRT = AFIRT)
  
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
    mutate(Sfree.pct = S1/init["S1"],
             Mfree.pct = M3/init["M3"],
             dose.nmol = dose.nmol)
  return(out)
}
































