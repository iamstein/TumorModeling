# Helper function that returns a range of variable when performing
# sensitvity analysis

read.param.file = function(filename) {
  d = read_excel(filename, 1)
  param.as.double        = as.numeric(d$Value)
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
                                    tau=tau,
                                    soluble = FALSE){
    # Arguments:
    #   params_file_path: full path of the parameters file.
    #   dose.nmol: dosing amout in nmol
    #   tau: dosing interval in days
    #   soluble: flag saying whether or not drug is soluble. Default to FALSE.
    # Return:
    #   A data frame of lumped parameters calculated from theory

    p    = as.data.frame(t(param.as.double))

    # Calculate Mtot3.ss, M30, Stot3.ss and S30
    numerator.DM3   = with(p, k13DM*(VD1/VD3)* ksynM1                      + (keDM1 + kshedDM1 + k13DM)* ksynM3)
    numerator.M3    = with(p, k13M *(VD1/VD3)* ksynM1                      + (keM1  + kshedM1  + k13M )* ksynM3)
    numerator.DM1   = with(p, k31DM*(VD3/VD1)* ksynM3                      + (keDM3 + kshedDM3 + k31DM)* ksynM1)
    numerator.M1    = with(p, k31M *(VD3/VD1)* ksynM3                      + (keM3  + kshedM3  + k31M )* ksynM1)

    denomenator.DM_1or3 = with(p, (keDM1 + kshedDM1 + k13DM)*(keDM3 + kshedDM3 + k31DM)-k31DM*k13DM)
    denomenator.M_1or3  = with(p, (keM1  + kshedM1  + k13M )*(keM3  + kshedM3  + k31M )-k31M *k13M )
    
    Mtot1.ss = numerator.DM1 / denomenator.DM_1or3
    Mtot3.ss = numerator.DM3 / denomenator.DM_1or3
    M10      = numerator.M1  / denomenator.M_1or3
    M30      = numerator.M3  / denomenator.M_1or3   
    
    #note that this aligns with the numerator columns above and can be copied and pasted for comparison
    numerator.DS3  = with(p, k13DS*(VD1/VD3)*(ksynS1 + kshedDM1*Mtot1.ss) + (keDS1             + k13DS)*(ksynS3 + kshedDM3*Mtot3.ss))
    numerator.S3   = with(p, k13S *(VD1/VD3)*(ksynS1 + kshedM1*M10)       + (keS1              + k13S) *(ksynS3 + kshedM3 *M30))
    
    denomenator.DS3= with(p, (keDS1            + k13DS)*(keDS3            + k31DS)-k31DS*k13DS)
    denomenator.S3 = with(p, (keS1             + k13S )*(keS3             + k31S )-k31S *k13S )
    
    Stot3.ss       = numerator.DS3 / denomenator.DS3
    S30            = numerator.S3  / denomenator.S3


    if (!soluble){
      Kssd = with(p, (koff3 + keDM3 + kshedDM3 + k31DM)/kon3)
      Kss  = with(p, (koff3 + keDM3 + kshedDM3        )/kon3)
      Kd   = with(p,  koff3                            /kon3)
      Tacc.tum = Mtot3.ss / M30
      
    } else {
      Kssd = with(p, (koff3 + keDS3 +          + k31DS)/kon3)
      Kss  = with(p, (koff3 + keDS3                   )/kon3)
      Kd   = with(p,  koff3                            /kon3)
      Tacc.tum = Stot3.ss / S30
    }

    # Biodistribution coefficient (reference: ModelF_Appendix)
    B = with(p, (k13D/(keD3 + k31D) * (VD1/VD3)))

    # Clearance
    CL = with(p, (keD1*VD1))

    # Average drug concentration in the central compartment
    Cavg1 = dose.nmol/(CL*tau)

    # Compute various AFIRTs
    AFIRT.Kssd = Kssd*Tacc.tum/(B*Cavg1)
    AFIRT.Kss  = Kss *Tacc.tum/(B*Cavg1)
    AFIRT.Kd   = Kd  *Tacc.tum/(B*Cavg1)

    lumped_parameters_theory = data.frame(type = "theory",
                                          M30=M30,
                                          Mtot3.ss=Mtot3.ss,
                                          Tacc.tum=Tacc.tum,
                                          B = B,
                                          Cavg1 = Cavg1,
                                          Cavg3 = B*Cavg1,
                                          AFIRT.Kssd = AFIRT.Kssd,
                                          AFIRT.Kss  = AFIRT.Kss,
                                          AFIRT.Kd   = AFIRT.Kd)
    return(lumped_parameters_theory)
 }


# Function simulates the lumped parameters

lumped.parameters.simulation = function(model=model, param.as.double=param.as.double,
                                        dose.nmol=dose.nmol, tmax=tmax, tau=tau, compartment, soluble = FALSE){

    # Arguments:
    #   model_name: name of the model
    #   params_file_path: full path of the parameters file.
    #   dose.nmol: dosing amount in nmol
    #   tmax: maximum doing period in days
    #   tau: dosing interval in days
    #   compartment: compartment to which dosing is applied
    #   (in model F case, compartment=2)
    #   soluble: flag saying whether or not drug is soluble. Default to FALSE.
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
        mutate(Sfree.pct = S3/init["S3"],
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

    # Soluble case
    if (soluble) {
      AFIRT = mean(steady_state$Sfree.pct)
      Tacc.tum = mean(steady_state$Stot3)/initial_state$S3
    }

    lumped_parameters_sim = data.frame(type = "simulation",
                                     M30=M30,
                                     Mtot3.ss=Mtot3.ss,
                                     Tacc.tum=Tacc.tum,
                                     Cavg1 = Cavg1,
                                     Cavg3 = Cavg3,
                                     B     = Cavg3/Cavg1,
                                     AFIRT = AFIRT,
                                     AFIRT.sim = AFIRT) #having one named sim will be helpful later on in Task01, Task02, etc.

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
































