# Helper function that returns a range of variable when performing
# sensitvity analysis

read.param.file = function(filename) {
  d                      = read_excel(filename, 1)
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

lumped.parameters.theory = function(param.as.double = param.as.double,
                                    dose.nmol       = dose.nmol,
                                    tau             = tau,
                                    soluble         = FALSE){
    # Arguments:
    #   params_file_path: full path of the parameters file.
    #   dose.nmol: dosing amout in nmol
    #   tau: dosing interval in days
    #   soluble: flag saying whether or not drug is soluble. Default to FALSE.
    # Return:
    #   A data frame of lumped parameters calculated from theory

    pars    = as.data.frame(t(param.as.double))

    # Calculate Mtot3.ss
    numerator.DM   = with(pars, k13DM*(VD1/VD3)*ksynM1+(keDM1+kshedDM1+k13DM)*ksynM3)
    denomenator.DM = with(pars, (keDM1+kshedDM1+k13DM)*(keDM3+kshedDM3+k31DM)-k31DM*k13DM)
    Mtot3.ss = numerator.DM / denomenator.DM

    # Calculate M30
    numerator.M    = with(pars, k13M *(VD1/VD3)*ksynM1+(keM1 +kshedM1 +k13M) *ksynM3)
    denomenator.M  = with(pars, (keM1 +kshedM1 +k13M) *(keM3 +kshedM3+k31M) -k31M *k13M)
    M30      = numerator.M  / denomenator.M

    # Target accumulation in the tumor compartment
    Tacc.tum = Mtot3.ss / M30

    if (!soluble){
      Kssd = with(pars, (koff3 + keDM3 + kshedDM3 + k31DM)/kon3)
      Kss  = with(pars, (koff3 + keDM3 + kshedDM3)        /kon3)
      Kd   = with(pars,  koff3                            /kon3)
    } else {
      Kssd = with(pars, (koff3 + keDS3 + kshedDM3 + k31DS)/kon3)
      Kss  = with(pars, (koff3 + keDS3 + kshedDM3)        /kon3)
      Kd   = with(pars,  koff3                            /kon3)

      # Mtot1.ss needed for S3tot.ss calculation.
      numerator.DS   = with(pars, (kshedDM3 + k31DM + keDM3)*ksynM1 + (VD3/VD1)*k31DM*ksynM3)
      denomenator.DS = with(pars, (kshedDM1 + k13DM + keDM1)*(kshedDM3 + k31DM + keDM3) - k13DM*k31DM)
      Mtot1.ss = numerator.DS / denomenator.DS

      # M10 needed for S30 calculation.
      numerator.S    = with(pars, k13M *(VD1/VD3)*ksynM1+(keM1 +kshedM1 +k13M) *ksynM3)
      denomenator.S  = with(pars, (keM1 +kshedM1 +k13M) *(keM3 +kshedM3+k31M) -k31M *k13M)
      M10 = numerator.S / denomenator.S

      # Calculate S3tot.ss
      numerator   = with(pars, k13DS*(VD1/VD3)*(ksynS1 + kshedDM1*Mtot1.ss)+(keDS1+k13DS)*(ksynS3 + kshedDM3*Mtot3.ss))
      denomenator = with(pars, (keDS1+k13DS)*(keDS3+k31DS)-k31DS*k13DS)
      Stot3.ss = numerator / denomenator

      # Calculate S30
      numerator    = with(pars, k13S*(VD1/VD3)*(ksynS1 + kshedM1*M10)+(keS1 + k13S)*(ksynS3 + kshedM3*M30))
      denomenator = with(pars, (keS1 + k13S)*(keS3 + k31S) - k31S*k13S)
      S30      = numerator  / denomenator

      # Target accumulation in the tumor compartment
      Tacc.tum = Stot3.ss / S30
    }

    # Biodistribution coefficient (reference: ModelF_Appendix)
    B = with(pars, (k13D/(keD3 + k31D) * (VD1/VD3)))

    # Clearance
    CL = with(pars, (keD1*VD1))

    # Average drug concentration in the central compartment
    Cavg1 = dose.nmol/(CL*tau)

    # Compute various AFIRTs
    AFIRT.Kssd = Kssd*Tacc.tum/(B*Cavg1)
    AFIRT.Kss  = Kss *Tacc.tum/(B*Cavg1)
    AFIRT.Kd   = Kd  *Tacc.tum/(B*Cavg1)

    Q2 = with(pars, k12D * VD1)
    Q3 = with(pars, (k13D/k31D) * VD1)

    a0 = with(pars, (CL/VD1)*(Q2/VD2)*(Q3/VD3))
    a1 = with(pars, (CL/VD1)*(Q3/VD3) + (Q2/VD2)*(Q3/VD3) + (Q2/VD2)*(Q3/VD1) + (CL/VD1)*(Q2/VD2) + (Q3/VD3)*(Q2/VD1))
    a2 = with(pars, (CL/VD1)+ (Q2/VD1) + (Q3/VD1) + (Q2/VD2) + (Q3/VD3))

    p  = a1 - (a2^2)/3
    q  = 2*(a2^3)/27 - a1*a2/3 + a0 
    r1 = (-(p^3)/27)^0.5
    r2 = 2*(r1^(1/3))
    
    phi = acos(-q/(2*r1))/3
    

    alpha = -(cos(phi)         *r2 - a2/3)
    beta  = -(cos(phi + 2*pi/3)*r2 - a2/3)
    gamma = -(cos(phi + 4*pi/3)*r2 - a2/3)

    V = with(pars, VD1)
    A = with(pars, (1/V) * ((k21D - alpha)/(alpha - beta)) * ((k31D - alpha)/(alpha - gamma)))
    B = with(pars, (1/V) * ((k21D - beta )/(beta - alpha)) * ((k31D - beta )/(beta - gamma )))
    C = with(pars, (1/V) * ((k21D - gamma)/(gamma - beta)) * ((k31D - gamma)/(gamma - alpha)))

    D = dose.nmol
    Cmin = D*((A*exp(-alpha*tau))/(1 - exp(-alpha*tau)) + 
              (B*exp(-beta *tau))/(1 - exp(-beta *tau)) + 
              (C*exp(-gamma*tau))/(1 - exp(-gamma*tau)))

    if(!soluble){
        Tfold = Mtot3.ss/M30
    }else{
        Tfold = Stot3.ss/S30
    }
     
    TFIRT.Kssd = Kssd*Tfold/(B*Cmin)
    TFIRT.Kss  = Kss *Tfold/(B*Cmin)
    TFIRT.Kd   = Kd  *Tfold/(B*Cmin)

    lumped_parameters_theory = data.frame(type       = "theory",
                                          M30        = M30,
                                          Mtot3.ss   = Mtot3.ss,
                                          Tacc.tum   = Tacc.tum,
                                          B          = B,
                                          Cavg1      = Cavg1,
                                          Cavg3      = B*Cavg1,
                                          AFIRT.Kssd = AFIRT.Kssd,
                                          AFIRT.Kss  = AFIRT.Kss,
                                          AFIRT.Kd   = AFIRT.Kd,
                                          TFIRT.Kssd = TFIRT.Kssd,
                                          TFIRT.Kss  = TFIRT.Kss,
                                          TFIRT.Kd   = TFIRT.Kd)
    return(lumped_parameters_theory)
 }


# Function simulates the lumped parameters

lumped.parameters.simulation = function(model           = model, 
                                        param.as.double = param.as.double,
                                        dose.nmol       = dose.nmol, 
                                        tmax            = tmax, 
                                        tau             = tau, 
                                        compartment,
                                        soluble         = FALSE){

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

    ## Average drug concentration in central compartment
    dose_applied = out %>%
        filter(time > 0)
    Cavg1 = mean(dose_applied$D1)

    # Average drug concentration in tumor compartment
    Cavg3 = mean(dose_applied$D3)

    # AFIRT and target accumulation
    if (soluble) {
      AFIRT = mean(steady_state$Sfree.pct)
      Tacc.tum = mean(steady_state$Stot3)/initial_state$S3
    } else {
      AFIRT = mean(steady_state$Mfree.pct)
      Tacc.tum = Mtot3.ss / M30
    }
    
    # Simulation of TFIRT
    if (soluble) {
        TFIRT = max(steady_state$Sfree.pct)
    } else {
        TFIRT = max(steady_state$Mfree.pct)
    } 


    lumped_parameters_sim = data.frame(type      = "simulation",
                                       M30       = M30,
                                       Mtot3.ss  = Mtot3.ss,
                                       Tacc.tum  = Tacc.tum,
                                       Cavg1     = Cavg1,
                                       Cavg3     = Cavg3,
                                       B         = Cavg3/Cavg1,
                                       AFIRT     = AFIRT,
                                       AFIRT.sim = AFIRT,
                                       TFIRT     = TFIRT,
                                       TFIRT.sim = TFIRT) #having one named sim will be helpful later on in Task01, Task02, etc.

    return(lumped_parameters_sim)
}


simulation = function(model           = model, 
                      param.as.double = param.as.double,
                      dose.nmol       = dose.nmol, 
                      tmax            = tmax, 
                      tau             = tau){
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

# This function does the sensitivity analysis on the user inputted parameter and compares the theoretical result to the simulated result.

# Input:
# model - model system of ODE's solved with RxODE. In this project, it is 'ivsc_4cmtct_shedct'.
# param.as.double - read parameters from Excel file. read.param.file("file directory").
# dose.nmol - dose in nmol.
# tmax - time of treatment in days
# tau - frequency of administering dose in days
# compartment - compartment where drug is administered
# param.to.change - parameter on which to do SA. This must be a string.
# param.to.change.range - range of parameter on which to do SA. The range must be symmetric in fold change. This must be a vector of odd length.
# soluble - boolean that is true/false if the drug is soluble/insoluble. Need this since soluble and insoluble are treated differently.

# Output:
# Data frame of AFIRT vs parameter value

compare.thy.sim = function(model = model,
                           param.as.double = param.as.double,
                           dose.nmol = dose.nmol,
                           tmax = tmax,
                           tau = tau,
                           compartment = compartment,
                           param.to.change = param.to.change,
                           param.to.change.range = param.to.change.range,
                           soluble = FALSE) {

  # --------------------------------------------------------------------------------
  # Simulation
  # --------------------------------------------------------------------------------

  df_sim = data.frame()
  # Iterate through values in range.
  if (param.to.change == 'dose'){
    for (param.iter in param.to.change.range){
      row = lumped.parameters.simulation(model, param.as.double, param.iter, tmax, tau, compartment, soluble)
      df_sim = rbind(df_sim, row)
    }
  } else {
    for (param.iter in param.to.change.range){
      param.as.double[param.to.change] = param.iter
      row = lumped.parameters.simulation(model, param.as.double, dose.nmol, tmax, tau, compartment, soluble)
      df_sim = rbind(df_sim, row)
    }
  }
  df_sim = df_sim %>% mutate(param.to.change = param.to.change.range,
                             fold.change = param.to.change.range/median(param.to.change.range))

  # --------------------------------------------------------------------------------
  # Theory
  # --------------------------------------------------------------------------------

  df_thy = data.frame()

  # Iterate through values in range.
  if (param.to.change == 'dose'){
    for (param.iter in param.to.change.range){
      row = lumped.parameters.theory(param.as.double, param.iter, tau, soluble)
      df_thy = rbind(df_thy, row)
    }
  } else{
    for (param.iter in param.to.change.range){
     param.as.double[param.to.change] = param.iter
     row = lumped.parameters.theory(param.as.double, dose.nmol, tau, soluble)
    df_thy = rbind(df_thy, row)
    }
  }

  df_thy = df_thy %>% mutate(param.to.change = param.to.change.range,
                             fold.change = param.to.change.range/median(param.to.change.range))

  # --------------------------------------------------------------------------------
  # Arrange theory and simulation in single data frame.
  # --------------------------------------------------------------------------------

  # I am tired of that "Unequal factor levels" error. This fixes it.
  levels(df_thy$type) = c("theory", "simulation")
  levels(df_sim$type) = c("theory", "simulation")
  df_compare = bind_rows(df_thy,df_sim)
  param = param.to.change
  df_compare = df_compare %>%
    mutate(param = param) %>%
    arrange(param.to.change,type) %>%
    mutate_if(is.numeric,signif,2)
  return(df_compare)
}