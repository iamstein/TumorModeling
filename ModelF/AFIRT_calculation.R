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

    # Calculate Mtot3.ss
    numerator.DM   = with(p, k13DM*(VD1/VD3)*ksynM1+(keDM1+kshedDM1+k13DM)*ksynM3)
    denomenator.DM = with(p, (keDM1+kshedDM1+k13DM)*(keDM3+kshedDM3+k31DM)-k31DM*k13DM)
    Mtot3.ss = numerator.DM / denomenator.DM

    # Calculate M30
    numerator.M    = with(p, k13M *(VD1/VD3)*ksynM1+(keM1 +kshedM1 +k13M) *ksynM3)
    denomenator.M  = with(p, (keM1 +kshedM1 +k13M) *(keM3 +kshedM3+k31M) -k31M *k13M)
    M30      = numerator.M  / denomenator.M

    # Target accumulation in the tumor compartment
    Tacc.tum = Mtot3.ss / M30

    if (!soluble){
      Kssd = with(p, (koff3 + keDM3 + kshedDM3 + k31DM)/kon3)
      Kss  = with(p, (koff3 + keDM3 + kshedDM3)        /kon3)
      Kd   = with(p,  koff3                            /kon3)
    } else {
      Kssd = with(p, (koff3 + keDS3 + kshedDM3 + k31DS)/kon3)
      Kss  = with(p, (koff3 + keDS3 + kshedDM3)        /kon3)
      Kd   = with(p,  koff3                            /kon3)

      # Mtot1.ss needed for S3tot.ss calculation.
      numerator.DS   = with(p, (kshedDM3 + k31DM + keDM3)*ksynM1 + (VD3/VD1)*k31DM*ksynM3)
      denomenator.DS = with(p, (kshedDM1 + k13DM + keDM1)*(kshedDM3 + k31DM + keDM3) - k13DM*k31DM)
      Mtot1.ss = numerator.DS / denomenator.DS

      # M10 needed for S30 calculation.
      numerator.S    = with(p, k13M *(VD1/VD3)*ksynM1+(keM1 +kshedM1 +k13M) *ksynM3)
      denomenator.S  = with(p, (keM1 +kshedM1 +k13M) *(keM3 +kshedM3+k31M) -k31M *k13M)
      M10 = numerator.S / denomenator.S

      # Calculate S3tot.ss
      numerator   = with(p, k13DS*(VD1/VD3)*(ksynS1 + kshedDM1*Mtot1.ss)+(keDS1+k13DS)*(ksynS3 + kshedDM3*Mtot3.ss))
      denomenator = with(p, (keDS1+k13DS)*(keDS3+k31DS)-k31DS*k13DS)
      Stot3.ss = numerator / denomenator

      # Calculate S30
      numerator    = with(p, k13S*(VD1/VD3)*(ksynS1 + kshedM1*M10)+(keS1 + k13S)*(ksynS3 + kshedM3*M30))
      denomenator = with(p, (keS1 + k13S)*(keS3 + k31S) - k31S*k13S)
      S30      = numerator  / denomenator

      # Target accumulation in the tumor compartment
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

    Q_2 = with(p, k12D * VD1)
    Q_3 = with(p, (k13D/k31D) * VD1)

    a_0 = with(p, (CL/VD1)*(Q_2/VD2)*(Q_3/VD3))
    a_1 = with(p, (CL/VD1)*(Q_3/VD3) + (Q_2/VD2)*(Q_3/VD3) + (Q_2/VD2)*(Q_3/VD1) + (CL/VD1)*(Q_2/VD2) + (Q_3/VD3)*(Q_2/VD1))
    a_2 = with(p, (CL/VD1)+ (Q_2/VD1) + (Q_3/VD1) + (Q_2/VD2) + (Q_3/VD3))

    P = a_1 - (a_2^2)/3
    Q = 2*((a_2^3)/27) - a_1*a_2/3 + a_0 
    r_1 = (-((P^3)/27))^0.5
    r_2 = 2*(r_1^(1/3))
  
    
    phi = acos(-Q/(2*r_1))/3
    

    alpha = -(cos(phi)*r_2 - a_2/3)
    beta = -(cos(phi + 2*pi/3)*r_2 - a_2/3)
    gamma = -(cos(phi + 4*pi/3)*r_2 - a_2/3)

    V = with(p, VD1)
    A = with(p, (1/V)*((k21D - alpha)/(alpha - beta))*((k31D - alpha)/(alpha - gamma)))
    B = with(p, (1/V)*((k21D - beta)/(beta-alpha))*((k31D - beta)/(beta - gamma)))
    C = with(p, (1/V)*((k21D - gamma)/(gamma - beta))*((k31D - gamma)/(gamma - alpha)))

    D = dose.nmol
    C_min = D*((A*exp(-alpha*tau))/(1 - exp(-alpha*tau)) + (B*exp(-beta*tau))/(1 - exp(-beta*tau)) + (C*exp(-gamma*tau))/(1 - exp(-gamma*tau)))

    if(!soluble){
        T_fold = Mtot3.ss/M30
    }else{
        T_fold = Stot3.ss/S30
    }
     
    TFIRT.Kssd = (Kssd*T_fold)/C_min
    TFIRT.Kss  = (Kss *T_fold)/C_min
    TFIRT.Kd   = (Kd  *T_fold)/C_min

    lumped_parameters_theory = data.frame(type = "theory",
                                          M30=M30,
                                          Mtot3.ss=Mtot3.ss,
                                          Tacc.tum=Tacc.tum,
                                          B = B,
                                          Cavg1 = Cavg1,
                                          Cavg3 = B*Cavg1,
                                          AFIRT.Kssd = AFIRT.Kssd,
                                          AFIRT.Kss  = AFIRT.Kss,
                                          AFIRT.Kd   = AFIRT.Kd,
                                          TFIRT.Kssd = TFIRT.Kssd,
                                          TFIRT.Kss  = TFIRT.Kss,
                                          TFIRT.Kd   = TFIRT.Kd)
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
    
    # Simulation of TFIRT
    if (soluble){
        TFIRT = max(steady_state$Sfree.pct)
    }else {
        TFIRT = max(steady_state$Mfree.pct)
    } 


    lumped_parameters_sim = data.frame(type = "simulation",
                                     M30=M30,
                                     Mtot3.ss=Mtot3.ss,
                                     Tacc.tum=Tacc.tum,
                                     Cavg1 = Cavg1,
                                     Cavg3 = Cavg3,
                                     B     = Cavg3/Cavg1,
                                     AFIRT = AFIRT,
                                     AFIRT.sim = AFIRT,
                                     TFIRT = TFIRT,
                                     TFIRT.sim = TFIRT) #having one named sim will be helpful later on in Task01, Task02, etc.

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
































