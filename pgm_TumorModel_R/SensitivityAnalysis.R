# The function basic.sensitivity.analysis returns a data frame of three columns:
# column 1 = time, 
# column 2 = variable
# column 3 = quantity that the sensitivity analysis is performed on
# To see an example of the use of the function basic.sensitivity, 
# you can run the script "BasicSensitivityAnalysisExample.R" located in the same directory 

lseq = function(from, to, length.out){
    sequence = seq(log(from), log(to), length.out=length.out)
    sequence = exp(sequence)
    return(sequence)
}

basic.sensitivity.analysis = function(parameter.file, # string, full path of the parameter file
                                      model.file, # string, full path of the model file
                                      model, #function, function name of the model
                                      quantity, # string, the quantity you want to perform sensitity analysis on
                                      variable, # string, the variable you want to perform sensitity analysis against
                                      variable.from, # double, initial variable value
                                      variable.to, # double, final variable value
                                      fold.number, # integer, number of folds for the sensitity analysis
                                      dose.nmol, # double, dose amount in nmol
                                      tau, # integer, dosing interval
                                      tmax, # integer, total number of days the observation is conducted
                                      compartment # integer, the number of the compartment to which dosing applied
                                      ){
  source("ams_initialize_script.R")
  source(model.file)

  d <- read.csv(parameter.file)
  param.as.double <- d$Value
  names(param.as.double) <- d$Parameter
  
  model = model
  
  ev = eventTable(amount.units="nmol", time.units = "days")
  sample.points = c(seq(-7, tmax, 0.1), 10^(-3:0)) # sample time, increment by 0.1
  sample.points = sort(sample.points)
  sample.points = unique(sample.points)
  ev$add.sampling(sample.points)
  ev$add.dosing(dose=dose.nmol, nbr.doses=floor(tmax/tau)+1, dosing.interval=tau,
                dosing.to=compartment)
  
  variable.range = seq(log(variable.from), log(variable.to), length.out=fold.number)
  variable.range = exp(variable.range)
  OUT = data.frame()
  for (value in variable.range){
    param.as.double[variable]=value
    init = model$init(param.as.double)
    out = model$rxode$solve(param.as.double, ev, init)
    out = model$rxout(out)
    out = out %>% 
      mutate(variable = value)
    OUT = rbind(OUT, out)
  }
  
  OUT = OUT[,c("time", "variable", quantity)]
  colnames(OUT)=c("time", variable, quantity)
  return(OUT) 
}
    
 
plot.basic.sensitivity.anlysis = function(data){
    names = names(data)
    g = ggplot(data, ase(x=time, y=data[,2], color=data[,1])) +
        scale.x.log10() +
        sclae.y.log10() + 
        geom_point() + 
        ylab(names[3]) + 
        xlab(names[1])
    return(g)
}

                        
  




