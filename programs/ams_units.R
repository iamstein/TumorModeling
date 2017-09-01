#unitConvert = function(d) {
  WT = c(Human=70,Cyno=3.5,Mouse=.02) #weight in kg

  #convert units
  d$value = as.numeric(d$value)
  value   = d$value
  valunit = d$valunit %>%
    str_replace("mL","ml") #replace mL with ml to simplify code

  #stop if there is any missing units
  if(any(is.na(d$valunit))) {
    d %>% filter(is.na(valunit)) %>% print()
    stop("some units are missing, need to be filled in")
  }
  
  #find which units are not standard
  standard.unit = c("-","1/d","1/nM/d","kDa","L","L/d","nM","nM/d")
  ii = which(!(valunit %in% standard.unit))
  print(sort(unique(valunit[ii])))
  
  for (i in ii) {
    valuniti = valunit[i] 
    valuei   = value[i]
    parami   = d$param[i]
    drugi    = d$drug[i]
    refi     = d$ref[i]
    Vci      = NULL
    MWi      = NULL
    print(sprintf("%d,%s,%s,%1.2f,%s",i,drugi,parami,valuei,valuniti))
    if (!is.na(parami)) {
    
      #for any molar conversion, we will need the molecular weight
      if (str_detect(valuniti,"^[munp]g")) {
        di  = filter(d,drug==drugi)
        if (parami %in% c("Vm","Km","Cminss","IC50","EC50")) {
          MWi = with(di,value[param=="MWD" & !is.na(param)])
        } else if (parami %in% c("KdS","ksynS","S0")) {
          MWi= with(di,value[param=="MWS" & !is.na(param)]) 
        } else if (parami %in% c("KdM","ksynM","M0")) {
          MWi= with(di,value[param=="MWM" & !is.na(param)])
        } else {
          stop("not sure which molecular weight to use")
        }
        stopifnot(!(is.null(MWi)==1 | length(MWi)>1 | length(MWi)==0))
      }
      
      #for per kg units or Vm, we'll need the central volume
      if (str_detect(valuniti,"/kg") | parami=="Vm") {
        di  = filter(d,ref==refi,drug==drugi)
        WTi = WT[di$species[1]]
        Vci = with(di,value[param=="Vc"])
        Vci.unit = with(di,valunit[param=="Vc"])
        if (str_detect(Vci.unit,"ml")) { #make sure Vci is in L;
          Vci = Vci/1000  }
        if (str_detect(Vci.unit,"/kg")) { #make sure Vci is not per kg
          Vci = Vci*WTi }
        stopifnot(length(Vci)==1)
      }
      
    #percentage    
      if (valuniti=="%") {
        valuei   = valuei/100
        valuniti = "-"
      }
    
    #time      
      if (str_detect(valuniti,"/h")) {
        valuei   = valuei*24
        valuniti = str_replace(valuniti,"/h","/d")
      } else if (str_detect(valuniti,"/min")) {
        valuei   = valuei*24*60
        valuniti = str_replace(valuniti,"/min","/d")      
      } else if (str_detect(valuniti,"/s")) {
        valuei   = valuei*24*60*60
        valuniti = str_replace(valuniti,"/s","/d")      
      }
      
    #per kg
      if (str_detect(valuniti,"/kg")) {
        valuei   = valuei*WTi
        valuniti = str_replace(valuniti,"/kg","")
      }      
      
    #volume
      if (str_detect(valuniti,"/ml")) {
        valuei = valuei*1000
        valuniti = str_replace(valuniti,"/ml","/L")
      } else if (str_detect(valuniti,"ml")) {
        valuei = valuei/1000
        valuniti = str_replace(valuniti,"ml","L")
      }
      
    #concentration
      if (str_detect(valuniti,"nM")) {
        #do nothing
      } else if (str_detect(valuniti,"/pM")) {
        valuei = valuei*1e3
        valuniti = str_replace(valuniti,"/pM","/nM")
      } else if (str_detect(valuniti,"pM")) {
        valuei = valuei/1e3
        valuniti = str_replace(valuniti,"pM","nM")
      } else if (str_detect(valuniti,"/M")) {
        valuei = valuei/1e9
        valuniti = str_replace(valuniti,"/M","/nM")
      } else if (str_detect(valuniti,"M")) {
        valuei = valuei*1e9
        valuniti = str_replace(valuniti,"M","nM")
      } 
      
    #mass
      if (str_detect(valuniti,"mg")) {
        valuei = valuei*1000
        valuniti = str_replace(valuniti,"mg","ug")
      } else if (str_detect(valuniti,"ng")) {
        valuei = valuei/1000
        valuniti = str_replace(valuniti,"ng","ug")
      } else if (str_detect(valuniti,"pg")) {
        valuei = valuei/1e6
        valuniti = str_replace(valuniti,"pg","ug")
      }
  
    #convert ug/L or ug to nM
      if (str_detect(valuniti,"ug/L")) {
        valuei = valuei/MWi
        valuniti = str_replace(valuniti,"ug/L","nM")
      } else if (str_detect(valuniti,"ug")) {
        valuei = valuei/MWi/Vci
        valuniti = str_replace(valuniti,"ug","nM")
      }
      
    #final checks
      stopifnot((length(valuniti)==1))
      stopifnot(valuniti %in% standard.unit)    
      
    #put into vector  
      value[i]   = valuei
      valunit[i] = valuniti
    }
  }
  d$valunit = valunit
  d$value   = value

  print(sort(unique(d$valunit)))

#  return(d)
#}
  
  
  
  
  