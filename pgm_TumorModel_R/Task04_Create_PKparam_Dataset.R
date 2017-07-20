source("ams_initialize_script.R")
dirs = get.dirs(sys.calls(),dirs)

#read in the data file
  d = xlsx::read.xlsx("../../../Data/Databases/GeneralizedDatasets_Long/Bx_PK_Binding.xlsx",1,stringsAsFactors=FALSE)
  d        = arrange(d,drug,param)
  d        = mutate(d,value=as.numeric(value))
  
#if there is a PREFREF column equal to 1 for that particular drug, keep only those rows
  #find drugs with a prefref==1
  df = d %>%
    filter(prefref==1)
  drug1 = unique(df$drug)
  
  #where a drug has a preferred ref (prefref==1), remove all other entries with prefref==0
  d = filter(d,!(drug %in% drug1 & prefref==0))
  
#convert all /h units to /d units
  d$valunit = as.character(d$valunit)
  
  ii = which(str_detect(d$valunit,"/h$"))
  d$value[ii]   = d$value[ii]*24
  d$valunit[ii] = str_replace_all(d$valunit[ii],"/h","/d")
  
#convert all pM units to nM
  d$valunit = as.character(d$valunit)
  
  ii = which(str_detect(d$valunit,"pM"))
  d$value[ii]   = d$value[ii]/1000
  d$valunit[ii] = str_replace_all(d$valunit[ii],"pM","nM")  
  
#convert all mL units to L units
  ii = which(str_detect(d$valunit,"^mL") | str_detect(d$valunit,"^ml"))
  d$value[ii]   = d$value[ii]/1000
  d$valunit[ii] = str_replace_all(d$valunit[ii],"mL","L")  
  d$valunit[ii] = str_replace_all(d$valunit[ii],"ml","L")  
  
#if units are /kg, multiply by 70
  ii = which(str_detect(d$valunit,"/kg"))
  d$value[ii]   = d$value[ii]*70
  d$valunit[ii] = str_replace_all(d$valunit[ii],"/kg","")  
  
  print("UNITS")
  print(sort(unique(d$valunit)))  
  
#remove any NA rows
  d = filter(d,!is.na(drug))
  dall = d
  
#average data together when there are multiple entries
  d = dall %>%
    select(ref,drug,param,value,valunit) %>%
    group_by(drug,param,valunit) %>%
    summarize(value = mean(value)) %>%
    ungroup()
    
#spread data so one row per drug
  ds = d %>%
    ungroup() %>%
    select(drug,param,value) %>%
    filter(param %in% c("ka","F","CL","Q","Vc","Vp","T0")) %>%
    spread(param,value) %>%
    rename(CL.L_d=CL,ka._d=ka,Q.L_d=Q,Vc.L=Vc,Vp.L=Vp,T0.nM=T0)
  
#summarize result - histogram ----
  g = ggplot(ds,aes(x=CL.L_d))
  g = g + geom_histogram(bins=10)
  g = g + scale.x.log10(.5)
  gg = saveplot(8,3.3,dirs,"CL_hist",draft.flag)
  grid.arrange(gg)
  
#summarize result - histogram ----
  g = g + aes(x=T0.nM)
  g = g + scale.x.log10(1)
  gg = saveplot(8,3.3,dirs,"T0_hist",draft.flag)
  grid.arrange(gg)
  
#save dataset
  write.csv(ds,paste0("../data/",dirs$outpref,"PKparam_dataset.csv"))
