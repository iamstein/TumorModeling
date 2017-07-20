source("ams_initialize_script.R")
dirs = get.dirs(sys.calls(),dirs)

#read in the data file
  d = xlsx::read.xlsx("../../../Data/Databases/GeneralizedDatasets_Long/Bx_kon_koff_Kd.xlsx",1,stringsAsFactors=FALSE)
  d        = arrange(d,drug,param)
  d        = mutate(d,value=as.numeric(value))
  
#in cases where only a range was provided, use the mean  
  ii = which(with(d,is.na(value) & varstnam=="range"))
  d$value[ii] = (d$varstat1[ii]+d$varstat2[ii])/2
  
#for places where Kd is not available, calculate it
  ii = which(with(d,is.na(value) & param=="Kd"))
  for (i in ii) {
    r   = d[ii,]
    kon = d$value[which(d$drug   ==r$drug &
                        d$comment==r$comment &
                        d$species==r$species &
                        d$param  =="kon")]
    koff= d$value[which(d$drug   ==r$drug &
                        d$comment==r$comment &
                        d$species==r$species &
                        d$param  =="koff")]
    d$value[i] = koff/kon
  }
  
#if there is a PREFREF column equal to 1 for that particular drug, keep only those rows
  #find drugs with a prefref==1
  df = d %>%
    filter(prefref==1)
  drug1 = unique(df$drug)
  
  #where a drug has a preferred ref (prefref==1), remove all other entries with prefref==0
  d = filter(d,!(drug %in% drug1 & prefref==0))
  
#keep only the Kd data  
  d = filter(d,param=="Kd")
  
#convert all M to nM
  d$valunit = as.character(d$valunit)
  
  ii = which(d$valunit=="M")
  d$value[ii]   = d$value[ii]*1e9
  d$valunit[ii] = "nM"
  
  ii = which(d$valunit=="pM")
  d$value[ii]   = d$value[ii]/1e3
  d$valunit[ii] = "nM"
  
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
    select(drug,value) %>%
    rename(Kd.nM = value)
  
#summarize result - histogram ----
  g = ggplot(ds,aes(x=Kd.nM))
  g = g + geom_histogram(bins=10)
  g = g + scale.x.log10(1)
  gg = saveplot(8,3.3,dirs,"kon_koff_hist",draft.flag)
  grid.arrange(gg)
  
#save dataset
  write.csv(ds,paste0("../data/",dirs$outpref,"Kd_dataset.csv"))
