#initialize
  source("ams_initialize_script.R")
  d = read.csv("../data/Chen16_CNTO345_GenDataset.csv",stringsAsFactors = FALSE)
  
  d = d %>% 
    filter(NAME!="dose") %>%
    mutate(TYPENAME = factor(TYPENAME,levels=c("Serum","Joint Lavage")),
           NAME     = mapvalues(NAME,c("pk","total target","free target"),c("Drug","Total Target","Free Target")),
           NAME     = factor(NAME    ,levels=c("Drug","Total Target","Free Target")),
           DOSE     = factor(DOSE    ,levels=c(30,3,.3)))
  
  
  
# plot data ----  
  g = ggplot(d,aes(x=STARTIME,y=VALUE,color=DOSE))
  g = g + geom_point(size=.5)
  g = g + stat_summary(fun.y=median,geom="line")
  #g = g + geom_smooth()
  g = g + scale.y.log10()
  g = g + facet_grid(TYPENAME~NAME,switch="y")
  g = g + geom_hline(aes(yintercept=LLIM),linetype="dotted")
  g = xscale("d21",increment=7)
  g = g + labs(y="Conc. (nM)",color="Dose\nmg/kg")
  print(g)
  
  ggsave("../results/Task02_CNTO_IL6.pdf",width=9,height=4)
  