source("ams_initialize_script.R")

#read in the data files
  source("Task03_Create_Kd_Dataset.R")
  source("Task04_Create_PKparam_Dataset.R")

  dirs = get.dirs(sys.calls(),dirs)


  d.Kd = read.csv("../data/Task03_Kd_dataset.csv",stringsAsFactors = FALSE)
  d.PK = read.csv("../data/Task04_PKparam_dataset.csv",stringsAsFactors = FALSE)
  d.info = xlsx::read.xlsx("../../../Data/Databases/GeneralizedDatasets_Long/mAb_Approved.xlsx",1,stringsAsFactors=FALSE)
  d.info = d.info %>%
    select(drug,target,target_type=type)
  
#read in the dosing file (for doses and dosing intervals)  
  d = xlsx::read.xlsx("../../../Data/Databases/Tables/mAb_Dosing_Indication_Table.xlsx",1,stringsAsFactors=FALSE)
  names(d) = tolower(names(d))
  d.dose  = d %>%
    arrange(drug,indicat) %>%
    filter(age=="Adult") %>%
    mutate(dose.mg = ifelse(doseunit=="mg/kg",
                            dose*70,
                            dose),
           dose.mg = ifelse(doseunit=="mg/m^2",
                            dose*1.75,
                            dose.mg),
           dose.nmol = dose.mg*scale.mg2nmol) %>%
    select(drug,route,indicat,indshort,dose.mg,dose.nmol,tau.d=doseint_day)

#merge dosing regimen, PK, and Kd params    
  d = d.dose %>%
    left_join(d.info,by="drug") %>%
    left_join(d.PK,by="drug") %>%
    left_join(d.Kd,by="drug")
  
#compute AFIRT ----
  d = d %>%
    mutate(dose.mg_weekly = dose.mg/tau.d*7,
           Cavg.nM = ifelse(route=="iv",
                            dose.nmol/(CL.L_d*tau.d),
                          F*dose.nmol/(CL.L_d*tau.d)),
           AFIRT   = 3*Kd.nM/Cavg.nM,
           Target_Drug = paste0(target,":",drug),
           TargetEnterDrug = paste0(target,"\n",drug)) %>%
    rename(Indication=indicat)
  
#plot AFIRT ----
  dplot = d %>% 
    filter(target_type=="soluble") %>%
    mutate(T0_flag = T0.nM*10<Cavg.nM)
  
  g = ggplot(dplot,aes(x=Target_Drug,y=AFIRT,shape=T0_flag))
  g = g + geom_jitter(width=.3,height=0,alpha=.5,size=3)
  g = g + xlab("Target and Drug")
  g = g + ylab("AFIRT: Predicted\nFree Fraction In Tissue ")
  g = g + labs(shape="Cavg>>T0")
  g = g + theme(axis.text.x = element_text(angle = 60, hjust = 1))
  g = g + scale.y.log10(.5,limits=c(min(c(1e-4,d$AFIRT)),max(1e-2,d$AFIRT)))
  g = g + scale_color_discrete(guide=FALSE)
  g = g + scale_shape_manual(values=c(`FALSE`=17,`TRUE`=16))
  g = g + ggtitle("AFIRT = Kd/(.33*Cavg)\n= Kd/(.33*F*Dose/(CL*tau))")
  gg = saveplot(10,6,dirs,"AFIRT",draft.flag)
  grid.arrange(gg)
  
#plot AFIRT with Disease Labels ----
  dplot = d %>% 
    filter(target_type=="soluble") %>%
    mutate(T0_flag = T0.nM*10<Cavg.nM,
           hline   = .1) %>%
    filter(T0_flag=="TRUE") %>%
    arrange(Target_Drug)
  
  g = ggplot(dplot,aes(x=Target_Drug,y=AFIRT,label=indshort))
  g = g + annotate("text",x=dplot$Target_Drug[3],y=.15,label=paste0((1-dplot$hline[1])*100,"% inhib"),color="blue")
  g = g + geom_hline(aes(yintercept=hline),color="blue",linetype="dashed")
  
  g = g + geom_point(alpha=.5,size=3)
  g = g + geom_text_repel(size=2.5,show.legend = FALSE)#width=.3,height=0,alpha=.5,size=3)
  g = g + xlab("Target and Drug")
  g = g + ylab("AFIRT (% Inhib)")
  #g = g + labs(color="Cavg>>T0")
  g = g + theme(axis.text.x = element_text(angle = 60, hjust = 1))
  
  breaks = 10^(-5:0)
  labels = paste0(breaks," (",100-breaks*100,"%)")
  g = g + scale_y_log10(limits=c(min(c(1e-4,d$AFIRT)),1),breaks=breaks,labels=labels)#max(1e-2,d$AFIRT)))
  #g = g + scale_color_discrete(guide=FALSE)
  #g = g + scale_color_manual(values=c(`FALSE`="red",`TRUE`="black"))
  g = g + ggtitle("Average Free target to \nInitial target Ratio in Tissue (AFIRT)")
  
  gg = saveplot(6,5,dirs,"AFIRT_Indicat",draft.flag)
  grid.arrange(gg)  

#print the indshort labels ----
  d1 = dplot %>%
    subset(!duplicated(indshort)) %>%
    arrange(indshort) %>%
    select(indshort,Indication)
  
  gg = tableGrob(d1,rows=NULL,cols = NULL,
                 theme = ttheme_default(core=list(fg_params=list(hjust=0,x=0)),
                                        colhead=list(fg_params=list(hjust=0,x=0))))
  grid.arrange(gg)
  ggsave("../results/Task05_IndicationKey.png",plot=gg,width=7,height=7)
  
#plot weekly dose ----
  g = g + aes(y=dose.mg_weekly)
  g = g + scale_y_continuous()
  g = g + scale.y.log10()
  #g = g + ylim(c(0,max(c(80,d$dose.mg_weekly))))
  g = g + scale.y.log10(.5)
  g = g + ylab("Dose per Week (mg/week)\nFor 70 kg patient")
  gg = saveplot(7,7,dirs,"Dose",draft.flag)
  grid.arrange(gg)
  
#plot average drug concentration vs T0   ----
  d1 = subset(d,!duplicated(d$drug))
  
  lim = c(min(c(d1$T0.nM,d1$Cavg.nM),na.rm=TRUE),
          max(c(d1$T0.nM,d1$Cavg.nM),na.rm=TRUE))
  
  g = ggplot(d1,aes(x=T0.nM,y=Cavg.nM,label=TargetEnterDrug))
  g = g + annotate("segment",x=lim[1],xend=lim[2],y=lim[1],yend=lim[2])
  g = g + geom_text_repel(size=2.5)
  g = g + geom_point(size=3)
  
  g = g + scale.y.log10(1,limits=lim)
  g = g + scale.x.log10(1,limits=lim)
  g = g + no.legend()
  g = g + ggtitle("Baseline Conc. vs Average Drug Conc.")
  g = g + labs(x="Baseline Target Conc. (nM)",y="Average Drug Conc.")
  #g = g + theme(panel.grid=element_line(color="black"),axis.line=element_line(color="green"))
  gg = saveplot(5,5,dirs,"T0_vs_Cavg",draft.flag)
  grid.arrange(gg)
  
#save final dataset for browsing ----
  dsave = d %>%
    select(drug,target,target_type,route,indication=Indication,indshort,dose.mg,tau.d,dose.mg_weekly,F,ka._d,CL.L_d,Q.L_d,Vc.L,Vp.L,Kd.nM,T0.nM,Cavg.nM,AFIRT)
  
  write.csv(dsave,"../data/Task05_Param_Summary_AFIRT.csv")
  
#NEXT STEPS - ADD XOLAIR, AND TALK TO PHIL
