source("ams_initialize_script.R")
source("ivsc_3cmtct_shed3.R")
dirs = get.dirs(sys.calls(),dirs)

mod      = ivsc_3cmtct_shed3()
d        = xlsx::read.xlsx("ivsc_3cmtct_shed3_param.xlsx",1,stringsAsFactors=FALSE)
d        = d %>% filter(!is.na(Parameter))
param    = d
row.names(param) = param$Parameter
p        = d$Value
names(p) = d$Parameter

#these functions aren't strictly needed here.  just testing because will be useful later
  p        = p[mod$pin]
  p        = mod$repar(p)
  
  print(t(data.frame(t(p))))

#simulate model and put into OUT
  
n.dose = 6 
dose.range = c(8,80,800) #mg
dose.tau= 2 #dose every tau weeks
OUT = list()
i   = 0
for (dose.mg in dose.range) {  
  ev       = eventTable(amount.units="nmol", time.units = "days")
  ev$add.sampling(unique(sort(c(seq(-7,13*7,.1),10^(-3:0)))))
  ev$add.dosing(dose=dose.mg*scale.mg2nmol,nbr.doses=n.dose,dosing.interval=dose.tau*7,dosing.to=2)
  
  init        = mod$init(p)
  out         = mod$rxode$solve(p, ev, mod$init(p))
  out         = mod$rxout(out,p)
  out = out %>%
    mutate(Sfree.pct = S1/init["S1"],
           Mfree.pct = M3/init["M3"],
           dose.mg   = dose.mg,
           dose.nM   = dose.mg/1000/150e3*1e9) #mg-->g-->mol-->nMol
  i=i+1
  OUT[[i]] = out
}
out = bind_rows(OUT)

#plot data ----
d = out %>%
  select(time,D1,D3,S1,S3,Stot1,Stot3,M3,Mtot3,Sfree.pct,Mfree.pct,dose.mg,dose.nM) %>%
  gather(variable,value,-c(time,dose.mg,dose.nM)) %>%
  mutate(value = ifelse(value>1e-8,value,NA))
  
var2nice = function(d) {
  d = d %>%  mutate(
         varnice = variable,
         varnice = str_replace_all(varnice,"1",".pla"),
         varnice = str_replace_all(varnice,"3",".tum"),
         varnice = str_replace_all(varnice,"D","Drug."),
         varnice = str_replace_all(varnice,"S","s-Target"),
         varnice = str_replace_all(varnice,"M","m-Target"),
         varnice = str_replace_all(varnice,".pct","%"),
         molecule= paste0(str_extract(varnice,"^\\w+."),"[nM]")
      )
  d$molecule[str_detect(d$varnice,"free")]="free.ratio"
  d$molecule = factor(d$molecule,levels=unique(d$molecule))
  return(d)
}
  
popt      = data.frame(variable = c("D1"   ,"D3"    ,"S1"   ,"S3"    ,"M3"    ,"Stot1","Stot3"  ,"Mtot3" ,"Sfree.pct" ,"Mfree.pct"),
                       color    = c("black","gray50","red"  ,"orange","blue"  ,"red4" ,"orange4","blue4" ,"pink"       ,"lightblue"),
                       linetype = c("solid","solid" ,"solid","solid" ,"solid" ,"solid","solid"  ,"solid" ,"dashed"     ,"dashed"),
                       stringsAsFactors = FALSE)
d    = var2nice(d)
popt = var2nice(popt)

color     = popt$color
names(color) = popt$varnice

linetype  = popt$linetype
names(linetype) = popt$varnice

#create additional plot I need for Mike Roy ---- 
dd = d %>%
  filter(variable=="D1") %>%
  mutate(dose.mgstr = paste0(dose.mg," mg q",dose.tau,"w"))
dd$dose.mgstr = factor(dd$dose.mgstr,levels=sort(unique(dd$dose.mgstr),decreasing = TRUE))
g = ggplot(dd,aes(x=time,y=value,group=dose.mgstr,color=dose.mgstr))
g = g + geom_line(size = 1)
g = xscale("d100",increment = 14)
g = g + scale.y.log10(1)
g = g + ylab("Conc [nM]")
gg = saveplot(4,4,dirs,"PK_Sim",draft.flag)
grid.arrange(gg)

#plot of PK-Target Engagement ----
pp  = as.data.frame(t(p))
CL  = pp$keD1*pp$VD1 #L/d
tau = 14 #days 
ABCisf = param["ABCisfD","Value"]
Kd  = pp$koff1/pp$kon1
Tacc.tum = 1

dose = d %>%
  subset(!duplicated(dose.mg)) %>%
  select(dose.nM,dose.mg)

Cavg = d %>%
  filter(varnice == "Drug..pla",time>28*2 & time<=28*3) %>%
  group_by(dose.nM) %>%
  summarise(Cavg.pla.ss = mean(value)) %>%
  ungroup() %>%
  mutate(Cavg.pla.ss.thry = dose.nM/CL/tau,
         Cavg.tum.ss.thry = ABCisf*Cavg.pla.ss.thry)

Cavgtum = d %>%
  filter(varnice == "Drug..tum",time>28*2 & time<=28*3) %>%
  group_by(dose.nM) %>%
  summarise(Cavg.tum.ss = mean(value))

AFIRT = d %>%
  filter(varnice == "m-Targetfree%",time>28*2 & time<=28*3) %>%
  group_by(dose.nM) %>%
  summarise(AFIRT.m = mean(value)) %>%
  ungroup()

dsumm = Cavg %>%
  left_join(Cavgtum,by="dose.nM") %>%
  left_join(dose,by="dose.nM") %>%
  left_join(AFIRT,by="dose.nM") %>%
  mutate(AFIRT.m.thry = Kd*Tacc.tum/(Cavg.pla.ss*ABCisf),
         molecule = factor("free.ratio")) 

#print(signif(dsumm,2))
  

#plot results ----   
g = ggplot(d,aes(x=time,y=value,group=varnice,color=varnice))
  g = g + facet_grid(molecule~dose.mg,scales = "free_y",switch = "y")
  g = g + geom_line(size=1)
  g = g + geom_hline(data=dsumm,aes(yintercept=AFIRT.m.thry),color="dodgerblue")
  #g = g + geom_line(aes(linetype=variable))
  #g = g + geom_point(aes(y=value),size=3)
  g = g + scale.y.log10()
  #g = g + scale.y.log10(limits=c())
  g = g + scale_color_manual(   values=color)
  #g = g + scale_linetype_manual(values=linetype)
  g = xscale("m3")
  g = g + ggtitle(paste0("Dosing (mg): q", dose.tau, "w"))
  #g = g + scale_size_manual(values=popt$linesize)
  g = g + ylab("Conc (nM) or percent")
  gg = saveplot(5,5,dirs,"3cmtct_shed3_explore",draft.flag)
  grid.arrange(gg)
  

  
