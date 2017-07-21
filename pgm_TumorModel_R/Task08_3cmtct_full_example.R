source("ams_initialize_script.R")
dirs = get.dirs(sys.calls(),dirs)

mod      = ivsc_3cmtct_full()
d        = xlsx::read.xlsx("../data/Bx_DTN_example.xlsx",1) abc
p        = d$Value
names(p) = d$Parameter

#these functions aren't strictly needed here.  just testing because will be useful later
  p        = p[mod$pin]
  p        = mod$repar(p)
  
  p["keT3"]= p["keT1"]
  
  print(t(data.frame(t(p))))

#simulate model and put into OUT
  ndose     = 1
  ev       = eventTable(amount.units="nmol", time.units = "days")
  ev$add.sampling(unique(sort(c(seq(-7,13*7,.1),10^(-9:0)))))
  ev$add.dosing(dose=300000000/ndose*scale.mpk2nmol,nbr.doses=4*ndose,dosing.interval=21/ndose,dosing.to=2)
  
  out         = mod$rxode$solve(p, ev, mod$init(p))
  out         = mod$rxout(out,p)
  out = out %>%
    mutate(ABCfree = D3/D1,
           ABCtot  = Dtot3/Dtot1)

#plot data ----
d = out %>%
  select(time,D1,D3,T1,T3,DT1,DT3,ABCtot) %>%
  gather(variable,value,-time) %>%
  mutate(value = ifelse(value>1e-8,value,NA))
  
popt      = data.frame(variable = c("D1","D3","T1","T3","DT1","DT3","ABCtot"),
                       color    = c("blue","blue4","red","red4","purple","purple4","black"),
                       linetype = c("solid","dashed","solid","dashed","solid","dashed","solid"),
                       stringsAsFactors = FALSE)

g = ggplot(d,aes(x=time,y=value,group=variable,color=variable))
  #g = g + facet_grid(~dose_mpk)
  g = g + geom_line(size=1)
  #g = g + geom_line(aes(linetype=variable))
  #g = g + geom_point(aes(y=dvalue),size=3)
  g = g + scale.y.log10()
  g = g + scale.y.log10(limits=c())
  g = g + scale_color_manual(   values=popt$color)
  g = xscale("d100",increment=7,t.start=-7)
  #g = g + scale_linetype_manual(values=popt$linetype)
  #g = g + scale_size_manual(values=popt$linesize)
  g = g + xlab("day") + ylab("Conc (nM)")
  gg = saveplot(7,3,dirs,"full",draft.flag)
  grid.arrange(gg)

# check formula for ABC ----
  #no elimination from tissue
  ABCthry.noTissElim = with(as.data.frame(t(p)),k13D*VD1/(k31D*VD3)) #doesn't account for elimination
  
  #elimination from tissue
  K    = with(as.data.frame(t(p))
              ,matrix(c(keD1+k13D    , -k13D*VD1/VD3,
                       -k31D*VD3/VD1 ,  keD3+k31D), nrow = 2))
  x    = solve(K,c(1,0))
  ABCthry.TissElim =  x[2]/x[1]
  
  ABCthry.TissElim2 = with(as.data.frame(t(p)),
                           k13D*VD1/VD3/(keD3+k31D))
  
  #numerics
  ABCnum = out %>%
    filter(time>21 & time<70) %>%
    summarize(mean(ABCtot)) %>%
    unlist()
  
  ABC_sum = data.frame(ABC.theory.noTissElim = ABCthry.noTissElim,
                       ABC.theory.TissElim   = ABCthry.TissElim,
                       ABC.theory.TissElim2  = ABCthry.TissElim2,
                       ABC.numerics = ABCnum)
  
  print(t(ABC_sum))
  
#check the formula for Tacc3 ----
  T03.thry   = with(as.data.frame(t(p)),
                    (ksyn1*VT1 /VT3 *k13T  + ksyn3*(keT1  + k13T )) / 
                    (keT1 *keT3  + keT1*k31T   + keT3*k13T))

  Ttot3.thry   = with(as.data.frame(t(p)),
                    (ksyn1*VDT1 /VDT3 *k13DT  + ksyn3*(keDT1  + k13DT )) / 
                      (keDT1 *keDT3  + keDT1*k31DT   + keDT3*k13DT))
  
    
  Tacc3.thry = with(as.data.frame(t(p)),
                   (ksyn1*VDT1/VDT3*k13DT + ksyn3*(keDT1 + k13DT)) / 
                   (ksyn1*VT1 /VT3 *k13T  + ksyn3*(keT1  + k13T )) *
                   (keT1 *keT3  + keT1*k31T   + keT3*k13T) /
                   (keDT1*keDT3 + keDT1*k31DT + keDT3*k13DT))
  
  T03.num    = out$T3[1]
  Ttot3.num  = out %>%
    filter(time>21 & time < 70) %>%
    summarize(mean(Ttot3)) %>%
    unlist()
  Ttot3_sum = data.frame(Ttot3.theory = Ttot3.thry,
                         Ttot3.numeric= Ttot3.num)
  print(t(Ttot3_sum))
  
  
  Tacc3.num = Ttot3.num/T03.num
  Tacc3_sum = data.frame(Tacc3.theory = Tacc3.thry,
                         Tacc3.numeric= Tacc3.num)
  print(t(Tacc3_sum))
  