# Setup.
#setwd("~/TumorModeling/pgm_TumorModel_R")
source("ams_initialize_script.R")
source("ivsc_4cmtct_shedct.R")
dirs = get.dirs(sys.calls(),dirs)

# Functions
lseq = function(from,to,length.out){exp(seq(log(from), log(to), length.out = length.out))}

# Creates readable version of variable names from input data.frame
var2nice = function(df) {
  df = df %>%  mutate(
         varnice = variable,
         varnice = str_replace_all(varnice,"1",".pla"), # plasma
         varnice = str_replace_all(varnice,"3",".tum"), # tumor
         varnice = str_replace_all(varnice,"D","Drug"),    
         varnice = str_replace_all(varnice,"S","s-Target"),
         varnice = str_replace_all(varnice,"M","m-Target"),
         varnice = str_replace_all(varnice,".pct","%"),
         molecule= paste0(str_extract(varnice,"^\\w+."),"[nM]")
      )
  df$molecule[str_detect(df$varnice,"free")]="free.ratio"
  df$molecule = factor(df$molecule,levels=unique(df$molecule))
  return(df)
}

# Model and parameters. 
mod      = ivsc_4cmtct_shedct()
d        = xlsx::read.xlsx("../data/ivsc_4cmtct_shedct_param.xlsx",1,stringsAsFactors=FALSE)
d        = d %>% filter(!is.na(Parameter))
param    = d
row.names(param) = param$Parameter
p        = d$Value
names(p) = d$Parameter

# Initial dosage and range of dosages to explore.
dose0 = 80
foldchange = lseq(0.01,100,20)
dose.range = dose0*foldchange

dose.tau= 14 # dose every tau weeks
n.dose = 600 
tmax   = 28*3

OUT = NULL
for (dose.mg in dose.range) {  
  # Create event table for each dose.
  ev       = eventTable(amount.units="nmol", time.units = "days")
  ev$add.sampling(unique(sort(c(seq(-7,tmax,.1),10^(-3:0)))))
  ev$add.dosing(dose=dose.mg*scale.mg2nmol,nbr.doses=floor(tmax/dose.tau)+1,dosing.interval=dose.tau,dosing.to=2)
  
  # Save initial parameter values in init.
  init        = mod$init(p)
  
  # Solve the ODE.
  out         = mod$rxode$solve(p, ev, mod$init(p))
  out         = mod$rxout(out,p)
  
  # Calculate simulated S0/S1,0 and M3/M3,0.
  out = out %>%
    mutate(Sfree.pct = S1/init["S1"],
           Mfree.pct = M3/init["M3"],
           dose.mg   = dose.mg,
           dose.nM   = dose.mg/1000/150e3*1e9) # mg-->g-->mol-->nMol
  OUT = rbind(OUT,out)
}

meanMfree = OUT %>% 
  filter(time>28*2 & time<=28*3 ) %>%
  group_by(dose.mg) %>% 
  summarize(Mavg = mean(Mfree.pct),
            dose.nM = dose.nM[1]) %>%
  ungroup()

d = OUT %>%
  dplyr::select(time,D1,D3,S1,S3,Stot1,Stot3,M1,Mtot1,M3,Mtot3,Sfree.pct,Mfree.pct,dose.mg,dose.nM) %>%
  gather(variable,value,-c(time,dose.mg,dose.nM)) %>%
  mutate(value = ifelse(value>1e-8,value,NA))

# Features for output graphs.
popt      = data.frame(variable = c("D1"   ,"D3"    ,"S1"   ,"S3" ,"M1"   ,"M3"    ,"Stot1","Stot3","Mtot1"  ,"Mtot3" ,"Sfree.pct" ,"Mfree.pct"),
                       color    = c("black","gray50","red"  ,"orange","blue"  ,"red4" ,"orange4","blue4" ,"pink"       ,"lightblue", "green", "brown"),
                       linetype = c("solid","solid" ,"solid","solid" ,"solid" ,"solid","solid"  ,"solid" ,"dashed"     ,"dashed", "solid", "solid"),
                       stringsAsFactors = FALSE)

d    = var2nice(d)
popt = var2nice(popt)

# Color and linetypes for graphs.
color     = popt$color
names(color) = popt$varnice
linetype  = popt$linetype
names(linetype) = popt$varnice

# Calculate constants necessary for AFIRT calculations.
pp  = as.data.frame(t(p))
CL  = pp$keD1*pp$VD1 #L/d
tau = dose.tau #days
Kd  = pp$koff3/pp$kon3
kss = with(pp, ( koff3 + keDM3 + kshedM3 ) / kon3 )

Mtot3ss = with(pp,
          (k13DM*ksynM1*(VD1/VD3)+(keDM1+kshedM1+k13DM)*ksynM3) /
          ((keDM1+kshedM1+k13DM)*(keDM3+kshedM3+k31DM)-k31DM*k13DM))
M03     = with(pp, 
          ((k13M*ksynM1*(VD1/VD3)+(keM1+kshedM1+k13M)*ksynM3)/
          ((keM1+kshedM1+k13M)*(keM3+kshedM3+k31M)-k31M*k13M)))
Tacc.tum = Mtot3ss/M03
B = with(pp, ( k13D*VD1/VD3 ) / ( keD3 + k31D ))

dose = d %>%
  subset(!duplicated(dose.mg)) %>%
  dplyr::select(dose.nM,dose.mg)

# AFIRT calculations.
AFIRT = d %>%
  filter(varnice == "m-Targetfree%",time>28*2 & time<=28*3) %>%
  group_by(dose.nM) %>%
  summarise(AFIRT.sim = mean(value),
            dose.mg   = dose.mg[1]) %>%
  ungroup() %>%
  mutate(AFIRT.theory.ss = kss*Tacc.tum*(CL*tau)/(B*dose.nM), 
        AFIRT.theory.d= Kd*Tacc.tum*(CL*tau)/(B*dose.nM))

# Gather makes a longer dataframe.
AFIRTlong = AFIRT %>%
  gather(key,value,-c(dose.mg,dose.nM))

# Calculate B simulated for graph.
B.sim = OUT %>%
  filter(time>28*2 & time<=28*3) %>%
  group_by(dose.nM) %>%
  summarise(B.sim = mean(Dtot3/Dtot1),
            dose.mg   = dose.mg[1],
            Dtotavg1.sim = mean(Dtot1)) %>%
  ungroup() %>%
  mutate(B.theory = B,
         Dtotavg1.theory = dose.nM/(pp$keD1*pp$VD1*dose.tau))

# Graphing.
theme = theme(text = element_text(size = 16), legend.title=element_blank())

g = ggplot(AFIRTlong, aes(dose.mg,value,color=key,shape=key)) +
  scale.x.log10()+scale.y.log10()+ 
  geom_line() +
  geom_point() +
  theme + 
  ylab("AFIRT") + 
  xlab("dose (mg)") + 
  guides(shape=FALSE) + 
  scale_color_manual(values = c("red", "green", "blue"), labels = c("Simulated","Theoretical,Kd", "Theoretical,Kss"))
print(g)

h = ggplot(B.sim, aes(dose.mg,B.sim)) +
  scale.x.log10()+
  geom_point(size=3) +
  geom_line(aes(y=B.theory)) + 
  theme
print(h)

h = ggplot(OUT,aes(time,D1,group=dose.mg))+geom_line()
print(h)