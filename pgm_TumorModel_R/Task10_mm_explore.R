source("ams_initialize_script.R")
source("ivsc_2cmtc_mm.R")
dirs = get.dirs(sys.calls(),dirs)

theme_set(theme_classic())

scale.nmol2mpk = 1e-9*150e3/1000*70 #nM->mg/kg
scale.mpk2nmol = 1/scale.nmol2mpk

mod      = ivsc_2cmtc_mm_CLV()
pin      = param.import("../data/panitumumab_Ma09_param.csv")
p        = mod$repar(pin)[mod$pin] #remove extra parameters (helps with redudancy and debugging)

#specify ranges of parameters to explore
  param   = c('dose' ,'CL' ,'Vm'    ,'Km')
  units   = c('mg'   ,'L/d','nmol/d','nM')
  nparam  = length(param)
  order   = 1:nparam             #order for parameters to be plotted
  title.scale = c(rep(1,nparam)) #in case one wants to change units in title

  explore = data.frame(param = param, units=units, title.scale=title.scale, order=1:nparam,
    fold.min = .1, fold.max = 10, fold.n   = 7,
    stringsAsFactors = FALSE)
  explore = explore[explore$order,]
  row.names(explore) = param
  
#specify dosing and sampling
  p["dose"]   = 420 #mg
  p["tau"]    = 1   #dosing interval
  ev = list(
    t.sample = c(1e-4,1e-3,1e-2,seq(.1,6*7,.1)), #d
    n.dose   = 1,
    cmt      = which(mod$cmtname=="Amt.central")
  )
  
#simulate model and put into OUT
#key variables are: ev-dosing event matrix, p0-baseline parameters, model, ex=exploration space
#function(model,p0,ev,ex)
time1 = proc.time()
p0    = p
OUT   = data.frame()
lseq = function(from,to,length.out){exp(seq(log(from), log(to), length.out = length.out))}
for (par in explore$param) {
  ex = explore[par,]
  for (foldpar in lseq(ex$fold.min,ex$fold.max,ex$fold.n)) {
    #set up parameters
    p      = p0
    p[par] = p0[par]*foldpar
    
    #set up dosing (where dose and dosing interval are parameters that can vary)
    event  = eventTable()
    event$add.sampling(ev$t.sample)
    event$add.dosing(dose=p["dose"], nbr.doses=ev$n.dose, dosing.interval=p["tau"], dosing.to=ev$cmt)
  
    #solve ode
    out      = mod$rxode$solve(mod$repar(p), event, mod$init(p))
    out      = mod$rxout(out,p)
  
    #compute Cnonlin for plotting
    Cnonlin = mod$fun.Cnonlin(p)
    tnonlin = approx(x=out$D,y=out$time,xout=Cnonlin)
    out$Cnonlin = Cnonlin
    out$tnonlin = tnonlin$y
    
    #dose normalized
    out$Dnorm = out$D/p["dose"]
    out$Cnonlin.norm = Cnonlin/p["dose"]
    
    #add useful key quantities for plotting
    out$pchg = par
    out$pval = p0[par]
    out$pfold= foldpar
    out$order= ex$order
    out$label.group= paste0(par," (",ex$units,")\n",
                            signif(ex$title.scale*p0[par]*ex$fold.min,2),"-",signif(ex$title.scale*p0[par]*ex$fold.max,2))
    out$label.group= factor(out$label.group,levels=unique(out$label.group)) #this bit fixes the order for ggplot
    out$label.ind  = paste(par,signif(p[par]))
    
    #create giant output matrix
    if   (length(OUT)==0) OUT=out 
    else OUT = rbind(OUT,out)  
  }
  
}

time2 = proc.time()

#plot data ----
g = ggplot(data=OUT,aes(x=time,y=Dnorm,group=label.ind,color=pfold))
g = g + facet_grid(~label.group)
g = g + geom_line()
g = g + scale_colour_gradient2(trans="log",low="blue",mid="grey50",high="red",guide="legend",breaks=c(10,1,.1))
g = g + scale_y_log10() 
g = g + geom_point(aes(x=tnonlin,y=Cnonlin.norm))
g = g + xlab("day") + ylab("Conc-Norm (ug/ml)/mg")
g = xscale("w6")
gg = saveplot(7,3,dirs,"pani_mm_CLV",draft.flag)
grid.arrange(gg)
time3 = proc.time()

print("code time =")
print(unname(time2[3]-time1[3]))
print("plot time =")
print(unname(time3[3]-time2[3]))


