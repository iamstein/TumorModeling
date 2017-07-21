#header
source("ams_initialize_script.R")
source("ivsc_2cmtc_keq.R")
dirs = get.dirs(sys.calls(),dirs)

#get model and parameters    
  mod      = ivsc_2cmtc_keqct_CLV_T0()
  pin      = param.import("../data/efalizumab_bauer99_keq_steinfit_param.csv")
  p        = mod$repar(pin)[mod$pin] #remove extra parameters (helps with redudancy and debugging)
  p["T0"]  = 10
  
#test function
  event    = eventTable()
  event$add.sampling(c(1e-4,1e-3,1e-2,seq(.1,6*7,.1)))
  event$add.dosing(dose=3*scale.mpk2nmol, nbr.doses=1, dosing.to=2)
  out      = mod$rxode$solve(mod$repar(p), event, mod$init(p))

#specify ranges of parameters to explore
  param   = c('dose' ,'CL' ,'Q'  ,'ksyn','Keq','T0')
  units   = c('nmol' ,'L/d','L/d','nM/d','nM' ,'nM')
  nparam  = length(param)
  order   = 1:nparam             #order for parameters to be plotted
  title.scale = c(rep(1,nparam)) #in case one wants to change units in title
  
  explore = data.frame(param = param, units=units, title.scale=title.scale, order=1:nparam,
                       fold.min = .01, fold.max = 100, fold.n   = 5,
                       stringsAsFactors = FALSE)
  explore = explore[explore$order,]
  row.names(explore) = param
  
  explore2= data.frame(param="T0",units="nM",fold.min=.001,fold.max=10,fold.n=5,stringsAsFactors = FALSE)
  
#specify dosing and sampling
  p["dose"]   = 1*scale.mpk2nmol #nmol
  p["tau"]    = 1   #dosing interval
  event = list(
    t.sample = c(1e-4,1e-3,1e-2,seq(.1,6*7,.1)), #d
    n.dose   = 1,
    cmt      = which(mod$cmtname=="Amt.central")
  )
  
#simulate model and put into OUT ----
#key variables are: ev-dosing event matrix, p0-baseline parameters, model, ex=exploration space
#function(model,p0,ev,ex)
time1 = proc.time()
OUT   = multi.solve(mod,p,event,explore,explore2)
OUT1  = OUT[!duplicated(OUT$label.ind),]
time2 = proc.time()


#OUTs  = select(OUT,time,D,Dnorm,pchg,label.group,label.group2,label.ind,pfold)
#OUTm  = melt(OUTs,id.vars = c("time","pchg","label.group","label.group2","label.ind","pfold"))

OUT$Dnorm[OUT$Dnorm<1e-4] = NA
OUT$D[OUT$D<1e-5] = NA
OUT$T0 = OUT$pval2.ex2

#plot result ----
g = ggplot(data=OUT,aes(x=time,y=D,group=label.ind,color=pfold))
g = g + facet_grid(label.group2~label.group,scales="free_y")
g = g + geom_line()
g = g + scale_colour_gradient2(trans="log",low="blue",mid="grey50",high="red",guide="legend",breaks=10^signif(-2:2))
g = g + scale_y_log10()
g = g + geom_hline(aes(yintercept=T0,linetype="dotted"),linetype="dotted")
#g = g + scale_y_log10(lim=c(1e-5,1)) 
#g = g + geom_point(aes(x=tnonlin,y=Cnonlin.norm))
g = g + xlab("day") #+ ylab("Dnorm (nM/nmol)")

g = xscale("w6")
gg = saveplot(8,8,dirs,"keqct_explore2d_D",draft.flag)
grid.arrange(gg)
time3 = proc.time()

print("code time =")
print(unname(time2[3]-time1[3]))
print("plot time =")
print(unname(time3[3]-time2[3]))

# add dnorm
gnorm = g + aes(y=Dnorm)
gnorm$layers[[2]] = NULL
gg = saveplot(8,8,dirs,"keqct_explore2d_Dnorm",draft.flag)
grid.arrange(gg)

