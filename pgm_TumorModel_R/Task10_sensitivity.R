# Setup.
setwd("~/GitHub/TumorModeling/pgm_TumorModel_R")
source("ams_initialize_script.R")
dirs = get.dirs(sys.calls(),dirs)
theme_set(theme_classic())

# Create fold change function.
lseq = function(from,to,length.out){exp(seq(log(from), log(to), length.out = length.out))}
scale.nmol2mpk = 1e-9*150e3/1000*70 #nM->mg/kg
scale.mpk2nmol = 1/scale.nmol2mpk

# Load data and model.
mod      = ivsc_3cmtct_full()
pin        = xlsx::read.xlsx("../data/Bx_DTN_example.xlsx",1)
# Assumption change value.
pin[14:17,3] = 0

#'1/(nM*d)'n) Create parameter vector.
p = pin$Value
names(p) = pin$Parameter

# Specify ranges of parameters to explore
# param   = c('kon3', 'koff3','dose')
# units   = c('1/(nM*d)','1/(nM*d)','mg')
param   = c('dose')
units   = c('mg')
nparam  = length(param)
order   = 1:nparam             #order for parameters to be plotted
title.scale = c(rep(1,nparam)) #in case one wants to change units in title

# Setup and exploration data.frame to iterate through for each of the analysis to perform.
explore = data.frame(param = param, units=units, title.scale=title.scale, order=1:nparam,
fold.min = .001, fold.max = 10, fold.n   = 6,
stringsAsFactors = FALSE)
explore = explore[explore$order,]
row.names(explore) = param
  
# Specify dosing and sampling. I am pretty sure I rewrite all of these values later when I set up the event table.
p["dose"]   = 420 #mg
p["tau"]    = 1   #dosing interval
ev = list(
t.sample = c(1e-4,1e-3,1e-2,seq(.1,6*7,.1)), #d
n.dose   = 1,
# cmt      = which(mod$cmtname=="Amt.central")
cmt      = which(mod$cmtshort=="D1")
)

# Simulate model and put into OUT
# key variables are: ev-dosing event matrix, p0-baseline parameters, model, ex=exploration space
# function(model,p0,ev,ex)
p0    = p
OUT   = data.frame()

for (par in explore$param) {
  ex = explore[par,]
  foldchange = lseq(ex$fold.min,ex$fold.max,ex$fold.n)
  for (foldpar in foldchange) {
    # Set up parameters
    p      = p0
    p[par] = p0[par]*foldpar
    T_30 = p["koff3"]/p["kon3"]
    # Set up dosing (where dose and dosing interval are parameters that can vary)
    ndose     = 1
    event     = eventTable(amount.units="nmol", time.units = "days")
    event$add.sampling(unique(sort(c(seq(0,13*7,.1),10^(-9:0)))))
    event$add.dosing(dose=p["dose"]/ndose*scale.mpk2nmol,nbr.doses=4*ndose,dosing.interval=21/ndose,dosing.to=2)
    
    # Solve the ODEs.
    out      = mod$rxode$solve(p, event, mod$init(p))
    out      = mod$rxout(out,p)
 
    #add useful key quantities for plotting
    out$pchg = par  # which parameter changed
    # out$pval = p0[par] # I think this should tell me what the parameter changed to but it isn't working that way.
    out$pval = p[par]
    out$pfold= foldpar
    out$order= ex$order
    out$label.group= paste0(par," (",ex$units,")\n",
                            signif(ex$title.scale*p0[par]*ex$fold.min,2),"-",signif(ex$title.scale*p0[par]*ex$fold.max,2))
    # out$label.group= factor(out$label.group,levels=unique(out$label.group)) # this bit fixes the order for ggplot
    out$label.ind  = paste(par,signif(p[par]))
    out$ratio = out$T3/T_30
    
    # Sanity check. 
    # g = ggplot(data=out,aes(x=time,y=ratio)) + labs(title = p[par]) + geom_line() + scale_y_log10() + xlab("day") 
    # print(g)
    
    # Append to data.frame()
    OUT = rbind(OUT,out)
  }
    # Sensitivity Plot
    g = ggplot(data=OUT,aes(x=time,y=Dtot1,group=pval, color=pval)) + geom_line() + facet_grid(~label.group) + 
        scale_y_log10() + scale_colour_gradient2(trans="log",name = par,guide="legend",breaks=foldchange*p0[par], low="blue",high="red" )
    print(g)
}
