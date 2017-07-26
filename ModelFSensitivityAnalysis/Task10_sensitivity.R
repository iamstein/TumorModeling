# Setup.
source("ams_initialize_script.R")
source("ivsc_4cmtct_shedct.R")
dirs = get.dirs(sys.calls(),dirs)
theme_set(theme_classic())

# Create fold change function.
lseq = function(from,to,length.out){exp(seq(log(from), log(to), length.out = length.out))}
scale.nmol2mpk = 1e-9*150e3/1000*70 #nM->mg/kg
scale.mpk2nmol = 1/scale.nmol2mpk

# Load data and model.
mod      = ivsc_4cmtct_shedct()
pin        = xlsx::read.xlsx("ivsc_4cmtct_shedct_param.xlsx",1)
# Assumption change value.

#'1/(nM*d)'n) Create parameter vector.
p = pin$Value
names(p) = pin$Parameter

# Specify ranges of parameters to explore

param   = c('dose','kon3')
units   = c('mg/kg', 'L/d')
nparam  = length(param)
order   = 1:nparam             #order for parameters to be plotted
title.scale = c(rep(1,nparam)) #in case one wants to change units in title

# Setup and exploration data.frame to iterate through for each of the analysis to perform.
explore = data.frame(param = param, units=units, title.scale=title.scale, order=1:nparam,
fold.min = 0.01, fold.max = 1, fold.n = 6,
stringsAsFactors = FALSE)

explore = explore[explore$order,]
row.names(explore) = param
  

p["dose"]   = 42 #mg/kg
p["kon3"]    = 100  
ev = list(
t.sample = c(1e-4,1e-3,1e-2,seq(.1,6*7,.1)), #d
n.dose   = 1,
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
  
    # Set up dosing (where dose and dosing interval are parameters that can vary)
    ndose     = 1
    event     = eventTable(amount.units="nmol", time.units = "days")
    event$add.sampling(unique(sort(c(seq(0,13*7,.1),10^(-9:0)))))
    event$add.dosing(dose=p["dose"]/ndose*scale.mpk2nmol,nbr.doses=8*ndose,dosing.interval=12/ndose,dosing.to=2)
  
    
    # Solve the ODEs.
    out      = mod$rxode$solve(p, event, mod$init(p))
    out      = mod$rxout(out,p)
 
    #add useful key quantities for plotting
    out$pchg = par  # which parameter changed

    out$pval = p[par]
    out$pfold= foldpar
    out$label.group= paste0(par," (",ex$units,")\n",
                            signif(ex$title.scale*p0[par]*ex$fold.min,2),"-",signif(ex$title.scale*p0[par]*ex$fold.max,2))

    out$label.ind  = paste(par,signif(p[par]))
    
    OUT = rbind(OUT,out)
  }
    # Sensitivity Plot
  
    g = ggplot(data=OUT,aes(x=time, y=Dtot3, group=pval, color=pval))
    g = g + geom_line() 
    g = g + facet_grid(~label.group) 
    g = g + scale_y_log10() 
#    g = g + scale_colour_gradient2(trans="log",name = par,guide="legend",breaks=foldchange*p0[par], low="blue",high="red" )
    print(g)
    
}
ggsave("SensitivityAnalysis.png", plot=g)