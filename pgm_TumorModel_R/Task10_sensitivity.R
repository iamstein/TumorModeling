# Setup.
setwd("~/TumorModeling/pgm_TumorModel_R")
source("ams_initialize_script.R")
dirs = get.dirs(sys.calls(),dirs)
theme_set(theme_classic())

# Create fold change function.
lseq = function(from,to,length.out){exp(seq(log(from), log(to), length.out = length.out))}

# Load data and model.
mod = ivsc_3cmtct_full()
pin = xlsx::read.xlsx("../data/Bx_DTN_example.xlsx",1)
# Assumption made for model C. Comment out if model changes to F.
pin[14:17,3] = 0

# Create parameter vector.
p = pin$Value
names(p) = pin$Parameter

# Specify parameters to explore
param   = c('kon3', 'koff3','dose','keD1','k13D','keDT3')
units   = c('1/(nM*d)','1/(nM*d)','mg','1/d','1/d','1/d')
nparam  = length(param)
order   = 1:nparam             # order for parameters to be plotted
title.scale = c(rep(1,nparam)) # in case one wants to change units in title

# Setup and exploration data.frame to iterate through for each of the analysis to perform.
explore = data.frame(param = param, units=units, title.scale=title.scale, order=1:nparam,
  fold.min = .01, fold.max = 10, fold.n   = 5,
  stringsAsFactors = FALSE)
explore = explore[explore$order,]
row.names(explore) = param
  
# Specify dosing and sampling.
p["dose"]   = 420 #mg
p["tau"]    = 1   #dosing interval
ev = list(
  t.sample = c(1e-4,1e-3,1e-2,seq(.1,6*7,.1)), #d
  n.dose   = 1,
  cmt      = which(mod$cmtshort=="D1")
)

# Simulate model and put into OUT
# key variables are: ev-dosing event matrix, p0-baseline parameters, model, ex=exploration space
p0    = p
OUT   = data.frame()

for (par in explore$param) {
  # Changes which parameter is currently being explored.
  ex = explore[par,]
  foldchange = lseq(ex$fold.min,ex$fold.max,ex$fold.n)
  for (foldpar in foldchange) {
    # Set up parameters
    # Set up local verion of parameter list for exploration.
    p      = p0
    p[par] = p0[par]*foldpar
    T_30 = p["koff3"]/p["kon3"]
    
    # Set up dosing (where dose and dosing interval are parameters that can vary)
    ndose     = 1
    event     = eventTable(amount.units="nmol", time.units = "days")
    event$add.sampling(unique(sort(c(seq(0,13*7,.1),10^(-9:0)))))
    event$add.dosing(dose=p["dose"]/ndose*scale.mpk2nmol,nbr.doses=4*ndose,
                     dosing.interval=21/ndose,dosing.to=2)
    
    # Solve the ODEs using RxODE.
    out      = mod$rxode$solve(p, event, mod$init(p))
    out      = mod$rxout(out,p)
 
    # Add useful key quantities for plotting.
    out$pchg = par  # which parameter changed
    out$pval = p[par]
    out$pfold= foldpar
    out$order= ex$order
    out$label.group= paste0(par," (",ex$units,")\n",
                            signif(ex$title.scale*p0[par]*ex$fold.min,2),"-",
                            signif(ex$title.scale*p0[par]*ex$fold.max,2))
    out$label.ind  = paste(par,signif(p[par]))
    out$ratio = out$T3/T_30
    
    # Append to data.frame()
    OUT = rbind(OUT,out)
  } # end param[i]

}# end param list

# Graphing.
theme = theme(text = element_text(size = 16))

# Vector of which variables from OUT you would like to be plotted.
plotme = c("Dtot1", "Ttot1", "ratio")

for (feature in plotme){
    # aes_string needed to do this in a loop
    print(ggplot(data=OUT,aes_string(x="time",y=feature,group="label.ind", color="pfold")) + 
        geom_line() + facet_grid(~label.group) + scale_y_log10() + theme + 
        scale_colour_gradient2(trans="log", guide="legend",low="blue", mid = "purple", high="red",breaks=c(10,1,.1)))
}