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