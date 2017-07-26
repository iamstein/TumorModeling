#header
source("ams_initialize_script.R")
source("ivsc_2cmtc_keq.R")


dirs <- get.dirs(sys.calls(), dirs)

mod <- ivsc_2cmtc_keqct_CLV_T0()

# input parameter
pin <- param.import("../data/efalizumab_bauer99_keq_steinfit_param.csv")
p <- mod$repar(pin)[mod$pin] # reparametrization

p["T0"] <- 10

# test function
event <- eventTable()
event$add.sampling(c('1e-4', '1e-3', '1e-2', seq(.1, 6*7, .1)))
event$add.dosing(dose<-3*scale.mpk2nmol, nbr.doses<-1, dosing.to<-2)
out <- mod$rxode$solve(mod$repar(p), event, mod$init(p))

# specify ranges of parameters to exploreparam <- c("dose", "CL", "Q", "ksyn", "Keq", "T0")
units <- c("nmol", "L/d", "L/d", "nM/d", "nM", "nM")
nparam <- length(param)
order <- 1:nparam
title.scale <- c(rep(1, nparam))

explore <- data.frame(param=param, units=units, title.scale=title.scale, order=order,
                    fold.min=0.01, fold.max=100, fold.n=5, StringAsFactors=FALSE)





explore2 <- data.frame(param="T0", units="nM", fold.min=0.001, fold.max=10, fold.n=5,
                    StringAsFactors = FALSE)


# specify dosing and sampling

p["dose"] = 1*scale.mpk2nmol
p["tau"] = 1
event = list(
    t.sample = c(1e-4, 1e-3, 1e-2, seq(.1, 6*7, 0.1)), 
    n.dose = 1,
    cmt = which(mod$cmtnae="Amt.central")
)

# simulate model and put into OUT
# key variables are: ev-dosing, event matrix, p0-baseline parameter, model, ex=exploration space
# function(model, p0, ev,ex)

time1 = proc.time() # take the initial time

OUT = multi.solve(mod, p, event, explore, explore2)















