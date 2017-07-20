source("ams_initialize_script.R")
source("ivsc_2cmtc_full.R")
dirs = get.dirs(sys.calls(),dirs)

mod      = ivsc_2cmtc_full()
p        = param.import("../data/siltuximab_mayer15_steintweak_param.csv")

#these functions aren't strictly needed here.  just testing because will be useful later
p        = p[mod$pin]
p        = mod$repar(p)

#read data
d        = read.csv("../data/Puchalski10_Pooled.csv",stringsAsFactors = FALSE)
ddose    = d[!duplicated(d$dose_nmol),]

#simulate model and put into OUT
for (i in 1:nrow(ddose)) {
  ev       = eventTable(amount.units = "nmol", time.units = "days")
  ev$add.sampling(seq(0,13*7,.1))
  ev$add.dosing(dose = ddose$dose_nmol[i],nbr.doses = 4, dosing.interval = 21, dosing.to = 2)
  
  out      = mod$rxode$solve(p, ev, mod$init(p))
  out      = mod$rxout(out,p)
  out$dose_mpk = paste(ddose$dose_mpk[i], "mg/kg")
  if (i==1) OUT=out else OUT = rbind(OUT,out)  
}

#create a melted dataset for plotting using ggplot
OUT         = select(OUT,time,dose_mpk,D,T,DT,Dtot,Ttot)
OUTm        = melt(OUT,id.vars = c("time","dose_mpk"))
OUTm$dvalue = NA
dm          = transmute(d,time=t, dose_mpk=paste(dose_mpk,"mg/kg"), variable=type, value=NA, dvalue=y)
allm        = rbind(OUTm,dm)

#set some plotting options
popt      = data.frame(variable = c("D"         ,"T"    ,"DT"    ,"Dtot"  ,"Ttot"),
                       color    = c("blue"      ,"red"  ,"purple","blue4" ,"red4"),
                       linetype = c("solid"     ,"solid","solid" ,"dashed","dashed"),
                       linesize = c( .5         , .5    , .5    , 1       , 1    ),
                       stringsAsFactors = FALSE)

#plot data
g = ggplot(data=allm,aes(x=time,y=value,group=variable,color=variable))
g = g + facet_grid(~dose_mpk)
g = g + geom_line(aes(linetype=variable,size=variable))
g = g + geom_point(aes(y=dvalue),size=3)
g = g + scale_y_log10()
g = g + scale_color_manual(   values=popt$color)
g = g + scale_linetype_manual(values=popt$linetype)
g = g + scale_size_manual(values=popt$linesize)
g = g + xlab("day") + ylab("Conc (nM)")

gg = saveplot(7,3,dirs,"full",draft.flag)
grid.arrange(gg)
