#HEADER 
# Need this to get rxode working. I don't know why I keep having to run this.
Sys.setenv(PATH = paste(Sys.getenv("PATH"), "C:/RBuildTools/3.4/bin/",
                        "C:/RBuildTools/3.4/mingw_64/bin", sep = ";"))
Sys.setenv(BINPREF = "C:/RBuildTools/3.4/mingw_64/bin/")

# To be called at the top of every Rmd file. Initialization code and some useful constants.
suppressMessages(source("ams_initialize_script.R"))

# Note: Insoluble and soluble drugs are treated differently. 
# For insoluble, Tacc = M3totss/M30. For soluble, Tacc = S3totss/S30. 
# We are also interested in different parameters in the SA for either case.
# See writeup for details.

# Initialize ----

# Load model.
model = ivsc_4cmtct_shedct(target = TRUE)
# Drugs to explore. 
drugs = c("Atezolizumab", "Trastuzumab", "Pembrolizumab", "Bevacizumab")
# List of parameters of interest.
parameters = c("dose")
# names(parameters) = parameters
# Create vector of soluble drugs.
soluble = "Bevacizumab"

# Create paths to data files for each drug.
paths = NULL
for (i in 1:length(drugs)) {
  paths[i] = paste0("../data/ModelF_", drugs[i],"_Params.xlsx")
}

# Dose time, frequency, compartment, nominal dose
tmax = 26*7 #days
tau  = 21   #days
compartment = 2
dose.nmol = scale.mpk2nmol
joined = NULL
isSoluble = FALSE

# Iterate over all the drugs. -----------------
for (i in 1:length(drugs)){
  # Load parameters.
  param.as.double =  read.param.file(paths[i])
  df_param =  as.data.frame(t(param.as.double))
  
  # Check if drug has soluble target.
  if (drugs[i] %in% soluble) {
    # Change any non-shed parameters including M to S.
    parameters = sapply(parameters, 
                        function(x) {
                          if (!grepl("shed|dose",x)) {return(gsub("M","S", x))}
                          return(x)
                        })
    isSoluble = TRUE
  }
  
  # Set range for parameters of interest in SA.
  # Check which parameters are nonzero, not including dose which isn't in df_param.
  n.iter = 21
  params.to.iterate = data.frame(dose = lseq(scale.mpk2nmol*1e-5,scale.mpk2nmol*1e5,length.out = n.iter))

  dfs = list()
  # Iterate all of the parameters for a single drug.
  for (j in 1:ncol(params.to.iterate)){
    dfs[[j]] = compare.thy.sim(model = model,
                               param.as.double = param.as.double,
                               dose.nmol = dose.nmol,
                               tmax = tmax,
                               tau  = tau,
                               compartment = compartment,
                               param.to.change = names(params.to.iterate)[j],
                               param.to.change.range = params.to.iterate[[j]],
                               soluble = isSoluble)
    dfs[[j]] = dfs[[j]] %>% mutate(drug = drugs[i], isSol = isSoluble)
  }
  joined = bind_rows(joined,dfs)
  
  # Reset isSoluble to false since the default in compare.thy.sim is false.
  isSoluble = FALSE
}

# Plot the output AFIRT ----

# I only wanted a subset of the parameters for gather.
joined.AT = joined[c("fold.change.param","AFIRT.Kssd_", "AFIRT.Kss_", "AFIRT.Kd_", "AFIRT.sim", "TFIRT.Kssd_", "TFIRT.Kss_", "TFIRT.Kd_", "TFIRT.sim","param","drug", "isSol","M30","S30","Thiele")] %>%
  rename(Thiele.Modulus=Thiele)

names(joined.AT) = names(joined.AT) %>%
  str_replace("_","")

# Use gather to make the long data frame for ggplot.
plots.AT = joined.AT %>% 
  gather(key, AFIRT.value, -param, -fold.change.param,-drug, -isSol, -M30, -S30) %>%
  filter(!is.na(AFIRT.value)) %>%
  mutate(Target = ifelse(isSol,"Soluble","Membrane-Bd"),
         Target = factor(Target,levels=c("Soluble","Membrane-Bd")),
         AT     = NA,
         AT     = ifelse(str_detect(key,"AFIRT"),"AFTIR",AT),
         AT     = ifelse(str_detect(key,"TFIRT"),"TFTIR",AT),
         AT     = ifelse(str_detect(key,"Thiele"),"Thiele.Modulus",AT),
         key    = str_replace(key,"\\wFIRT.",""),
         key    = str_replace(key,"sim","simulation"),
         key    = str_replace(key,"K","theory.K"),
         key    = factor(key,levels=c("simulation","theory.Kssd","theory.Kss","theory.Kd","Thiele.Modulus")),
         T30    = ifelse(isSol,S30,M30),
         T30    = paste("T0 =",signif(T30,1),"nM")) %>%
  mutate(drug   = str_replace(drug,"Herceptin","Trastuzumab"),
         drug = factor(drug,levels=c("Bevacizumab","Pembrolizumab","Atezolizumab","Trastuzumab")))

plots.A = filter(plots.AT,AT=="AFTIR" | AT=="Thiele.Modulus")
plots.T = filter(plots.AT,AT=="TFTIR" | AT=="Thiele.Modulus")

Thiele.Thresh = 0.2
plot.Thiele.hline = plots.A %>%
  filter(AT=="AFTIR") %>%
  mutate(AT = "Thiele.Modulus") %>%
  filter(!duplicated(paste(param,drug)))

plot.Thiele.vline = plots.A %>%
  filter(AT == "Thiele.Modulus") %>%
  mutate(Thiele.Diff = abs(AFIRT.value-Thiele.Thresh)) %>%
  arrange(Thiele.Diff) %>%
  group_by(drug) %>%
  slice(1) %>%
  filter(!(fold.change.param==100)) 

plot.Thiele.vline.AFTIR = plot.Thiele.vline %>%
  mutate(AT = "AFTIR") %>%
  bind_rows(plot.Thiele.vline)

plot.Thiele.vline.TFTIR = plot.Thiele.vline %>%
  mutate(AT = "TFTIR") %>%
  bind_rows(plot.Thiele.vline)

plot.Thiele.Thresh.Label = plot.Thiele.hline %>%
  filter(drug=="Bevacizumab")
#----
dose.window.A = plots.A %>%
  select(AT,drug,Target,T30) %>%
  filter(!duplicated(paste0(AT,drug))) %>%
  mutate(xmin = 0.1, 
         xmax = 10, 
         ymin = 1e-7,
         ymax = c(1,1,1,1,1e7,1e7,1e7))

g = ggplot(plots.A, aes(fold.change.param, AFIRT.value, color=key, shape=key, linetype=key))
g = g + scale.x.log10(2,breaks = c(1e-4,1e-2,1,1e2,1e4),labels=c("1e-4","0.01","1","100","1e4"))
g = g + scale.y.log10(2)
g = g + geom_line(mapping=aes(size=key))
g = g + geom_hline(data=plot.Thiele.hline,aes(yintercept=Thiele.Thresh))
g = g + geom_label(data=plot.Thiele.Thresh.Label,aes(x=1,y=Thiele.Thresh,label=Thiele.Thresh,color=NULL,size=NULL),show.legend = FALSE)
g = g + facet_grid(AT~drug + Target + T30,scales = "free_y",switch = "y")
g = g + labs(x     = "Dose (mg/kg) every 3 weeks",
             y     = "Value",
             #y    = "Average Free Tissue target to\nInitial target Ratio* (AFTIR*)",
             color = "",
             shape = "",
             size  = "",
             linetype="")
g = g + theme(axis.title  = element_text(size=15),
              strip.text  = element_text(size=12),
              legend.text = element_text(size=12))
g = g + scale_color_manual(values = c(simulation = "black",
                                      theory.Kd  = "red",
                                      theory.Kss = "green4",
                                      theory.Kssd= "blue",
                                      Thiele.Modulus = "grey50"))
#g = g + scale_shape_manual(values = c(simulation = 46,
#                                      theory.Kd  = 15,
#                                      theory.Kss = 17,
#                                      theory.Kssd= 16))
g = g + scale_size_manual(values  = c(simulation = 1.5,
                                      theory.Kd  = 1.5,
                                      theory.Kss = 1.5,
                                      theory.Kssd= 1.5,
                                      Thiele.Modulus = 2))
g = g + scale_linetype_manual(values  = c(simulation = "solid",
                                      theory.Kd = "dotdash",
                                      theory.Kss = "dashed",
                                      theory.Kssd = "longdash",
                                      Thiele.Modulus = "solid"))
g.keep = g
g = g + geom_vline(data=plot.Thiele.vline.AFTIR,aes(xintercept=fold.change.param))
#g = g + geom_rect(data=dose.window.A,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
#                  fill = "grey50", alpha = .5, inherit.aes = FALSE)
print(g)
width = 10
height = 6
ggsave("../results/Task09f_AFTIR_Kssd_Kss_Kd_wThiele.pdf",width = width,height=height)

#----

g = g.keep %+% plots.T
g = g + geom_vline(data=plot.Thiele.vline.TFTIR,aes(xintercept=fold.change.param))
dose.window.T = dose.window.A %>%
  mutate(AT = ifelse(AT=="AFTIR","TFTIR",AT))
#g = g + geom_rect(data=dose.window.T,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
#                  fill = "grey50", alpha = .5, inherit.aes = FALSE)
#g = g + labs(y="Trough Free Tissue target to\ninitial target Ratio* (TFTIR*)")
print(g)
ggsave("../results/Task09f_TFTIR_Kssd_Kss_Kd_wThiele.pdf",width = width,height=height)

