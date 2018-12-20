#HEADER 

# To be called at the top of every Rmd file. Initialization code and some useful constants.
suppressMessages(source("ams_initialize_script.R"))

# Note: Insoluble and soluble drugs are treated differently. 
# For insoluble, Tacc = M3totss/M30. For soluble, Tacc = S3totss/S30. 
# We are also interested in different parameters in the SA for either case.
# See writeup for details.

# Initialize ----

# Load model.
model = ivsc_4cmtct_shedct(target = TRUE)
drugs = c("Atezolizumab", "Trastuzumab", "Pembrolizumab", "Bevacizumab")


# Create paths to data files for each drug.
param = NULL
for (i in 1:length(drugs)) {
  filename= paste0("../data/ModelF_", drugs[i],"_Params.xlsx")
  param[[i]] = read_excel(filename, 1) %>%
    mutate(Drug = drugs[i],
           Value = as.numeric(Value)) %>%
    select(Order,Drug,ParamType,Molecule,Description,Parameter,Value,Units) %>%
    filter(Parameter %in% model$pin) %>%
    filter(!(Parameter %in% c("F","ka")))
}
d = bind_rows(param) %>%
  select(Drug,Parameter,Value) %>%
  mutate(Value = as.character(signif(Value,2))) %>%
  spread(Drug,Value)

units = param[[1]] %>%
  select(Parameter,Units)
d = left_join(d,units,by="Parameter") %>% 
  select(Parameter,Units,everything())

file.prefix = paste0("../results/Task14_ParamTable")
readr::write_csv(d,paste0(file.prefix,".csv")) #to csv file
xtab = xtable::xtable(d) #to latex file
print(xtab,file=paste0(file.prefix,".tex"),include.rownames=FALSE,floating = FALSE)
kable(d)

