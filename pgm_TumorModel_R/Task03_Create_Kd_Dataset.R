source("ams_initialize_script.R")
dirs = get.dirs(sys.calls(),dirs)

#read in the data file
din = read_excel("../data/Bx_PK_Binding.xlsx",1)
d   = din %>%
  arrange(drug,param) %>%
  mutate(value=as.numeric(value)) %>%
  filter(!str_detect(param,"_custom"))  

source("ams_units.R")

d = d %>%
  mutate(param.unit = paste(param,valunit))

#explore parameters
x = filter(d,param %in% c("kon","koff"))

#plot parameters ----
g = ggplot(x,aes(x=value))
g = g + geom_histogram(bins=10)
g = g + facet_wrap(~param.unit,scales="free")
g = g + scale.x.log10()
gg = saveplot(6,4,dirs,"kon_koff_hist",draft.flag)
grid.arrange(gg)

