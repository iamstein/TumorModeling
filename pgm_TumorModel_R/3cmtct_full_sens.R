library(grid)
library(gridExtra)
library(tidyverse) 
setwd("~/GitHub/TumorModeling/pgm_TumorModel_R")
source("ams_initialize_script.R")
dirs = get.dirs(sys.calls(),dirs)

koff3 = seq(1,1000,200)
koff1 = seq(1,1000,200)
kon1  = seq(1,1000,200)
kon3  = seq(1,1000,200)

params <- list(koff3 = koff3, koff1 = koff1, kon1 = kon1, kon3 = kon3)
drugs <- c("Atezo", "Pembro", "Herceptin")

n = length(params)

# Start going through parameter lists.
for (i in 1:4){
    print(i)
    T1Total = NULL
    D1Total = NULL
    T3Total = NULL
    D3Total = NULL
    T3_realtot = NULL
    T3oT30 = NULL
    # Change which paramater I am modifying.
    k = length(params[[i]])
    labels = as.character(params[[i]])
    for (j in 1:k){
        mod      = ivsc_3cmtct_full()
        d        = xlsx::read.xlsx("../data/Bx_DTN_example.xlsx",1)
        # d        = xlsx::read.xlsx("./miandra/Atezolizumab.xlsx",1)
        # Finding the row in params which has the target parameter.
        row <- which(as.character(d$Parameter) == names(params)[i]) 
        # Change the parameter of interest. 
        d[row,3] = params[[i]][j]
        p        = d$Value
        names(p) = d$Parameter
        
        # Store T_30 from data file.
        # T_30 = d[which(as.character(d$Parameter) == "T30"),3]
        T_30 = d[19,3]/d[21,3]
            
        # These functions aren't strictly needed here.  just testing because will be useful later
        p        = p[mod$pin]
        p        = mod$repar(p)
        
        p["keT3"]= p["keT1"]
        
        # Simulate model and put into OUT
        ndose     = 1
        ev       = eventTable(amount.units="nmol", time.units = "days")
        ev$add.sampling(unique(sort(c(seq(-7,13*7,.1),10^(-9:0)))))
        ev$add.dosing(dose=300/ndose*scale.mpk2nmol,nbr.doses=4*ndose,dosing.interval=21/ndose,dosing.to=2)
        
        out         = mod$rxode$solve(p, ev, mod$init(p))
        out         = mod$rxout(out,p)
        out = out %>% mutate(ABCfree = D3/D1,ABCtot  = Dtot3/Dtot1)
        
        # plot data ----
        d = out %>%
            dplyr::select(time,D1,D3,T1,T3,DT1,DT3,ABCtot) %>%
            gather(variable,value,-time) %>%
            mutate(value = ifelse(value>1e-8,value,NA))
        
        # Split data to calculate total values in central and tumor.
        output <- split(d, d$variable)
        n <- nrow(output$ABCtot)
        # Create an data.frame to collect the sums.
        T1 <- data.frame(time = d$time, variable = rep(labels[j], n), value = output$T1$value + output$DT1$value)
        T1Total <- rbind(T1Total, T1)
       
        D1 <- data.frame(time = d$time, variable = rep(labels[j], n), value = output$D1$value + output$DT1$value)
        D1Total <- rbind(D1Total, D1)
        
        T3 <- data.frame(time = d$time, variable = rep(labels[j], n), value = output$T3$value + output$DT3$value)
        T3Total <- rbind(D3Total, T1)
        
        D3 <- data.frame(time = d$time, variable = rep(labels[j], n), value = output$D3$value + output$DT3$value)
        D3Total <- rbind(D3Total, D3) 
        
        T3_real <- data.frame(time = d$time, variable = rep(labels[j], n), value = output$T3$value/T_30)
        T3oT30 <- rbind(T3oT30, T3_real)
    }
    
    if (TRUE){
        g1 = ggplot(T1Total, aes(x = time,y=value,color=variable)) + geom_line(size=1) + scale.y.log10() + ylab("[ T1_tot] (nM)") + 
            theme(legend.position="none",axis.title.x=element_blank(),axis.text.x=element_blank())
        g3 = ggplot(T3Total, aes(x = time,y=value,color=variable)) + geom_line(size=1) + scale.y.log10() + ylab("[ T3_tot] (nM)") + 
            theme(legend.position="none",axis.title.x=element_blank(),axis.text.x=element_blank())
        g2 = ggplot(D1Total, aes(x = time,y=value,color=variable)) + geom_line(size=1) + scale.y.log10() + ylab("[ D1_tot] (nM)") +
            theme(legend.position="none",axis.title.x=element_blank(),axis.text.x=element_blank())
        g4 = ggplot(D3Total, aes(x = time,y=value,color=variable)) + geom_line(size=1) + scale.y.log10() + ylab("[ D3_tot] (nM)") +
            theme(legend.position="none",axis.title.x=element_blank(),axis.text.x=element_blank())
        g5 = ggplot(T3oT30, aes(x = time,y=value,color=variable)) + geom_line(size=1) + scale.y.log10() + ylab("T3/T30") + theme(legend.position="none")
        graph <- grid.arrange(g2, g4, g1, g3, g5, ncol = 1)
        
        # print(graph)
        # filename <- paste0("./miandra/output/graphDTN",names(params)[i],".png")
        filename <- paste0("./miandra/output/graph",drugs[i],names(params)[i],".png")
        ggsave(file = filename,graph)
        
        # Test to see if ratios are different for T3/T30. They are not.
        # l <- c("801", "601")
        # dat <- T3oT30 %>% filter(!(variable %in% l))
        
        # Test to see why there is a change in T3_Tot
        # l <- c("601","401", "201", "1")
        # dat <- T3Total %>% filter(!(variable %in% l))
        # ggplot(dat, aes(x = time,y=value,color=variable)) + geom_line(size=1) + scale.y.log10() + ylab("[ T3_tot] (nM)")
    }
}
