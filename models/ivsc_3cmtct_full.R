library(RxODE)
library(dplyr)

ivsc_3cmtct_full = function() {
  model           = list()
  model$name      = as.character(sys.calls()[[sys.nframe()]])
  
  #COMPARTMENTS AND INITIAL CONDITIONS
  model$cmtshort  = c('AmtD0','D1','D2','D3','T1','T3','DT1','DT3')
  model$init      = function(p){
             init = c(AmtD0=0,D1=0,D2=0,D3=0,T1=0,T3=0,DT1=0,DT3=0)
             p    = p %>% t() %>% as.data.frame()
                
                
             ksyn = with(p,c(ksyn1,ksyn3))
             K    = with(p,matrix(c(keT1+k13T    , -k13T*VT1/VT3,
                                   -k31T*VT3/VT1 ,  keT3+k31T), nrow = 2, byrow=TRUE))
             x    = solve(K,ksyn)
                
             init["T1"] = unlist(x[1])
             init["T3"] = unlist(x[2])
             return(init)
  }
  
  #PARAMEETRS IN MODEL
  model$pin       = c('F','ka','VD1','VD2','VD3','VT1','VT3','VDT1','VDT3',
                      'k12D','k21D','k13D','k31D','k13T','k31T','k13DT','k31DT',
                      'ksyn1','ksyn3','keD1','keD3','keT1','keT3','keDT1','keDT3',
                      'kon1','koff1','kon3','koff3'); #input parameters
  model$pode      = model$pin
  model$repar     = function(p){return(p)} #reparameterization function

  model$rxode.str = '
     D1           = AmtD1/VD1;
     d/dt(AmtD0)  =  -ka *AmtD0;
     d/dt(AmtD1)  =(F*ka *AmtD0/VD1- k13D *D1 +    k31D *VD3/VD1*D3  - keD1 *D1  - kon1*D1*T1 + koff1*DT1 - k12D*D1 + k21D*VD2/VD1*D2)*VD1;
     d/dt(D2)     =                                                                                         k12D*VD1/VD2*D1 - k21D*D2;
     d/dt(D3)     =                  k13D *VD1/VD3*D1     - k31D*D3  - keD3 *D3  - kon3*D3*T3 + koff3*DT3;
     d/dt(T1)     = ksyn1          - k13T *T1 +     k31T*VT3/VT1*T3  - keT1 *T1  - kon1*D1*T1 + koff1*DT1;
     d/dt(T3)     = ksyn3          + k13T *VT1/VT3*T1     - k31T*T3  - keT3 *T3  - kon3*D3*T3 + koff3*DT3;
     d/dt(DT1)    =                - k13DT*DT1 + k31DT*VDT3/VDT1*DT3 - keDT1*DT1 + kon1*D1*T1 - koff1*DT1;
     d/dt(DT3)    =                  k13DT*VDT1/VDT3*DT1 - k31DT*DT3 - keDT3*DT3 + kon3*D3*T3 - koff3*DT3;
  '
  model$rxode     = RxODE(model = model$rxode.str, modName = model$name)
  
  model$rxout     = function(result,p)    {
    result        = as.data.frame(result)
    result = mutate(result,
                    Dtot1 = D1+DT1,
                    Ttot1 = T1+DT1,
                    Dtot3 = D3+DT3,
                    Ttot3 = T3+DT3)
  }
  
  return(model)
}

