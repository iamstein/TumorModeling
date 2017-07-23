library(RxODE)
library(dplyr)

# This model is for Model F

ivsc_3cmtct_shedct = function() {
  model           = list()
  model$name      = as.character(sys.calls()[[sys.nframe()]])
  
  #COMPARTMENTS AND INITIAL CONDITIONS
  model$cmtshort  = c('AmtD0','D1','D2','D3','S1','S3','M1', 'DM1','M3','DS1','DS3','DM3')
  model$init      = function(p){
             init = c(AmtD0=0,D1=0,D2=0,D3=0,S1=0,S3=0,M1=0,M3=0,DS1=0,DS3=0,DM1=0,DM3=0)
             p    = p %>% t() %>% as.data.frame()
                
                
             ksyn = with(p,c(ksynS1,ksynS3,ksynM1,ksynM3))
             K    = with(p,matrix(c(keS1+k13S         , -k13S*VS1/VS3       ,  -kshedM1,         0,
                                        -k31S*VS3/VS1 ,  keS3 + k31S        ,  0       ,  -kshedM3,
                                         0            ,  0                  , kshedM1+keM1+k13, -k31M*VD3/VD1,
                                         0            ,  0                  , -K13M*VD3/VD1,   kshedM3+keM3+k31M), 
                                nrow = 4, byrow=TRUE))
             x    = solve(K,ksyn)
                
             init["S1"] = unlist(x[1])
             init["S3"] = unlist(x[2])
             init["M1"] = unlist(x[3])
             init["M3"] = unlist(x[4])
             return(init)
  }
  
  #PARAMEETRS IN MODEL
  model$pin       = c('F','ka','VD1','VD2','VD3','VS1','VS3','VDS1','VDS3',
                      'k12D','k21D','k13D','k31D','k13S','k31S','k13DS','k31DS','k13M','k31M',
                      'ksynS1','ksynS3','ksynM3','keD1','keD3','keS1','keS3','keDS1','keDS3','keM3','keDM3',
                      'kon1','koff1','kon3','koff3','kshedM1', 'kshedDM1',
                      'kshedM3','kshedDM3'); #input parameters
  model$pode      = model$pin
  model$repar     = function(p){return(p)} #reparameterization function

                   #INPUT/SYNTHESIS/SHED   DISTRIBUTION (CENTRAL/TUMOR)      BINDING       DISTRIBUTION (CENTRAL/PERIPH)
  model$rxode.str = '
     D1           = AmtD1/VD1;
     d/dt(AmtD0)  =  -ka *AmtD0;
     d/dt(AmtD1)  =(F*ka *AmtD0/VD1 - k13D *D1 + k31D *VD3/VD1*D3 - keD1 *D1  - kon1*D1*S1 + koff1*DS1 - k12D*D1 + k21D*VD2/VD1*D2)*VD1;
     d/dt(D1)     = F*ka *AmtD0/VD1 - k13D *D1 + k31D *VD3/VD1*D3 - keD1 *D1  - kon1*D1*S1 + koff1*DS1 - k12D*D1 + k21D*VD2/VD1*D2; 
     d/dt(D2)     = k12D*VD1/VD2*D1 - k21D*D2;
     d/dt(D3)     = k13D *VD1/VD3*D1     - k31D*D3  - keD3 *D3  - kon3*D3*(S3+M3) + koff3*(DS3+DM3);
     d/dt(S1)     = ksynS1 - k13S *S1 +     k31S*VS3/VS1*S3 + kshedM1*M1  - keS1 *S1  - kon1*D1*S1      + koff1*DS1;
     d/dt(S3)     = ksynS3 +kshedM3*M3   + k13S *VS1/VS3*S1     - k31S*S3  - keS3 *S3  - kon3*D3*S3      + koff3*DS3;
     d/dt(M3)     = ksynM3 -kshedM3*M3 - keM3 *M3 -k31M*M3 + k13M*VD1/VD3*M1- kon3*D3*M3 + koff3*DM3; 
     d/dt(M1)     = ksynM1 - kshedM1*M1 - keM1*M1 - k13M*M1 + k31M*VD3/VD1*M3 - kon1*D1*M1 + koff1*DM1;
     d/dt(DS1)    = - k13DS*DS1 + k31DS*VDS3/VDS1*DS3 - keDS1*DS1 + kon1*D1*S1 - koff1*DS1 +kshedDM1*DM1;
     d/dt(DS3)    =         kshedDM3*DM3 + k13DS*VDS1/VDS3*DS1 - k31DS*DS3 - keDS3*DS3 + kon3*D3*S3      - koff3*DS3;
     d/dt(DM3)    =        -kshedDM3*DM3                                   - keDM3*DM3 + kon3*D3*M3      - koff3*DM3;
     d/dt(DM1)    = -keDM1*DM1 - kshedDM1*DM1 + kon1*D1*M1 - koff1*DM1;
  '

# Q: In d/dt(M1) and d/dt(M3) can we use VD3/VD1? Do We need to introduce VM3/VM1?
# R: k13DS and k31DS in this file are k13DT and k31DT in the Model F part of Andy's slide.




  model$rxode     = RxODE(model = model$rxode.str, modName = model$name)
  
  model$rxout     = function(result,p)    {
    result        = as.data.frame(result)
    result = mutate(result,
                    Dtot1 = D1+DS1,
                    Stot1 = S1+DS1,
                    Dtot3 = D3+DS3,
                    Stot3 = S3+DS3,
                    Mtot3 = M3+DM3)
  }
  
  return(model)
}

