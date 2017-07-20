function[model] = obj_iv_2cmtc_full()
%model object

model.name  = get_tmdd_name(mfilename);

model.cmtname =   {'Ac','Ap','T',           'DT'};                    
model.init  = @(p)([0    0    p.ksyn/p.keT   0]);
model.dcmt  = 1; 				 				%compartment for dosing

model.pin   = {'Vc','keD','k12','k21','ksyn','keT','keDT','koff','kon'}; %input params
model.pode  = {'Vc','keD','k12','k21','ksyn','keT','keDT','koff','kon'}; %params used in ode
model.repar = @(p)(p);                          %reparameterizing function to go from pin to pode

model.ode   = @ode_iv_2cmtc_full; 				%differential equation
model.out   = @out_iv_2cmtc_full; 				%convert ode output into informative structure
model.plot  = @plot_iv_2cmtc_qe;  				%default plotting code for model output
model.check = @check_iv_2cmtc_full_Elmeliegy15; %check to make sure model is implemented properly

model.fun.AFI   = @(p)   (p.koff/p.kon)*(p.keT/p.keDT)*(p.CL*p.tau/p.dose);
model.fun.TFI   = @(p,q) (p.koff/p.kon)*(p.keT/p.keDT)/...
                      (p.dose*(q.C1*exp(-q.lam1*p.tau)/(1-exp(-q.lam1*p.tau))) + ...
                       p.dose*(q.C2*exp(-q.lam2*p.tau)/(1-exp(-q.lam2*p.tau)))); 
model.fun.Macro = @micro2macro_bolus_2cmt;                   