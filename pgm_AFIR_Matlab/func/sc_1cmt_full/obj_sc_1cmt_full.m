function[model] = obj_sc_1cmt_full()
%model object

model.name  = get_tmdd_name(mfilename);

model.cmtname =   {'Ad','Ac','T',          'DT'};
model.init  = @(p)([0    0    p.ksyn/p.keT  0]);
model.dcmt  = 1; 					   %compartment for dosing

model.pin   = {'ka','F','V','keD','ksyn','keT','keDT','koff','kon'}; %input params
model.pode  = {'ka','F','V','keD','ksyn','keT','keDT','koff','kon'}; %params used in ode
model.repar = @(p)(p);                          %reparameterizing function to go from pin to pode

model.ode = @ode_sc_1cmt_full;   %differential equation
model.out = @out_sc_1cmt_full;   %convert ode output into informative structure

model.plot= @plot_iv_2cmtc_qe;  %default plotting code for model output
model.check=@check_sc_1cmt_full_Hayashi07; %check to make sure model is implemented properly

model.fun.AFI   = @(p)      (p.koff/p.kon)*(p.keT/p.keDT)*(p.CL*p.tau/(p.F*p.dose));
model.fun.TFI   = @(p,q) abs(p.koff/p.kon)*(p.keT/p.keDT)/...
                        (p.F*p.dose*(q.C*exp(-q.ke*p.tau)/(1-exp(-q.ke*p.tau))) - ...
                         p.F*p.dose*(q.C*exp(-q.ka*p.tau)/(1-exp(-q.ka*p.tau))));
model.fun.Macro = @micro2macro_oral1_1cmt;