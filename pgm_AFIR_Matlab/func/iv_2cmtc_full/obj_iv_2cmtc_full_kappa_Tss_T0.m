function[model] = obj_iv_2cmtc_full_kappa_Tss_T0()
%model object

model       = obj_iv_2cmtc_full;
model.name  = get_tmdd_name(mfilename);
model.pin   = model.pin   = {'Vc','keD','k12','k21','Tss','T0','keDT','kon','kappa'}; %input params
model.repar =  @(p)(repar_kappa2koff((repar_Txx2kxx(p))));