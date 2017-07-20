function[model] = obj_iv_2cmtc_full_AFI_CLV_T0()
%model object

model       =  obj_iv_2cmtc_full;
model.name  = get_tmdd_name(mfilename);
model.pin   = {'Vc','CL','Q','Vp','T0','keT','keDT','AFI','koff'}; %input params
model.repar =  @(p)(repar_AFI2kon(repar_T02ksyn(repar_CLV2kxx(p))));