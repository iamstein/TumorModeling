function[model] = obj_iv_2cmtc_full_AFI_CLV()
%model object

model       =  obj_iv_2cmtc_full;
model.name  = get_tmdd_name(mfilename);
model.pin   = {'Vc','CL','Q','Vp','ksyn','keT','keDT','AFI','koff'}; %input params
model.repar =  @(p)(repar_AFI2kon(repar_CLV2kxx(p)));