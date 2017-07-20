function[model] = obj_iv_2cmtc_full_CLV()
%model object

model       = obj_iv_2cmtc_full;
model.name  = get_tmdd_name(mfilename);
model.pin   = {'Vc','CL','Q','Vp','ksyn','keT','keDT','koff','kon'}; %input params
model.repar = @repar_CLV2kxx;
