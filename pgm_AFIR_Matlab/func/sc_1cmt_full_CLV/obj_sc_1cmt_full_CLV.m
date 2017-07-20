function[model] = obj_sc_1cmt_full_CLV()
%model object

model       = obj_sc_1cmt_full;
model.name  = get_tmdd_name(mfilename);
model.pin   = {'ka','F','V','CL','ksyn','keT','keDT','koff','kon'}; %input params
model.repar = @repar_CLV2keD1;
