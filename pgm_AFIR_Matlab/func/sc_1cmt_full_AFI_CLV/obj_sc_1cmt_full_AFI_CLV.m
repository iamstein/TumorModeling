function[model] = obj_sc_1cmt_full_AFI_CLV()
%model object


model       = obj_sc_1cmt_full;
model.name  = get_tmdd_name(mfilename);
model.pin   = {'ka','F','V','CL','ksyn','keT','keDT','AFI','koff'}; %input params
model.repar =  @(p)(repar_AFI2kon(repar_CLV2keD1(p)));
