function[model] = obj_sc_1cmt_full_AFI_CLV_Tss_T0_koff()
%model object

model       = obj_sc_1cmt_full;
model.name  = get_tmdd_name(mfilename);
model.pin   = {'ka','F','V','CL','Tss','T0','keDT','koff','AFI'}; %input params
model.repar =  @(p)(repar_AFI2kon(repar_Txx2kxx(repar_CLV2keD1(p))));
