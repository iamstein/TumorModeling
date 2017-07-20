function[model] = obj_iv_2cmtc_qe()
%model object

model.info= @info_iv_2cmtc_qe;  %information about model compartments 
model.init= @init_iv_2cmtc_qe;  %default initial condition setting
model.dcmt= 1; 					   %compartment for dosing
model.ode = @ode_iv_2cmtc_qe;   %differential equation
model.out = @out_iv_2cmtc_qe;   %convert ode output into informative structure
model.plot= @plot_iv_2cmtc_qe;  %default plotting code for model output
model.check=@check_iv_2cmtc_qe_Meno05; %check to make sure model is implemented properly
