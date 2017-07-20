function[info] = info_iv_2cmt_qe_Cen_Wang14()
%information about the model

info.cmt  = (1:3)';
info.name = {'AmtDtot';'AmtDptot';'Ttot'};
info.units= {'amt'    ;'amt'     ;'conc'};
info.ksyn = {''       ;''        ;'conc/t'};
info.longname = {'Total Amount Drug in Central Compartment';...
    'Total Amount of Drug in Periph. Compartment';...
    'Total Concentration of Target in Central Compartment'};

info = dataset(info);