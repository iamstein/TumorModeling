%% headers
    clear
    close all

    scale.mpk2nmol = 70*1e-3/150e3*1e9;
    scale.nmol2mpk = 1/scale.mpk2nmol; %nM->mg/kg

    scale.nM2ngml.mAb  = 150;
    scale.nM2ugml.mAb  = 150/1000;
    scale.ngml2nM.mAb  = 1/150;
    scale.ugml2nM.mAb  = 1000/150;
    
%% load in the pk dataset
    mfn        = mfilename('fullpath');
    datadir    = mfn(1:(find(mfn=='/',1,'last')));    
    din        = dataset('File',[datadir 'Digitized_PKBinding_Beva_Elmeleigy.csv'],'Delimiter',',');

%% expand pk to generalized dataset
    o = [];
    o.InputDataset = din;
    o.REF      = 'elmeleigy15';
    o.DRUG     = 'bevacizumab';
    o.SPECIES  = 'human';
    o.INDNAME  = 'healthy males';

    d          = general_dataset(o);       
    d.ID       = d.DOSEGRP;
    
    dorig      = d; %original dataset before scaling
    
%% output Generalized Dataset        
    filename = [datadir 'Elmeleigy15_Bevacizumab_GenDataset.csv'];
    export(d,'File',filename,'Delimiter',',');    
    
%% convert to nmol and nM
    idose = strcmp(d.TYPENAME,'dose');
    ipk   = strcmp(d.TYPENAME,'pk');
    ipd   = strcmp(d.TYPENAME,'pd');

    d.VALUE(idose) = d.VALUE(idose)*scale.mpk2nmol;
    d.VALUE(ipk)   = d.VALUE(ipk)  *scale.ugml2nM.mAb;
    d.LLIM(ipk)    = d.LLIM(ipk)   *scale.ugml2nM.mAb;
    d.ULIM(ipk)    = d.ULIM(ipk)   *scale.ugml2nM.mAb;    
    
%% create Nonmem Dataset
    %initialize dataset
        nonmem_fields = {'ID','AMTnmol','TINF','ADM','TIMEday','DVnM','CENS','MDV'};
        n = size(d,1);
        for i=1:length(nonmem_fields)
            dn.(nonmem_fields{i}) = zeros(n,1);
        end
        
    %ID
        dn.ID = d.ID;        
        
    %DOSING
        dn.AMTnmol(idose) = d.VALUE(idose);
        if all(~isnan(d.ENDTIME))
            dn.TINF(idose) = d.ENDTIME(idose);
        end
        dn.ADM(idose) = 2;
        dn.ADM(idose & ismember(d.ROUTE,{'oral','subcutaneous','intramuscular'})) = 1;
    
    %TIME and PK   
        LLIM = d.LLIM;
        LLIM(isnan(LLIM)) = -Inf;
    
        dn.TIMEday = d.STARTIME;
        dn.CMT        = ones(n,1);
        dn.CMT(idose) = 1;
        dn.CMT(ipk)   = 2;
        dn.CMT(ipd)   = 3;
        dn.DVnM(~idose) = d.VALUE(~idose);
        dn.MDV(idose) = 1;
        icens         = dn.DVnM<=LLIM & LLIM~=-Inf & ~idose;
        dn.DVnM(icens)= LLIM(icens);
        dn.CENS(icens)= 1;
        
    
        
%% output Generalized Dataset
    dn       = dataset(dn);
    filename = [datadir 'Elmeleigy15_Bevacizumab_NLME.csv'];
    export(dn,'File',filename,'Delimiter',',');

%% plot dataset
    o.field_x     = 'STARTIME';
    o.field_y     = 'VALUE';
    o.field_id    = 'NAME';
    o.field_color = 'NAME';     
    o.Axis.YScale = 'log';    
    
    dplot = dorig(ismember(dorig.TYPENAME,{'pk','pd'}),:);
    plot_ds(dplot,o)
    xlabel(dplot.TIMEUNIT{1});
    ylabel('pk (ug/ml), target (pg/ml)');
    %ylabel(dplot.VALUNIT{2});