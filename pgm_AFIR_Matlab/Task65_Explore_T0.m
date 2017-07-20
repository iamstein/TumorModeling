%% header
    clearvars -except OUT PP
    close all;
    startup;
    
    mfn = mfilename;
    ii  = find(mfn=='_',1,'first');
    mfnpref = mfn(1:ii);

%% run type
    REUSE   = 0; %to reuse previous run
    runtype = {'lum'}; %basic and lumped
    drugs   = {'si'}; %for omalizumab, bevacizumab, siltuximab
    
    %runtype = {'lum'};
    %drugs   = {'si'};
    
%% load in the parameters
    p.om = get_omalizumab_meno05_fig5a_param;
    p.be = get_bevacizumab_elmeliegy15_steintweak_param;
    p.si = get_siltuximab_mayer15_steintweak_param;            
    
%% set model    
    model.om.bas = obj_sc_1cmt_full_CLV;
    model.be.bas = obj_iv_2cmtc_full_CLV;
    model.si.bas = obj_iv_2cmtc_full_CLV;

    model.om.lum = obj_sc_1cmt_full_AFI_CLV_Tss_T0_koff;
    model.be.lum = obj_iv_2cmtc_full_AFI_CLV_Tss_T0_koff;
    model.si.lum = obj_iv_2cmtc_full_AFI_CLV_Tss_T0_koff;
    
    %create a string for the model name
        for idrug=1:length(drugs)
        for itype=1:length(runtype)
            sdrug = drugs{idrug};
            stype = runtype{itype};
            sm    = char(model.(sdrug).(stype).ode);
            model.(sdrug).(stype).name = sm(5:end);
        end
        end
    
    outvars       = {'Ctot'     ,'Ttot'     ,'Tpchg'};%,'DT'     ,'T0_T'     ,'Ttot_T0'     ,'Tpchg'};
    %opt.YLabel    = {'Dtot (nM)','Ttot (nM)','T/T0' };%,'DT (nM)','T0-T (nM)','Ttot-T0 (nM)','T/T0'};
    opt.YLabel = {' ',' ',' '};
    
%% set dosing
    p.si.tau = 21;
    p.be.tau = 14;
    p.om.tau = 28;
    
    dose1mg.om = 150;    %mg
    dose1mg.be = 5*70;   %mg
    dose1mg.si = 3*70; %mg
    
    LastMonth  = 14;    
    tlastd.om  = LastMonth*28;%d
    tlastd.be  = LastMonth*28;%d
    tlastd.si  = LastMonth*28;%d
    
    dt      = .5;
    sout    = '$\\Ttotavg/\\Ttotss = \\{';
    for i=1:length(drugs)
        s = drugs{i}; 
        ps= p.(s);

        t.(s) = [-21 -2 -1 10.^(-9:-1) .2:dt:tlastd.(s)]; %time in days     
    
        ps.dose     = dose1mg.(s)*1e6/150e3; %dose in nmol
        dose.(s).t  = 0;%0:ps.tau:tlastd.(s);
        dose.(s).d  = repmat(ps.dose,1,length(dose.(s).t));
        
        ps.AFI = ps.Kd*(ps.Tss/ps.T0)*(ps.CL*ps.tau/ps.dose);
        
        p.(s)  = ps;
        
        if ~isfield(ps,'F')
            F = 1;
        else
            F = ps.F;
        end
        Davg     = ps.dose*F/(ps.CL*ps.tau);        
        Tfrac    = ps.T0/ps.Tss;
        Ttotfrac = Tfrac*(1-(1-Tfrac)/(ps.Kd/Davg+1))^-1;
        fprintf('%s - Ttotavg/Ttotss = %1.2f\n',s,Ttotfrac);
        
        sout = [sout sprintf('%1.2f,',Ttotfrac)];
    end
    sout(end) = [];
    sout  = [sout '\\}$'];
    disp(sout);

    
%% BASIC - set parameters to explore
    s          = 'bas';
    
    u          = '\uparrow';
    d          = '\downarrow';
    AFI        = 'AFI';
    Kd         = 'Kd';
    kon        = 'kon';
    ksyn       = 'ksyn';
    koff       = 'koff';
    
    xx         = [];
    xx.p_names = {'dose'       ,'tau'  ,'CL'   ,'kon'     ,'koff' ,'ksyn' ,'keT'  ,'keDT' }';
    xx.p_units = {'mg/kg'      ,'d'    ,'L/d'  ,'1/(nM*d)','1/d'  ,'nM/d' ,'1/d'  ,'nM'   }';
    xx.titlesc = [150e3/1e6/70  1       1       1           1       1       1       1    ]';
    xx.p_deriv = {[AFI d]     ,[AFI u],[AFI u],[AFI u]    ,[AFI d],'-'    ,[AFI u],[AFI d]}';
    xx.order   = [ 1            2       3       4           5       8       6       7     ]';
    xx         = dataset(xx);
    xx         = sortrows(xx,'order');
    PAR.(s)    = xx;
    
    
    n             = 5;
    fold2         = 2.^(linspace(-1,1,n));
    fold3         = 3.^(linspace(-1,1,n));
    fold4         = 4.^(linspace(-1,1,n));
    fold10        = logspace(-1,1,n);
    foldx         = fold10;
    logrange.dose = foldx;
    logrange.tau  = fold4;
    logrange.CL   = fold4;
    logrange.kon  = foldx;
    logrange.koff = foldx;           
    logrange.ksyn = foldx;
    logrange.keT  = foldx;
    logrange.keDT = foldx;
    LOGRANGE.(s)  = logrange;

%% LUMPED - setparameters to explore    
    s = 'lum';

    xx         = [];
    xx.p_names = {'dose'      ,'tau'  ,'CL'   ,'AFI'  ,'T0'   ,'Tss'   ,'keDT'  ,'koff'}';
    xx.p_units = {'mg/kg'     ,'d'    ,'L/d'  ,'-'    ,'nM'   ,'nM'    ,'1/d'   ,'1/d' }';
    xx.titlesc = [150e3/1e6/70  1       1       1       1       1        1        1    ]';
    xx.p_deriv = {[kon u]     ,[kon d],[kon d],[kon d],[kon u],[ksyn u],[ksyn u],[kon u]}';
    xx.order   = [ 1            2       3       4       7       5        6        8    ]';
    xx         = dataset(xx);
    xx         = sortrows(xx,'order');
    PAR.(s)    = xx;
    
    PAR.(s)    = PAR.(s)(strcmp(PAR.(s).p_names,'T0'),:);
        
    logrange      = [];
    logrange.T0   = foldx;
    LOGRANGE.(s)  = logrange;
    
    
%% set up axes
    %omalizumab    
        Axis.XLim       = [-14 31];
        Axis.XTick      = 0:(28*3):1000;
        Axis.XTickLabel = Axis.XTick/28;
        
        Axis.YLim       = [.03 1e3];
        Axis.YTick      = 10.^(-10:10);
        Axis.YTickLabel = prettynum(Axis.YTick);
        AX{1}           = Axis;

        Axis.YLim       = [.5 50];
        AX{2}           = Axis;

        Axis.YScale     = 'log';
        Axis.YLim       = [.001 2];
        AX{3}           = Axis;
    
        AXIS.om = AX;
        
    %bevacizumab
        Axis            = AX{1};
        Axis.XLim       = [-28 12*28];
        Axis.XTick      = 0:(28*3):400;
        Axis.XTickLabel = Axis.XTick/28;                
        %Axis.YLim       = [30 3e4];
        AX{1}           = Axis;
        
        AX{2}.YLim      = [3e-4 4e-1];
        AX{2}.XLim      = Axis.XLim;
        AX{2}.XTick     = Axis.XTick;
        AX{2}.XTickLabel= Axis.XTickLabel;
        AX{3}.YLim      = [.01 2];
    
        AX{3}.XLim      = Axis.XLim;
        AX{3}.XTick     = Axis.XTick;
        AX{3}.XTickLabel= Axis.XTickLabel;
        AX{3}.YLim      = [0.05 1.2];
        AXIS.be = AX;
        
    %siltuximab
        AX{2}.YLim      = [3e-6 4e-1];
        AX{3}.YLim      = [.003 2];
        
        AXIS.si = AX;
        
    %position for title label
        yy          = 1.4e5;
        Position.om = [-450 yy 0];
        Position.be = [-130 yy 0];
        Position.si = [-130 yy 0];

    
%% set up the other plotting options    
    opt.SubplotFun = 'subax';    
    opt.TitleParam = 1;
    opt.LineThick  = 0;
    opt.XLabel     = 'Month';
    opt.XArrowFrac = .5;
    opt.YArrowFrac = .2;
    opt.SubAxis.MarginLeft = .2;    
    opt.SubAxis.MarginRight = .1;    
    
%% run model
    fprintf('\n');
    for itype=1:length(runtype)
    for idrug=1:length(drugs)
        sdrug = drugs{idrug};
        stype = runtype{itype};

        opt.Axis    = AXIS.(sdrug);
        
        %reuse previous run
            opt.OUT     = [];
            opt.P       = [];
            if REUSE && exist('OUT','var')
            if isfield(OUT,sdrug)
            if isfield(OUT.(sdrug),stype)
                opt.OUT  = OUT.(sdrug).(stype);
                opt.P    =  PP.(sdrug).(stype);
                fprintf('reusing previous run for %s-%s\n',sdrug,stype);
            end
            end
            end
        
        % run and plot
            figname([sdrug '-' stype]);
            [out P] = explore_tmdd_Nparamx(...
                                    PAR.(stype),...
                                    outvars,...
                                    LOGRANGE.(stype),...
                                    model.(sdrug).(stype),...
                                    t.(sdrug),...
                                    dose.(sdrug),...
                                    p.(sdrug),...
                                    opt);   

            OUT.(sdrug).(stype) = out;
             PP.(sdrug).(stype) = P;    
             
         % add box
            %h = add_title_label(Position.(sdrug));
         
    end
    end
    save(mfn,'OUT','P');
    
    h = get(gcf,'Children');
    set(gcf,'CurrentAxes',h(2))
    plot([-100 1000],3e-4*[1 1],'k:')
    
    % save figure
    siz = [2.5 3.7];
    fig2file2([mfnpref 'T0_silt_lumped_sens'],siz);

    
