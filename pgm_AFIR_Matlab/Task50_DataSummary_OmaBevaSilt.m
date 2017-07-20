%% header
    clear
    close all;
    startup; 
    
    mfn = mfilename;
    ii  = find(mfn=='_',1,'first');
    mfnpref = mfn(1:ii);

%% setup omalizumab
    s = 'omalizumab';
    [data.(s) dose.(s)] = load_Meno05_omalizumab_fig5a(0);    
    p.(s)               = get_omalizumab_meno05_fig5a_param;    
    model.(s)           = obj_sc_1cmt_full;    
    t.(s)               = 0:.1:63;
    
%% set up bevacizumab
    s = 'bevacizumab';
    [data.(s) dose.(s)] = load_Elmeliegy15_Fig2(0);
    p.(s)               = get_bevacizumab_elmeliegy15_steintweak_param;
    model.(s)           = obj_iv_2cmtc_qe;    
    t.(s)               = 0:.1:90;

%% set up siltuximab
    s = 'siltuximab';
    [data.(s) dose.(s)] = load_Puchalski10_Siltuximab(0,2);            
    p.(s)               = get_siltuximab_mayer15_steintweak_param;
    model.(s)           = obj_iv_2cmtc_qe;
    t.(s)               = 0:.1:90;

%% do simulations and plot data
    clf
    drugs = fieldnames(data);

    opt.Tot.LineStyle = '-';
    opt.Tot.LineWidth = 1;
    opt.PlotVars = {'Dtot','Ttot','T'};
    opt.Axis.XLim  = [-7 120];
    opt.Axis.XTick = 0:21:150;
    opt.Axis.XTickLabel = opt.Axis.XTick/7;
    opt.XLabel   = 'Week';

    %opt.Title    = 'Bevacizumab-VEGF';
    opt.YLabel   = 'Conc (nM)';
    opt.Axis.YLim= [1e-6 2e3];
    opt.Axis.YTick = 10.^(-5:1:3);
    opt.Axis.YTickLabel = prettynum(opt.Axis.YTick);
    opt.Axis.YScale='log';
    
    %nolegend positions
    YPos.omalizumab = [1 3 12];
    YPos.bevacizumab= [.0015 .005 20];
    YPos.siltuximab = [.7e-5 .2 280];
    
    rr = 1; cc = length(drugs);
    for i=1:length(drugs)
        subax(rr,cc,i,'MarginBottom',.2,'MarginTop',.1); hold on       
        s = drugs{i};

        out    = solve_tmdd2(model.(s),t.(s),dose.(s),p.(s));

        opt.Data = data.(s);    
        opt.Title = {[ ('a'+i-1) ') ' s],dose.(s).dstr};
        opt.NoLegend.YPos = YPos.(s);
        
        plot_iv_2cmtc_qe(out,opt);  
        prettyaxes(rr,cc,i);
    end

    fig2file2([mfnpref 'OmaBevaSilt_DataSummary'],[5.5 2.8])

    