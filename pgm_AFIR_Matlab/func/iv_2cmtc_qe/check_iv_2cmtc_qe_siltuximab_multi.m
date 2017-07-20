function[] = check_iv_2cmtc_qe_siltuximab_multi(model_obj)

if nargin==0
    model = obj_iv_2cmtc_qe;
else
    model = model_obj;
end


%% load in the data 
    [data Dose] = load_Puchalski10_Siltuximab(0);
            
%% parameters    
    p      = get_siltuximab_mayer15_steintweak_param;
    
    %p.kappa= 20;
    %p.Tss  = 8e-2;
    %p.T0   = 1e-3;
    %p.keDT = .04;
    
    t      = 0:.1:92; %time in days
    
%% simulate data    
    out = {};
    for i=1:length(data)
        out{i} = solve_tmdd2(model,t,Dose(i),p);
    end
    
%% plot figure
    clf
    rr = 1; cc = 2;
    
    Axis.XLim = [-5 100];
    Axis.XTick = 0:21:100;
    Axis.XTickLabel = Axis.XTick/7;
    opt.XLabel = 'Week';
    
    Axis.YScale = 'log';    
    opt.YLim   = 'Conc (nM)';
    
    opt.Axis = Axis;
    
    for i=1:length(data)
        opt.Data = data(i);
        subplot(rr,cc,i) 
            plot_iv_2cmtc_qe(out{i},opt);
            title(sprintf('%d mpk',Dose(i).dmpk));            
    end
    1;