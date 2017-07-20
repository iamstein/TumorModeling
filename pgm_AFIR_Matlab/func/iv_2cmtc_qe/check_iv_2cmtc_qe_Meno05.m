function[] = check_iv_2cmtc_qe_Meno05(model_obj)

if nargin==0
    model = obj_iv_2cmtc_qe;
else
    model = model_obj;
end

%% load in the data 
    [data Dose] = load_Charles99(0);
            
%% parameters    

%specify model properties
    p      = get_infliximab_Meno05_param; 
    t      = .01:.01:56; %time in days
    
%recreate figure from Wang14
figname('Meno05_Fig8');
hold on

opt.PlotVars = {'Ttot'};
col = col_redblue(length(data));
h = [];

rr = 1; cc = 2;
for i=1:length(data)
    opt.Data  = data(i);

% recreate Figure 3 and 4 from Wang14   
    out = solve_tmdd2(model,t,Dose(i),p);
    
    DisplayName = [num2str(Dose(i).dmpk) ' mpk'];
    
    subplot(rr,cc,1); hold on
        h1(i) = plot(out.t,out.Ctot,'-','Color',col(i,:),'DisplayName',DisplayName);
        title('Infliximab');
        xlabel('Time (Days)');
        ylabel('Conc (nM)');
        set(gca,'YScale','log')
    
    subplot(rr,cc,2); hold on 
        x = data(i).Ttot;
        h2(i) = plot(x.t,x.y,'o-','Color',col(i,:),'DisplayName',DisplayName);
        plot(out.t,out.Ttot,'-','Color',col(i,:));
        title('TNFa');
        xlabel('Time (Days)');
        ylabel('Conc (nM)');   
        
end
subplot(rr,cc,1)
    nolegend(h1,'Position','left','XLimTweakPct',0.4)
    1;
subplot(rr,cc,2)
    nolegend(h2)

    