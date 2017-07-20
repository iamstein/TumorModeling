%function[] = check_sc_1cmt_full_Hayashi07

model = obj_sc_1cmt_qe;
model2= obj_sc_1cmt_full;
           
%% parameters    
    p      = get_omalizumab_hayashi07_param; 
    p.koff = 10; %1/d
    p.kon  = p.koff/p.Kd;
    
    t      = 0:.1:80; %time in days
    dose.t = 0;
    dose.d = 0;
    
%% recreate figure 
figname('Hayashi07_Fig3');
rr = 3; cc = 4;
dose_range = [75 150 300 375]; %mg

AxisX.XLim  = [-10 90];
AxisX.XTick = 0:20:80; 

AxisY.Dtot.YLim = [0 60];
AxisY.Ttot.YLim = [0 2500];
AxisY.T.YLim    = [1 300];
AxisY.T.YScale  = 'log';
AxisY.T.YTick   = [5 10 50 100];

PlotVar = {'Dtot','Ttot','T'};
YLabel  = {'Dtot (ug/ml)','Ttot (ng/ml)',{'T (ng/ml)','TIMES 3 (?)'}};

for idose=1:length(dose_range)
    %solve model
        dose.d = dose_range(idose)*1e6/150e3 ;%(1mg)*(1e6ng/mg)*(MW mol/g) = nmol
        out = solve_tmdd2(model,t,dose,p);
        out2= solve_tmdd2(model2,t,dose,p);
                    
    %convert back from nM to ug/ml for omalizumab
    %convert back from nM to ng/ml for IgE and complex 
        out.Dtot = out.Dtot /1000*150;
        out2.Dtot= out2.Dtot/1000*150;
        out.Ttot = out.Ttot     *(190+150);
        out2.Ttot= out2.Ttot    *(190+150);        
        out.T    = out.T        *190 * 3; %<------ FUDGE FACTOR
        out2.T   = out2.T       *190 * 3; %<------ FUDGE FACTOR
        fprintf('Needed a "fudge factor" of 3x on T - not sure why\n');
        
    % loop through plot variables
    for ivar=1:length(PlotVar)
        subplot2(rr,cc,ivar,idose); 
        s = PlotVar{ivar};
        plot(t,out.(s),'k-','DisplayName','qe','LineWidth',2);
        hold on
        plot(t,out2.(s),'r-','DisplayName','full');
        
        if ivar==1
            title(sprintf('DOSE: %d mg',dose_range(idose)));
        end
        if ivar==3
            xlabel('Day');
        end
        if idose==1
            ylabel(YLabel{ivar});
        end
        set(gca,AxisX);
        set(gca,AxisY.(s));
    end         
end
nolegend


    