function[data Dose ds] = load_Puchalski10_Siltuximab(PlotFlag,ind)

if nargin==0
    PlotFlag=1;
end
if nargin<2
    ind = 1:2; %data indexes to output
end

MW.D  = 150;
MW.T  =  25;
MW.DT = MW.D + MW.T; 

%% load in the data 
    dose_load = [3 6]; %mpk
    
    mfn        = mfilename('fullpath');
    fdir = mfn(1:(find(mfn=='/',1,'last')));
    data = [];
    ds   = [];
    for i=1:length(dose_load)
        d1 = [];
        dosestr = num2str(dose_load(i));
        
        fname = [fdir 'Fig1b_PK_' num2str(dose_load(i)) 'mpk.csv'];
        dpk   = csvread(fname);

        fname = [fdir 'SuppFig2b_IL6_' num2str(dose_load(i)) 'mpk.csv'];
        dtarg  = csvread(fname);       
                
        Dtot.t      = dpk(:,1); %#ok<*AGROW>
        Dtot.tu     = 'd';
        Dtot.y_ugml = dpk(:,2);
        Dtot.y      = dpk(:,2)*1000/MW.D;
        Dtot.yu     = 'nM';
        npk         = size(dpk,1);
        
        Ttot.t      = dtarg(:,1);
        Ttot.tu     = 'd';
        Ttot.y_pgml = dtarg(:,2);        
        Ttot.y      = dtarg(:,2)/1e3/MW.DT;
        Ttot.yu     = 'nM';
        ntarg       = size(dtarg,1); 
        ntot        = npk+ntarg;
        
        Dose(i).t = [0 21 42 63];
        Dose(i).d = dose_load(i)*70*1e6/150e3*ones(size(Dose(i).t));
        Dose(i).du= 'nmol';
        Dose(i).dmpk = dose_load(i);
        Dose(i).dstr = sprintf('%d mg/kg every 3 weeks',dose_load(i));
        
        data(i).Dtot = Dtot;
        data(i).Ttot = Ttot;
                
        d1.id        = repmat(dose_load(i),ntot,1);
        d1.dose_mpk  = repmat(dose_load(i),ntot,1);
        d1.dose_nmol = d1.dose_mpk*70/1000/150e3*1e9;
        d1.t         = [Dtot.t; Ttot.t];
        d1.tu        = repmat({'day'},ntot,1);
        d1.y         = [Dtot.y; Ttot.y];
        d1.yu        = repmat({'nM'},ntot,1);
        d1.type      = [cellpop(npk,1,'Dtot'); cellpop(ntarg,1,'Ttot')];
        d1.dosetype  = [cellpop(npk,1,[dosestr 'Dtot']); cellpop(ntarg,1,[dosestr 'Ttot'])];
        d1           = dataset(d1);
        
        ds = [ds; d1];
    end
    
%% plot results
if PlotFlag==1
    opt.field_id        = 'dosetype';
    opt.field_x         = 't';
    opt.field_y         = 'y';
    opt.field_color     = 'type';
    opt.field_subplot   = 'id';
    opt.ColorMat        = [0 0 .5; .5 0 0];
    
    Axis.XLim = [-5 100];
    Axis.XTick = 0:7:100;
    Axis.XTickLabel = Axis.XTick/7;
    opt.XLabel = 'Week';
    
    Axis.YScale = 'log';    
    opt.Axis.YLim   = 'Conc (nM)';
    
    opt.Axis = Axis;
    
    plot_ds(ds,opt);
    1;
end

p    = mfilename('fullpath');
fdir = p(1:find(p=='/',1,'last'));
filename = [fdir 'Puchalski10_Pooled.csv'];
export(ds,'File',filename,'Delimiter',',');

data = data(ind);
Dose = Dose(ind);
    
    