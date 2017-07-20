function[data dose] = load_Elmeliegy15_Fig2(PlotFlag)

if nargin==0
    PlotFlag=1;
end

%% load in the data     
    mfn     = mfilename('fullpath');
    fdir    = mfn(1:(find(mfn=='/',1,'last')));
    str     = {'PK'  ,'VEGF'};
    strout  = {'Dtot','Ttot'};
    for i=1:length(str)
        fname = [fdir 'Fig2_' str{i} '.csv'];
        din   = csvread(fname);

        s     = strout{i};                
        switch s
            case{'Dtot'}
                din(:,2) = din(:,2)*1000/150; %convert ug/ml to nM
            case{'Ttot'}
                din(:,2) = din(:,2)/1000/(19.5); %convert pg/ml to nM
        end
        
        
        data.(s).t  = din(:,1);
        data.(s).tu = 'days';
        data.(s).y  = din(:,2);
        data.(s).yu = 'nM';
    end
    
    if PlotFlag==1
        figname('Elmeliegy_Fig2');
        
        opt.Data = data;
        opt.XLabel = 'Time (days)';
        opt.YLabel = 'Conc (nM)';
        plot_iv_2cmtc_qe([],opt);
    end
      
    d1     = 5*70*1e6/150e3; %5mg/kg*70kg*(1e6ng/mg)/150 = nmol
    dose.t = 0;
    dose.d = d1;    
    dose.dstr = '5 mg/kg';
    
    filename = [fdir 'Elmeliegy15_Fig2.mat'];
    save(filename,'data','dose');
    %export(ds,'File',filename,'Delimiter',',');
    