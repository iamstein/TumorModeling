function[data dose] = load_Meno05_omalizumab_fig5a(PlotFlag)

if nargin==0
    PlotFlag=1;
end

%% load in the data     
    mfn     = mfilename('fullpath');
    fdir    = mfn(1:(find(mfn=='/',1,'last')));
    str     = {'PK'  ,'totIgE','freeIgE'};
    strout  = {'Dtot','Ttot'  ,'T'};
    for i=1:length(str)
        fname = [fdir 'Fig5a_' str{i} '.csv'];
        din   = csvread(fname);

        s     = strout{i};                
                
        data.(s).t  = din(:,1);
        data.(s).tu = 'days';
        data.(s).y  = din(:,2);
        data.(s).yu = 'nM';
    end
    
    if PlotFlag==1
        figname('Meno05_Fig5a');
        
        opt.Data = data;
        opt.XLabel = 'Time (days)';
        opt.YLabel = 'Conc (nM)';
        plot_iv_2cmtc_qe([],opt);
    end
      
    d1              = 90*1e6/150e3; %mg*(1e6ng/mg)/150 = nmol
    dose.t          = 0;
    dose.d          = d1;  
    dose.dstr = '90 mg';


    
    filename = [fdir 'Meno05_Data_Omalizumab_fig5a.mat'];
    save(filename,'data','dose');
    %export(ds,'File',filename,'Delimiter',',');
    