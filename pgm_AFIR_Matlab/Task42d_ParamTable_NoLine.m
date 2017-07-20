%% explore TMDD model

%% header
    clear
    close all;
    startup
    mfn = mfilename;
    ii  = find(mfn=='_',1,'first');
    mfnpref = mfn(1:ii);
            
        
%% output table
    %pname = {'drug'        ,'target'        ,'ka'  ,'F','Vc' ,'kpc'  ,'kcp'  ,'keD'  ,'keT'  ,'keDT'  };
    %ptex  = {'\\text{Drug}','\\text{Target}','\\ka','F','V_c','\\kcp','\\kpc','\\keD','\\keT','\\keDT'};
    %units = {'units'       ,''              ,'1/d' ,'' ,'L'  ,'1/d'  ,'1/d'  ,'1/d'  ,'1/d'  ,'1/d'};

%% load in parameters    
    p.om = get_omalizumab_meno05_fig5a_param;
    p.om.drug= 'Omalizumab';
    p.om.k12 = NaN;
    p.om.k21 = NaN;
    p.om.Vc  = p.om.V;
    p.om.Vp  = NaN;
    p.om.Q   = NaN;
    
    p.be = get_bevacizumab_elmeliegy15_steintweak_param;
    p.be.drug= 'Bevacizumab';
    p.be.F   = NaN;
    p.be.ka  = NaN;
    
    p.si = get_siltuximab_mayer15_steintweak_param;
    p.si.drug= 'Siltuximab';
    p.si.F   = NaN;
    p.si.ka  = NaN;
    
    
%% load in table format file    
    df = dataset('File','Task42_TableFormat.csv','Delimiter',',');
    
%% output table    
    shead = {};
    pstr = fieldnames(p);
    shead{end+1} = '\\begin{tabular}{|lclcccc|}';
    shead{end+1} = '\\hline';
    shead{end+1} = ['Type & Param. & Description & ' p.om.drug ' & ' p.be.drug ' & ' p.si.drug ' & Units \\\\']; 
    shead{end+1} = '\\hline';
    sbody = {};
    for ip=1:length(df)
        s = [df.ptype{ip} ' '];
        s = [s '& $' df.tex{ip} '$'];  
        s = [s ' & ' df.desc{ip}];

        pnamei = df.pname{ip};
        for id=1:length(pstr)
            if isempty(pnamei)
                s = [s ' & '];
            else
                num = p.(pstr{id}).(pnamei);
                if isnumeric(num) && ~isnan(num)                
                    s = [s ' & ' prettynum(num,2)]; %#ok<*AGROW>
                elseif isnumeric(num) && isnan(num)
                    s = [s ' & - '];
                else
                    s = [s ' & ' num];
                end
            end
        end
        s = [s ' & ' df.units{ip} '\\\\'];
        sbody{end+1} = s; %#ok<*SAGROW>
    end
    sfoot{1} = '\\end{tabular}';
    
    sout = [shead sbody sfoot];

    fname = ['./' mfnpref 'ParamTable.tex'];
    
    fid = fopen(fname,'w');
    for i=1:length(sout)
        fprintf(fid,sout{i});
        fprintf(fid,'\n');
    end
    fclose(fid);
    
%create a short table for main manuscript
    
    shead = {};
    shead{end+1} = '\\begin{tabular}{|l|c|l|}';
    shead{end+1} = '\\hline';
    shead{end+1} = 'Type & Param. & Description \\\\'; 
    shead{end+1} = '\\hline';
    
    sbodythin = {};
    for i=2:length(sbody)
        s = sbody{i};
        s = strsplit(s,'&');
        s = s([1 2 3]);
        s{end} = [s{end} '\\\\ \\hline'];
        s = strjoin(s,' & ');
        sbodythin{i-1} = s;
    end
    
    sout = [shead sbodythin sfoot];
    fname = ['./' mfnpref 'ParamTable_Thin.tex'];

    fid = fopen(fname,'w');
    for i=1:length(sout)
        fprintf(fid,sout{i});
        fprintf(fid,'\n');
    end
    fclose(fid);