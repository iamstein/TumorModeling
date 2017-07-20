function[AllSolutions AllParameters] = explore_tmdd_Nparamx(PAR,outvars,logrange,model,t,dose,p0,varargin)
%explore_tmdd_Nparamx(p_names,outvars,model,t,dose,p0)
    
    if isstruct(PAR) || isa(PAR,'dataset')
        p_names = PAR.p_names;
    else
        p_names = PAR;
    end
   
    o.p_units   = cellpop(1,length(p_names),'');
    o.p_deriv   = cellpop(1,length(p_names),'');
    o.titlesc   = ones(1,length(p_names));    
    if isstruct(PAR) || isa(PAR,'dataset')
        f = {'p_units','p_deriv','titlesc'};
        for i=1:length(f)        
            s = f{i};
            if isfield(PAR,s) || dataset_iscol(PAR,s)
                o.(s) = PAR.(s);
            end
        end
    end
    
    np      = length(p_names);
    noutput = length(outvars);

%default options
    o.Axis = [];
    o.SubplotFun = 'subplot'; %could be subplot or subax
    o.SubAxis   = [];
                
    o.XLabel = 'Time';
    if length(outvars)==1
        o.YLabel = '';
    else
        o.Ylabel = cellpop(noutput,1,'');
    end        
        
    
    o.ArrowFlag = 1;
    o.XArrowFrac= 0.8;
    o.YArrowFrac= 0.25;
    o.LineThick = 0;
    o.Ccrit     = 0;    
    o.CcritFun  = @(p) p.Vm/(p.ke*p.Vc);
    
    o.Ccrit2    = 0; %for transition to delta phase
    o.YMin      = NaN;
    
    o.OUT       = [];
    o.P         = [];
    
%replace the default values with anything set in the options structure
    o = parse_options(o,varargin{:});
   
    rr = noutput;
    cc = np;
    
    AllSolutions = cell(1,np);
    AllParameters= cell(1,np);
    for ip=1:np
        p_name  = p_names{ip};
        
        if isnumeric(logrange)
            p_range = logrange*p0.(p_name);
        elseif isstruct(logrange)
            p_range = logrange.(p_name).*p0.(p_name);
        end

        for io=1:noutput
            outvar = outvars{io};

            switch o.SubplotFun
                case 'subplot'
                    subplot2(noutput,np,io,ip)
                case 'subax'
                    o.SubAxis.PaddingTop = 0;
                    subax2(noutput,np,io,ip,o.SubAxis)
            end            
            
            if o.LineThick==1
                opt.LineThick = logrange.(p_name)==1;
            end    
            opt.Axis = o.Axis{io};
            if io==1 %first plot in row
                opt.OUT        = [];
                if ~isempty(o.OUT)
                    opt.OUT        = o.OUT{ip};
                end
                opt.p_units    = o.p_units{ip};
                opt.p_deriv    = o.p_deriv{ip};
                opt.titlesc    = o.titlesc(ip);
                opt.TitleParam = o.TitleParam; 
                opt.Ccrit      = o.Ccrit;
                opt.Ccrit2     = o.Ccrit2;
                opt.CcritFun   = o.CcritFun;   
                opt.YMin       = o.YMin;
                opt.XArrowFrac = o.XArrowFrac;
                opt.YArrowFrac = o.YArrowFrac;
                
                [~, OUT, P] = explore_tmdd_1paramx_SMALL(model,t,dose,p0,p_name,p_range,outvar,opt);                
                1;
            else
                opt.TitleParam = 0;
                opt.Ccrit      = 0;
                opt.YMin       = o.YMin;
                opt.OUT        = OUT; %solution from first time through (so don't need to solve again
                explore_tmdd_1paramx_SMALL(model,t,dose,p0,p_name,p_range,outvar,opt);
            end
            
            if noutput==1
                Axis   = o.Axis;
                YLabel = o.YLabel;
            else
                Axis = o.Axis{io};
                YLabel = o.YLabel{io};
            end
            
            
            if length(Axis)==1
                set(gca,Axis)
            else
                set(gca,Axis{io})
            end
            xlabel(o.XLabel);
            if isempty(YLabel)
                ylabel(outvar);
            else
                ylabel(YLabel);
            end
            
            switch o.SubplotFun
                case 'subax'
                    iplot = np*(io-1)+ip;
                    prettyaxes(rr,cc,iplot)
                    if ~isfield(Axis,'YLim')
                        error('Axis.YLim should be specified as its not shown on every single graph.  Or, you can use subplot insetad of subax');
                    end
                    if io>=2
                        title('');
                    end
            end
            drawnow
        end
        AllSolutions{ip} = OUT;
        AllParameters{ip} = P;
        
        %add Ccrit label
        %{
        if o.Ccrit==1 && ip==np
            s = func2str(opt.CcritFun);
            s = strrep(s,'p.','');
            s = s(5:end);
            s = ['Ccrit = ' s]; %#ok<*AGROW>
            set(gca,'YAxisLocation','right')
            ylabel(s)
        end
        %}
    end