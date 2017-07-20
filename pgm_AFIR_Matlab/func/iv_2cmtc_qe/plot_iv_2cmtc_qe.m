function[H] = plot_iv_2cmtc_qe(out,varargin)
%default plotting code

%plot type
    o.PlotType       = 1;  %1=all in one;  3=total conc, total targ, free targ

%data
    o.Data           = []; %contains data to plot with the model
    o.DataProps.Marker = 'o';
    
%variables to plot - if empty use default based on plot type (below)
    o.PlotVars       = [];

%line properties for PlotType==1
    o.Tot.LineStyle = '--';
    o.Tot.LineWidth = 2;

    o.Dtot.Color     = [0 0 .5];
    o.D.Color        = [0 .7 .7];
    o.Ttot.Color     = [.5 0 0];
    o.T.Color        = [1 .1 0];
    o.DT.Color       = [1  0 1];
    
%line properties for PlotType==3 (same)
    if isstruct(out)
        o.Model.Color    = [0 0 1];
        o.Model.LineWidth= 2;    
    elseif iscell(out)
        n                = length(out);
        o.Model.Color    = col_redblue(n);
        o.Model.LineWidth= 2*ones(n,1);
        o.Model = dataset(o.Model);
    end

%axis properties
    o.Axis.YScale    = 'log';   
    o.AxisInd        = {}; %Axis properties for indiv subplot
    o.XLabel         = 'Time';
    o.YLabel         = 'Conc';
    o.Title          = {};
    
%default values
    o = parse_options(o,varargin{:});
    if ~iscell(o.PlotVars) && ~isempty(o.PlotVars)
        o.PlotVars = {o.PlotVars};
    end
    
    data = o.Data;
    
switch o.PlotType
    case {1}  
        hold on
        set(gca,o.Axis)        

        %variabels to plot
            PlotVars = o.PlotVars;
            if isempty(PlotVars)
                PlotVars = {'Dtot','Ttot','D','T','DT'};
            end            

        %plot model
            H = [];
            if ~isempty(out)
                for i=1:length(PlotVars)
                    s = PlotVars{i};
                    h = plot(out.t,out.(s),'DisplayName',s);
                    set(h,o.(s));	            

                    if ~isempty(strfind(s,'tot'))
                        if ~isempty(o.Tot)
                            set(h,o.Tot);
                        end
                    end
                    H(i) = h;
                end
            end
        
        %plot data
            if ~isempty(data)
                for i=1:length(PlotVars)
                    s = PlotVars{i};
                    if isfield(data,s)
                        h = plot(data.(s).t,data.(s).y,'o','Color',o.(s).Color);
                        if isempty(out)
                            set(h,'DisplayName',s);
                            H(i) = h;
                        end
                    end
                end
            end   
            nolegend(H)%,o.NoLegend);
        
        %set axes
            xlabel(o.XLabel);
            ylabel(o.YLabel);
            title(o.Title);

    case {3}
        %variables to plot
            PlotVars = o.PlotVars;
            if isempty(PlotVars)
                PlotVars = {'Dtot','Ttot','Tpchg'};
            end
        
            rr = 1; cc = length(PlotVars);        
        
        %get variables into right from
            if isstruct(out) %single model
                OUT  = {out};
                DATA = data;
                MODEL= struct2dataset(o.Model);
            else
                OUT = out;
                DATA= data;
                MODEL = o.Model;
            end        

        %plot the variables
            for i=1:length(PlotVars)
                s = PlotVars{i};
                subplot(rr,cc,i)
                hold on

                %set axis properties
                    set(gca,o.Axis)
                    if ~isempty(o.AxisInd)
                        set(gca,o.AxisInd{i});
                    end

                for j=1:length(OUT)
                    out = OUT{j};
                    ModelProps = dataset2struct(MODEL(j,:));
                        
                    %plot data             
                        if ~isempty(data)
                            if j<=length(DATA) 
                                data = DATA(j);
                                if isfield(data,s)
                                    h = plot(data.(s).t,data.(s).y,'ko');
                                    set(h,o.DataProps);
                                    if length(MODEL)==length(DATA) && length(MODEL)>1
                                        set(h,'Color',ModelProps.Color)
                                    end
                                end
                            end
                        end

                    %plot model
                        h = plot(out.t,out.(s),'b');
                        H(j,i) = h;                                    
                        set(h,ModelProps);
                        1;
                end                                

                %set labels and titles
                    xlabel(o.XLabel);
                    if iscell(o.YLabel)
                        ylabel(o.YLabel{i})
                    else
                        ylabel(s);
                    end
                    if ~isempty(o.Title)
                        title(o.Title{i})
                    else
                        title(s);
                    end
            end
end