function[h OUT P] = explore_tmdd_1paramx(model,t,dose,p0,p_name,p_range,outvar,varargin)
%explore_tmdd_1param(model,t,dose,p0,p_name,p_range,outvar,OPTIONS)

n = length(p_range);
o.OUT = [];

o.FitProps.LineStyle = '-';
o.ColorIndiv         = diverging_map(linspace(0,1,n),[0 0 0.6],[0.6 0 0]);
if mod(n,2)==1
    o.ColorIndiv(n/2+.5,:) = .4*[1 1 1];
end
o.XLabel = 'Time';
o.YLabel = outvar;

o.Ccrit       = 0;
o.YMin        = NaN;

o.Axis        = [];

o.LineThick   = false(1,n);
o.p_units     = '';
o.p_deriv     = '';
o.titlesc     = 1;
o.TitleParam  = 0;
o.ArrowFlag   = 1;
o.XArrowFrac  = .8;
o.YArrowFrac  = .25;
o = parse_options(o,varargin{:}); 

dose0 = dose;
hold on    
for i=1:length(p_range)
    p = p0;
    
    p.(p_name) = p_range(i);
    switch p_name
        case 'dose'
            dose   = dose0;
            dose.d = ones(size(dose.t))*p.dose;
        case 'tau' 
            dose   = dose0;
            dose.t = 0:(p.tau):dose.t(end);
            ind    = mod(1:length(dose.t),length(dose0.t))+1;
            dose.d = dose0.d(ind);
    end

    %solve model
        if isempty(o.OUT)
            if ischar(model) %for solve_tmdd (v1)
                out{i}     = solve_tmdd(model,[],t,dose,p);
            else %for solve_tmdd2 using function
                out{i}     = solve_tmdd2(model,t,dose,p);
            end                    
        else
            out{i} = o.OUT{i};
        end
        y = out{i}.(outvar);      
        
        if ~isnan(o.YMin)
            y(y<o.YMin) = o.YMin;
        end
        
        if isfield(model,'pin')
            p = strip_param(p,model.pin);
        else
            warning('no pin specified, model names may be redundant');
        end
        P{i} = p;
        1;
        
    %plot model
        h(i) = plot(t,y,'Color',[.5 .5 .5]); %#ok<*AGROW>
        if ~isempty(o.ColorIndiv)
            set(h(i),'Color',o.ColorIndiv(i,:));
        end
        setrep(h(i),o.FitProps)
        
        if o.LineThick(i)
            set(h(i),'LineWidth',4);
        end
        
        if o.Ccrit==1
            if ~isstruct(p)
                p = dataset2struct(p);
            end
            
            Ccrit = feval(o.CcritFun,p);               
            ii    = y>1e-10;            
            tcrit = interp1(y(ii),t(ii),Ccrit);
            
            plot(tcrit,Ccrit,'ko','MarkerFaceColor',o.ColorIndiv(i,:));
            1;
        end
        %if o.Ccrit2==1            
        %    1; 
        %end        
        
    %set axes properties
        set(gca,'YScale','log');
        xlabel(o.XLabel);
        ylabel(o.YLabel);
        
        1;
end

set(gca,o.Axis);

if o.ArrowFlag==1
    XLim   = get(gca,'XLim');
    XArrow = (1-o.XArrowFrac)*XLim(1) + o.XArrowFrac*XLim(2);
    YLim   = get(gca,'YLim');
    YArrow = 10.^((1-o.YArrowFrac)*log10(YLim(1)) + o.YArrowFrac*log10(YLim(2))) ;
    
    %metricN = trapz(out{end}.t,(out{end}.(outvar)));
    %metric1 = trapz(  out{1}.t,(  out{1}.(outvar)));
        
    ind = out{1}.t >= XLim(1) & out{1}.t <= XLim(2);    
    tt = out{1}.t(ind);
    y1 = log10(out{1}.(outvar)(ind));
    y2 = log10(out{end}.(outvar)(ind));
    ydiff = trapz(tt,y2-y1)/(tt(end)-tt(1));
    direction_threshold = .08;
    
    ArrowText{1} = '';    
    low  = sprintf('\\color[rgb]{%1.1f,%1.1f,%1.1f}low' ,o.ColorIndiv(1  ,1),o.ColorIndiv(1  ,2),o.ColorIndiv(1  ,3));
    mid  = '\color[rgb]{0,0,0}';
    high = sprintf('\\color[rgb]{%1.1f,%1.1f,%1.1f}high',o.ColorIndiv(end,1),o.ColorIndiv(end,2),o.ColorIndiv(end,3));
    if ydiff > direction_threshold
        ArrowText{1} = high;
        ArrowText{2} = [mid '\uparrow'];
        ArrowText{3} = low;
    elseif ydiff < -direction_threshold
        ArrowText{1} = low;
        ArrowText{2} = [mid '\downarrow'];
        ArrowText{3} = high;
    end
    
    if ~isempty(ArrowText)   
        text(XArrow,YArrow,ArrowText...
            ,'HorizontalAlignment','Center','BackgroundColor','w','Margin',0.001...
        );
    end
    1;
end        

tstr = {p_name};
if o.TitleParam==1
    tstr{1} = [p_name ' ('  o.p_units ')'];
    
    %values of parameter
        s1 = sprintf('\\color[rgb]{%1.1f,%1.1f,%1.1f}%s',...
            o.ColorIndiv(1,1),o.ColorIndiv(1,2),o.ColorIndiv(1,3),prettynum(min(p_range*o.titlesc),2));

        s2 = sprintf('\\color[rgb]{%1.1f,%1.1f,%1.1f}%s',...
            o.ColorIndiv(end,1),o.ColorIndiv(end,2),o.ColorIndiv(end,3),prettynum(max(p_range*o.titlesc),2));

        bk = '\color[rgb]{0,0,0}';
        tstr{2} = [s1 bk '-' s2 bk];
    
    %derived parameter (value that changes so others can stay fixed)
        if ~isempty(o.p_deriv)
            tstr{3} = ['\color[rgb]{.5 .5 .5}' o.p_deriv];
        end    
end
title(tstr)
OUT = out;