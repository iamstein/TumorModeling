function[varargout] = nolegend(varargin)
%labels the lines in a graph (instead of using legend) as suggested by Tufte
%nolegend(labels)
%nolegend(h,labels)
%nolegend(labels,OPTIONS)
%nolegend(h,labels,OPTIONS)
%
%h      = the handles for each line - numeric (if not entered, default is 
%         get(gca,'Children')
%labels = the labels to output - cell of strings.  it must be entered
%OPTIONS= a structure containing data like SpacingPct, XOffsetPct,
%         Position, more info provided below.  This uses the input parser,
%         so options could also be entered as PropertyName, PropertyValue
%
%example code
%{
t = 0:.1:19.6;
rr = 3; cc = 1;
subplot(rr,cc,1)
plot(t,sin(t),t,cos(t))
nolegend({'sin','cos'});

subplot(rr,cc,2)
plot(t,sin(t),t,cos(t))
nolegend({'sin','cos'},'SpacingPct',.20);

subplot(rr,cc,3)
plot(t,sin(t),t,cos(t))
opt = [];
nolegend({'sin','cos'},'Position','left');

for i=1:rr
    subplot(rr,cc,i)
    ylim([-1.5 1.5])
end
%}
 

%Default Parameters    
    %TextProperties = [];
    o.SpacingPct     = .05; %minimum spacing as a fraction of the Y-size of the axis
    o.XOffsetPct     = .02; %X offset as a fraction of the size of the X axis
    o.Position       = 'right'; %could be left or right       
    o.XLimTweak      =   1; %Tweak the XLim a bit so text is more likely to be in plot
    o.XLimTweakPct   =  .1;
    o.FontWeight     = 'normal';
    o.XPos           = [];  %manually specify X position
    o.YPos           = [];  %manually specify the y positions
    o.MiddleOffsetAdd= 0; %additive offset
    o.MiddleOffsetMult=1; %multiplier offest
    o.TextProps       =[];
    
%parse the inputs
    %line handles
    h = findobj(gca,'Type','Line');
    h = flipud(fliplr(h));
    if ~isempty(varargin)
        if all(ishandle(varargin{1}))
            h = varargin{1};
            varargin(1) = [];
        end
    end
    
    %line labels
    labels = get(h,'DisplayName');
    if ~iscell(labels)
        labels = {labels};
    end
    if ~isempty(varargin)
        if iscell(varargin{1})
            labels = varargin{1};
            varargin(1) = [];
        end
    end 
    
    %Put default parameter names into a variable and setup input parser    
    o = parse_options(o,varargin{:});
        
%check if we're in linear space or logspace
%calculate the minimum Y spacing between each label
    XScale = get(gca,'XScale');
    YScale = get(gca,'YScale');
    XLim   = get(gca,'XLim');
    YLim   = get(gca,'YLim');
    
    switch YScale
        case {'linear'}
            dY           = o.SpacingPct*(YLim(2)-YLim(1));
        case {'log'}
            dY           = o.SpacingPct*(log10(YLim(2)) - log10(YLim(1)));
    end
    dX = o.XOffsetPct*(XLim(2)-XLim(1));
    if strcmp(XScale,'log')
        dX = o.XOffsetPct*(log10(XLim(2))-log10(XLim(1)));
    end
        
%first get the y position of each line
    xpos  = zeros(size(h));
    ypos  = zeros(size(h));
    col   = cell(size(h));
    for i=1:length(h)
        x = get(h(i),'XData');
        y = get(h(i),'YData');
        col{i} = get(h(i),'Color');
        switch o.Position
            case {'left','Left'}
                [xpos(i), ind] = min(x);
                if strcmp(XScale,'linear')
                    xpos(i) = xpos(i)-dX;
                else
                    xpos(i) = 10.^(log10(xpos(i)) - dX);
                end
                ypos(i)       = y(ind);
            case {'right','Right'}
                [xpos(i), ind]= max(x);
                if strcmp(XScale,'linear')
                    xpos(i) = xpos(i)+dX;
                else
                    xpos(i) = 10.^(log10(xpos(i)) + dX);
                end
                ypos(i)       = y(ind);
            case {'middle','Middle','center','Center'}
                if strcmp(XScale,'linear')
                    xpos(i) = (max(x)-min(x))/2;
                else
                    xpos(i) = (max(log10(x))-min(log10(x)))/2;
                    xpos(i) = 10^xpos(i);
                end
                [~, ind] = min(abs(x-xpos(i)));
                ypos(i)  = y(ind);
            otherwise
                error('invalid Position');
        end
    end
    
    if strcmp(YScale,'log')
        ypos = log10(ypos);
    end
    
%sort by position, starting at bottom first
    [ysort, isort] = sort(ypos);
    labelsort     = labels(isort);
    xsort         = xpos(isort);
    csort         = col(isort);
            
%calculate new y positions (called y2)     
    y2sort = zeros(size(ysort));
    for i=1:length(ysort)        
        if i==1
            y2sort(i) = ysort(i);
        else
            if ysort(i)-y2sort(i-1) < dY
                y2sort(i) = y2sort(i-1) + dY;
            else
                y2sort(i) = ysort(i);
            end
        end
    end
    
    if strcmp(YScale,'log')
        y2sort = 10.^(y2sort);
    end

    if ismember(o.Position,{'middle','Middle','center','Center'})
        y2sort = y2sort*o.MiddleOffsetMult + o.MiddleOffsetAdd;
    end    
    
    
%put the labels in the graph
    htxt = zeros(size(h));
    for i=1:length(y2sort)
        xtext = xsort(i);
        if ~isempty(o.XPos)
            xtext = o.XPos(i);
        end
        ytext = y2sort(i);
        if ~isempty(o.YPos) %manually set YPos
            ytext = o.YPos(i);
        end
        htxt(i) = text(xtext,ytext,labelsort{i},'Color',csort{i},'FontWeight',o.FontWeight);
        if ismember(o.Position,{'left','Left'})
            set(htxt(i),'HorizontalAlignment','Right');
        elseif ismember(o.Position,{'middle','Middle','center','Center'})
            set(htxt(i),'HorizontalAlignment','Center');
        end            
    end
    if ~isempty(o.TextProps)
        set(htxt,o.TextProps);
    end    
    
%tweak the XLimits
    if o.XLimTweak==1
        XLim = get(gca,'XLim');
        switch o.Position
            case {'right','Right'}                
                xx = max(x) + (max(x)-min(x))*o.XLimTweakPct;
                if strcmp(XScale,'log')
                    x = log10(x);
                    xx = max(x) + (max(x)-min(x))*o.XLimTweakPct;
                    xx = 10.^(xx);
                end
                xlim([XLim(1) max(xx,XLim(2))]);
            case {'left','Left'}
                xx = min(x) - (max(x)-min(x))*o.XLimTweakPct;
                if strcmp(XScale,'log')
                    x = log10(x);
                    xx = min(x) - (max(x)-min(x))*o.XLimTweakPct;
                    xx = 10.^(xx);
                end                
                xlim([min(xx,XLim(1)) XLim(2)]);
        end
    end
    
%outputs
if nargout==1
    varargout{1} = htxt;
end
