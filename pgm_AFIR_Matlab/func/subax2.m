function[] = subax2(Nr,Nc,ir,ic,varargin)
%a version of subaxis with default value (so i don't have to set
%everything all the time)

o.Spacing = 0;
o.PaddingLeft = 0.001;
o.PaddingRight = o.PaddingLeft;
o.PaddingTop  = .04;
o.PaddingBottom = 0;
o.MarginTop = .15;
o.MarginBottom = .15;
o.MarginLeft = .15;
o.MarginRight = .01;

o = parse_options(o,varargin{:});

if nargin==0
    Nr = 1;
    Nc = 1;
    ir = 1;
    ic = 1;
end
    
iplot = Nc*(ir-1)+ic;

subaxis(Nr,Nc,iplot,'Spacing',o.Spacing,...
            'MarginLeft',o.MarginLeft,'MarginRight',o.MarginRight,'MarginTop',o.MarginTop,'MarginBottom',o.MarginBottom,...
            'PaddingTop',o.PaddingTop,'PaddingBottom',o.PaddingBottom,'PaddingLeft',o.PaddingLeft,'PaddingRight',o.PaddingRight);        
        1;
