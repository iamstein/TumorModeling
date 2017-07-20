function[] = explore_tmdd_Nparamx_AFITFI(PAR,logrange,OUT,P,model,varargin)

o.Color1 = [0 0 .6];
o.Color2 = [.6 0 0];

o.dx     = .5;
o.dxtick = 2;

o.Axis   = [];
o.XType       = 'FoldChange'; %can either be FoldChange or Parameter

o = parse_options(o,varargin{:});

if isstruct(PAR) || isa(PAR,'dataset')
    p_names = PAR.p_names;
else
    p_names = PAR;
end


rr = 2;
cc = length(p_names);

XLabelFlag = 0;
nparam = length(p_names);
for i=1:nparam
    nsens = length(OUT{i});
    pname = p_names{i};
        
    Colorsi    = diverging_map(linspace(0,1,nsens),o.Color1,o.Color2);
    if mod(nsens,2)==0
        Colorsi(n/2+.5,:) = .4*[1 1 1];
    end

    AFIcalc = [];
    AFItheory= [];
    
    for j=1:nsens
        out     = OUT{i}{j};
        p       = P{i}{j};
        
        switch o.XType
        case 'FoldChange'
            xvals = logrange.(pname);
        case 'Parameter'
            xvals(j) = p.(pname);
        end

        tlast   = out.t(end);
        ind     = out.t>=tlast-p.tau & out.t<=tlast;
        dt      = max(out.t(ind)-min(out.t(ind)));            

        pode    = feval(model.repar,p);
        q       = feval(model.fun.Macro,p); %calculate the Macroconstants C, lam, etc.        

        optax.MarginBottom = .2;
        optax.MarginTop    = .1;
        optax.PaddingTop   = 0;
        
        subax2(rr,cc,1,i,optax); hold on           
            AFI     = trapz(out.t(ind),out.Tpchg(ind))/dt;
            plot(xvals(j),AFI,'ko','MarkerFaceColor',Colorsi(j,:),'MarkerSize',8);
            AFIcalc(j)   = AFI; %#ok<*AGROW,*SAGROW>
            AFItheory(j) = feval(model.fun.AFI,pode);
    
        subax2(rr,cc,2,i,optax); hold on           
            TFI          = max(out.Tpchg(ind));
            plot(xvals(j),TFI,'ko','MarkerFaceColor',Colorsi(j,:),'MarkerSize',8);
            TFIcalc(j)   = TFI;
            TFItheory(j) = feval(model.fun.TFI,pode,q);
    end 
    
    subax2(rr,cc,1,i,optax); hold on                   
        plot(xvals,AFItheory,':','Color',[ 0 .7 0],'LineWidth',2);      
        plot(xvals,AFIcalc  ,'-','Color',[.5 .5 .5]);      
        set(gca,o.Axis)
        ylabel('AFI ratio');
        title(pname)    
        prettyaxes(rr,cc,i)


    subax2(rr,cc,2,i,optax); hold on                   
        plot(xvals,TFItheory,':','Color',[ 0 .7 0],'LineWidth',2);      
        plot(xvals,TFIcalc  ,'-','Color',[.5 .5 .5]);      
        set(gca,o.Axis)
        ylabel('TFI ratio');
        if i==round(nparam/2) && XLabelFlag == 0;
            xlabel('Fold Change in Parameter')
            XLabelFlag = 1;
        end
        prettyaxes(rr,cc,i+cc)
end