function[] = fig2file2(filepref, size)
%fig2file - uses matlab - central export command to export figure to file
%size should be in inches

set(gcf,'Units','inches');
set(gcf,'PaperUnits','inches')

set(gcf,'PaperPosition',[0 0 size])
set(gcf,'PaperSize',size)

if exist('../results','dir')
    fdir = '../results/';
elseif exist('./results','dir')
    fdir = './results/';
else
	warning('figs being output to current working directory');
    fdir = './';
end

a = ver;
yr = str2double(a(1).Release(3:7));

if yr>2012
    col = get(gcf,'Color');
    set(gcf,'Color','none')
end

%print([fdir filepref '.eps'],'-depsc2');
print([fdir filepref '.pdf'],'-dpdf');

if yr>2012
    set(gcf,'Color',col);
end