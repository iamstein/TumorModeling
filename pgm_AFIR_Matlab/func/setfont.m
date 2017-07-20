function[] = setfont(FontSize,FontName)
%setfont - sets the font of the x and y axis and the title and the axes to
%FontSize

if nargin<2
    FontName = 'Helvetica';
end

a = get(gcf,'Children');
for i=1:length(a)
    b = get(a,'Children');
    if ~iscell(b)
        c{1} = b;
        b    = c;
    end
    for j=1:length(b{i})
        if isprop(b{i}(j),'FontSize')
            set(b{i}(j),'FontSize',FontSize)
        end
        1;
    end
    try h(1) = get(a(i),'XLabel'); catch, end
    try h(2) = get(a(i),'YLabel'); catch, end
    try h(3) = get(a(i),'ZLabel'); catch, end
    try h(4) = get(a(i),'Title'); catch, end
    try h(5) = a(i); catch, end
    try set(h,'FontSize',FontSize,'FontName',FontName); catch, end
end

if verLessThan('matlab','8.5')
    a = findobj(gcf,'-propertyname','FontSize');
else
    a = findobj(gcf,'-property','FontSize');
end
for i=1:length(a)
    set(a(i),'FontSize',FontSize,'FontName',FontName);
end
