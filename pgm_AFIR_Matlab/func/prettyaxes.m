function[] = prettyaxes(rr,cc,i,n)
%for tight plots (using subax or subaxis) get rid of some x/y labels
%rr = number of rows
%cc = number of columns
%i  = plot index
%n  = total number of plots (default rr*cc)

if nargin<4
    n = rr*cc;
end

if mod(i,cc)~=1 && cc>1
    set(gca,'YTickLabel',{});
    ylabel('');
end
if i<=n-cc && rr>1
    set(gca,'XTickLabel',{});
    xlabel('');
end