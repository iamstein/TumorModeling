function[S] = prettynum(NUM,nsig,nmin)
%prettynum(n) = makes a nice string with only nsig significant digits (default n is 1)
%below nmin, just use %1.(nsig-1)g format string

if nargin<=1
    nsig = 1;
end
if nargin<=2
    nmin = 1e-3;
end

S = cell(size(NUM));
for i=1:length(NUM)
    num = NUM(i);

    p        = floor(log10(num));
    numpref  = roundb(num/10.^p,-(nsig-1));
    numprint = numpref*10^p;

    if numprint>=nmin && 1/numprint>=nmin
        format   = ['%1.' num2str(max(0,nsig-p-1)) 'f'];
        s        = sprintf(format,numprint);
    else
        format1  = ['%1.' num2str(max(0,nsig-1)) 'f'];    
        s        = sprintf([format1 'e%1.0f'],numpref,p);
        1;
    end

    if s(end)=='0' && any(s=='.')
        s = s(1:end-1);
    end
    if s(end)=='.'
        s = s(1:end-1);
    end
    if num == 0
        s = '0';
    end    
    if s(end)=='e';
        s = s(1:end-1);
    end
    S{i} = s;   
end

if length(S)==1
    S = S{1};
end
