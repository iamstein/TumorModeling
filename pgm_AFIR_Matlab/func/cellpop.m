function[c] = cellpop(m,n,stuff)
%CELLPOP(M,N,stuff) - creates cell of size mxn populated with "stuff"

if ischar(m)
    v1 = m;
    v2 = n;
    v3 = stuff;
    
    m = v2;
    n = v3;
    stuff = v1;
end

c = cell(m,n);
c(:) = {stuff};