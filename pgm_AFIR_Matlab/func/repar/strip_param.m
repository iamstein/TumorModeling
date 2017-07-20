function[q] = strip_param(p,pnames)
%keeps only the essential pnames, specified by pnames

q = [];
for i=1:length(pnames)
    s = pnames{i};
    q.(s) = p.(s);
end

keep_anyway = {'dose','tau'};
for i=1:length(keep_anyway)
    s = keep_anyway{i};
    if isfield(p,s)
        q.(s) = p.(s);
    end
end
    