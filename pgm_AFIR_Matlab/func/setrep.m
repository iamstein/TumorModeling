function[] = setrep(h,s)
%SETREP(h,s) - replaces all fields of h with s
%h = a figure handle
%s = structure containing some figure properties

if ~isempty(s)
    set(h,s);
end
