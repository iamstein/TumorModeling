function[x] = dataset_iscol(ds,colname)
%dataset_iscol(ds,colname)
%check to see if colname is a column of dataset ds

x = false;
if strcmp(class(ds),'dataset')
    varnames = get(ds,'VarNames');
    x = ismember(colname,varnames);
end

