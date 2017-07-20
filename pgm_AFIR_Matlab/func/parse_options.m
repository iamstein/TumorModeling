function[opt] = parse_options(def,varargin)
%parses the options input into a function
%opt = parse_options(def,varargin{:});
%
%def is a structe containing the default options values
%varargin can either be a structure or a list of
%{ParameterName1, ParameterValue1, ... }
%
%
%example code
%{
    o = []; def = [];
    def.a = 1;
    def.b = 2

    opt1 = parse_options(def,'a','test1')

    o.b  = 'test2';
    opt2 = parse_options(def,o)
%}

if verLessThan('matlab','100.0')
    opt = def;
    if ~isempty(varargin)
        %for a structure
        if length(varargin)==1
            p = varargin{1};
            if isstruct(p)
                f = fieldnames(p);                    
                for i=1:length(f)
                    fi = f{i};                        
                    opt.(fi) = p.(fi);
                end
            elseif ~isempty(p)
                error('varargin of length one should be a structure');
            end
        else
            %for a list of VarName, VarValue
            if mod(length(varargin),2)~=0
                error('varargin of length >1 should be even ParamName, ParamValue');
            end
            for i=1:2:length(varargin)
                name = varargin{i};
                value= varargin{i+1};
                opt.(name) = value;
            end
        end
    end
else
    if ~isempty(varargin)
        p = inputParser;
        f = fieldnames(def);
        for i=1:length(f)
            p.addOptional(f{i},def.(f{i}));
        end
        if length(varargin)==1 %the options are a structure
            p.parse(varargin{1});
        else
            p.parse(varargin{:})
        end
        opt = p.Results;
    else
        opt = def;
    end
end