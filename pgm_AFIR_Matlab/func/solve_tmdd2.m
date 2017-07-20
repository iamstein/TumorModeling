function[out Ctot Y] = solve_tmdd2(model,tsim,dose,p)
%[out C] = solve_tmdd(model_string,init,tsim,dose,p);
%model        = model structure, containing model.ode, model.init, model.dcmt, and model.out
%				which specify the ODE, initial condition, dosing compartment, and output variables
%				model.init can either be a vector or a function of the model parameters
%tsim         = time of simulations
%dose         = structure containing dose information for n doses
%               dose.type= 'bolus' or 'infusion'  -- if empty, the default is bolus      
%				dose.t   = row vector of times (1xn)
%               dose.tstop=row vector of stop times for infusion
%			 	dose.d   = row vector of doses (1xn)
%				dose.cmt = compartments in which the dose goes (usually 1xn, but could be mxn, in case of complex absorption)
%                          if empty, the default is the first compartment
%				dose.frac= fraction of dose that goes into each of the compartments (usually ones(1,n) except for complex absorption)
%                          if empty, the default is all in one compartment
%p			  = structure of parameters
%               p.F    will affect the bioavailability of the dose in solve_tmdd2
%               p.tlag will affect the lag time of the dose in solve_tmdd2
%				all other parameters inform
%
% EXAMPLE CODE
%{
    %set up the model
    tsim   = -5:.1:56; %time for simulation
    dose.t = 0;       %dosing times
    dose.d = 4;       %dose amount
    P      = get_siltuximab_Wang14_param(); 
    p      = P(1,:);  %model parameters
    
    ax.XLim   = [-7   60];  %axis properties
    ax.YLim   = [1e-6 1e2];
    ax.XTick  = 0:7:56;
    
    subplot(1,2,1)
        model  = obj_iv_2cmt_qe_Cen;
        out    = solve_tmdd2(model,tsim,dose,p); %solve model
        feval(model.plot,out); %plot model
        title('Standard QE model');
        set(gca,ax);

    subplot(1,2,2)
        model  = obj_iv_2cmt_qe_Cen_Wang14;
        out    = solve_tmdd2(model,tsim,dose,p); %solve model
        feval(model.plot,out); %plot model
        title('Wang14 QE model');
        set(gca,ax);
%}

%reparameterize the model accordingly, in terms of rate constants
    if isfield(model,'pode')
        p = feval(model.repar,p);
        p = strip_param(p,model.pode);  
    else
        warning('no model repar function.  risk of redundant parameters\n');
    end

%set initial condition
    if isnumeric(model.init)
    	init = model.init;
    else
        init = feval(model.init,p);
    end

%set dosing compartment and dosing properties
    if ~isfield(dose,'cmt')
        dose.cmt  = ones(size(dose.t))*model.dcmt;
    end
    if ~isfield(dose,'frac')
        dose.frac = ones(size(dose.t));
    end
    if ~isfield(dose,'type')
        dose.type = 'bolus';
    end
    if strcmp(dose.type,'bolus') && ~isfield(dose,'tstop')
        dose.tstop = NaN;
    end
    if ~isfield(p,'dose_inf') %infusion dose, used in some ode functions
        p.dose_inf = 0;
    end
    if ~isfield(p,'fu') %free fraction of drug, used in some ode functions
        p.fu = 1;
    end

%get variable names in p
    switch(class(p))
        case {'dataset'}
            varnames = get(p,'VarNames');
        case {'struct'}
            varnames = fieldnames(p);
        otherwise
            error('invalid class for p');
    end

%make sure everything is in column vector form
    tsim       = tsim(:);
    dose.d     = dose.d(:);
    dose.t     = dose.t(:);
    dose.tstop = dose.tstop(:);
    
%modify dose matrix in case of infusion
%includes zero dosing after infusion stops
    if strcmp(dose.type,'infusion')
        t          = [dose.t' dose.tstop'];
        t          = t(:);
        dose.t     = t;
        d          = [dose.d' zeros(size(dose.d))'];
        dose.d     = d(:);
    else
        t          = dose.t(:);
    end

%offset dose time if there's a tlag parameter
    if ismember('tlag',varnames)
        dose.t = dose.t + p.tlag;
    end

%break time into sections based on dosing change
    ndose      = find(t<=tsim(end),1,'last'); %loop through doses (only those before the last time point sampled)
    TT{1}      = union(tsim(tsim<=t(1)),t(1));    
    TT{1}      = TT{1}(:);
    for i=1:ndose-1
        TT{i+1} = union(t([i i+1]), tsim(tsim>=t(i) & tsim<=t(i+1)));
    end
    TT{ndose+1}= union(t(ndose),tsim(tsim>=t(ndose)));
    
%ode solver options
    AbsTol = 1e-14;
    opt = odeset('NonNegative',ones(size(init)),'AbsTol',AbsTol);
        
%solve the model before the first dose
    T  = TT{1};
    y0 = init;    
    if length(TT{1})>1
        p.t0  = TT{1}(1);
        [~,y] = ode15s(model.ode,TT{1},y0,opt,p);
        Y     = y;
        y0    = Y(end,:); %reset initial condition for next section
    else
        Y     = y0(:)';
    end
    
%solve model after the first dose
    for i=1:length(TT)-1
        tspan    = TT{i+1};
        
        %for bolus dose, add dose to initial condition        
            if strcmp(dose.type,'bolus')
                for j=1:size(dose.cmt,1)
                    idose = dose.cmt(j,i);
                    y0(idose) = y0(idose) + dose.frac(j,i)*dose.d(i); %add dose
                end
            end
                
        if ~(i==ndose && length(tspan)==1)
            if strcmp(dose.type,'infusion')
                p.dose_inf = dose.d(i);
            end
            
            p.t0    = TT{i+1}(1);
            [~, y]  = ode15s(model.ode,tspan,y0,opt,p);
            
            %because we can get a warning when parameters are difficult
            %and ode will spit out a value for y shorter than tspan, we
            %fill the end of y with NaNs so the function doesn't break
                if size(y,1)<length(tspan)
                    y(end+1:length(tspan),:) = NaN;
                end
                        
            %because when length(tspan)==2, default is to return a bunch of
            %times in between, we need this bit to shrink y back down
                if length(tspan)==2
                    y = y([1 end],:);
                end
            
            %output variables
                Y       = [Y; y]; %#ok<*AGROW>
                T       = [T; TT{i+1}];
                y0      = Y(end,:);
        end
    end
    
%reduce output to be the size of the time vector T
    [~, ii] = unique(T);
    T       = T(ii,:);
    Y       = Y(ii,:);

    ii      = ismember(T,tsim);
    T       = T(ii,:);
    Y       = Y(ii,:);
        
%compute outputs based on model type
    out     = feval(model.out,Y,p);
    Ctot    = out.Ctot;
    out.t   = T;
