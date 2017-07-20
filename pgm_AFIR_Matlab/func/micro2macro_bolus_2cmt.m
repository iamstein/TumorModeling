function[pout] = micro2macro_bolus_2cmt(p)
%converts pk parameters (CL,V2,V3,Q) to equation parameters (C1,C2,lam1,lam2)

    if ~isstruct(p)
        p = dataset2struct(p);
    end

    %take parameters out of the structure
        f   = fieldnames(p);
        if all(ismember({'V2','V3'},f))
            V2 = p.V2;
            V3 = p.V3;
        elseif all(ismember({'Vc','Vp'},f))
            V2 = p.Vc;
            V3 = p.Vp;
        elseif all(ismember({'Vp','Vt'},f))
            V2 = p.Vp;
            V3 = p.Vt;
        elseif all(ismember({'V1','V2'},f)) && ~ismember('V3',f)
            V2 = p.V1;
            V3 = p.V2;
        elseif ismember('Vc',f)
            V2 = p.Vc;
        end
        if isfield(p,'FD')
            FD = p.FD;
        else
            FD = 1;
        end

    %convert parameters to rate constants if needed
    if isfield(p,'CL')
        CL  = p.CL;
    end
    if isfield(p,'Q')
        Q   = p.Q  ;
    end
    if isfield(p,'ke')
        ke  = p.ke;
    else
        ke  = CL/V2;      
    end
    if isfield(p,'k12') && isfield(p,'k21')
        k12 = p.k12;
        k21 = p.k21;
    elseif isfield(p,'kcp') && isfield(p,'kpc')
        k12 = p.kcp;
        k21 = p.kpc;
            
    else
        k12 = Q  /V2;
        k21 = Q  /V3;
    end
        
    %computer analytical solution constants
        lam1= .5*(k12+k21+ke + sqrt( (k12+k21+ke)^2 - 4*ke*k21));
        lam2= .5*(k12+k21+ke - sqrt( (k12+k21+ke)^2 - 4*ke*k21));
        
        pout.C1  = FD*(k21-lam1)/(V2*(lam2-lam1));
        pout.C2  = FD*(k21-lam2)/(V2*(lam1-lam2));
        pout.C1p = FD* k12/(V2*(lam2-lam1)); %peripheral 1
        pout.C2p = FD* k12/(V2*(lam1-lam2)); %peripheral 2
        pout.lam1= lam1;
        pout.lam2= lam2;
        pout.FD  = FD;
        1;