function[p q units] = get_bevacizumab_elmeliegy15_steintweak_param()

%parameters from table
    fprintf('parameters for bevacizumab\n');
    units.Volume    = 'L';
    units.Time      = 'day';
    units.Dose      = 'nmol';
    units.Clearance = 'L/day';
    units.Vmax      = 'nM/day';
    units.Conc      = 'nM';
    
    p.drug  = 'bevacizumab';
    p.target= 'VEGF';
    
    %Table 3 from panoilia15
        CL         = 0.2;    %L/day    
        V1         = 3.14;        %L
        V2         = 2.36;     %L
        Q          = .36;      %L/day

        p.CL       = CL;
        p.Vc       = V1;
        p.Vp       = V2;
        p.Q        = Q;
        
    %panoilia15 binding parameters    
        %BM0        = 212;      %ng/L
        %kout       = .401;     %1/d
        %Kss        = 267;      %nM;
        
        p.Vc       = V1;       %L
        p.ke       = CL/V1;    %1/day
        p.kcp      = Q/V1;     %1/day
        p.kpc      = Q/V2;     %1/day
    
        p.keD      = p.ke;
        p.k12      = p.kcp;
        p.k21      = p.kpc;
                
	%some binding parameters from elmeliegy
        %ksynVEGF   = 27.12;   %pg/ml/d
        %kdegVEGF   =  0.55;   %1/d
        %EC50       = 29.88;   %ug/ml
        
        T0   = .002;
        Tss  = .2;
        keDT = .07;
           
        p.keDT = keDT;
        p.ksyn = Tss*p.keDT;
        p.keT  = p.ksyn/T0;   
        
        p.koff  = 36;
        p.kon   = 2;
        p.Kd    = p.koff/p.kon; %nM
        
    %tweaking to fit data
        p.T0       = T0;
        p.Tss      = Tss;
        p.Acc      = Tss/T0;
        
    q      = micro2macro_bolus_2cmt(p);

        
