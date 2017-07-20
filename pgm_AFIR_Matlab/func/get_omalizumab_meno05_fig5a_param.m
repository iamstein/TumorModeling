function[p q units] = get_omalizumab_hayashi07_param()

%parameters from table
    fprintf('parameters for omalizumab from Hayashi 2007\n');
    units.Volume    = 'L';
    units.Time      = 'day';
    units.Dose      = 'nmol';
    units.Clearance = 'L/day';
    units.Vmax      = 'nmol/d';
    units.Conc      = 'mg/L=ug/ml';
    
	p.drug  = 'omalizumab';
    p.target= 'IgE';
    %{
    CL_X	= 7.32/1000*24; %ml/h --> L/d - clearance of drug
    dCL_C   = 5.86/1000*24; %ml/h --> L/d - delta clearance of complex (CL_C = CL_X + dCL_C)
    CL_C    = CL_X + dCL_C;
    CL_E    = 71.0/1000*24; %ml/h --> L/d - clearance of IgE
    
    V_X     = 5900/1000;    %ml--> L volume of drug
    V_E     = V_X;          %assumed to be same
    V_C     = 3630/1000;    %ml--> L volume of complex
    
    PE      = 30.3*1000*24/190e3/V_E; %ug/h --> nM/d
    
    ka      = .02*24;       %1/h --> 1/d
    Kd      = 1.07;         %nM
    %}
    
    %hand tuning parameters
    p.F     = .42;
    p.ka    = .35;
    p.V     = 3;
    p.ksyn  = 1.4;
    p.keD   = 0.005;
    p.keT   = .9333;
    p.keDT  = .2;
    p.koff  = 23.3;
    p.kon   = 10;
    p.Kd    = p.koff/p.kon;

    p.CL    = p.keD*p.V;
    
    p.T0   = p.ksyn/p.keT;
    p.Tss  = p.ksyn/p.keDT;
    p.Acc  = p.Tss/p.T0;

    q      = [];
