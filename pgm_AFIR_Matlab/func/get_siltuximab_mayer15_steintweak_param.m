function[p q units] = get_siltuximab_mayer15_steintweak_param()
% From Mayer et al., Cancer Chemother Pharmacol, 2015 siltuximab paper
fprintf('parameters for siltuximab\n');
%from Table 2 for Mayer15 on Siltuximab
    p.target = 'IL-6';

    p.keD  = 0.0584;  %1/d
    p.k12  = 0.138;   %1/d
    p.k21  = 0.187;   %1/d
    p.Vc   = 58.1/1000*70;  %ml/kg --> L
    p.V    = p.Vc;

    p.CL   = p.keD*p.Vc;
    p.Q    = p.k12*p.Vc;
    p.Vp   = p.Q/p.k21;    
    
%from fitting of Pulchalski10 and info from Wang14
    p.ksyn = .005; %For baseline levels of 5pg/ml = 2e-4 nM
    p.keT  = 40.00001; %review of cytokine litearture BiologicsWorkstreams/Work/Task12_Cytokines
    p.keDT = .03; %fit of Charoin10
    p.Kd   = .02;  %Wang14
    
    p.kon  = 10;   %guess
    p.koff = p.Kd*p.kon;
    
    p.T0   = p.ksyn/p.keT;
    p.Tss  = p.ksyn/p.keDT;
    p.Acc  = p.Tss/p.T0;
    
%units
    units.Time      = 'days';
    units.Dose      = 'nmol';
    units.Conc      = 'nM';
    
    units.Clearance = 'L/d';
    units.Volume    = 'L';
    
    units.kin       = 'mg/L/d CRP';
    units.IC50      = 'mg/L CRP';

    q      = micro2macro_bolus_2cmt(p);
