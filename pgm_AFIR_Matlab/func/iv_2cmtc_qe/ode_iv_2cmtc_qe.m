function[dY] = ode_iv_2cmtc_qe(~,Y,p)

Actot  = Y(1);
Ap     = Y(2);
Ttot   = Y(3);

Dtot   = Actot/p.Vc;  %concentration in central compartment

B      = Dtot  - Ttot - p.Kd;
D      = 0.5*( B + sqrt(B^2 + 4*p.Kd*Dtot) ); % free drug concentration in central
Ac     = D*p.Vc;

DT     = Ttot*D/(p.Kd+D); %bound receptor concentration
T      = Ttot - DT;       %free receptor concentration

         %TRANSIT               %DOSE/SYN    %CLEARANCE TERMS
dActot = -p.k12*Ac + p.k21*Ap              - p.keD*Ac - p.keDT*DT*p.Vc;
dAp    =  p.k12*Ac - p.k21*Ap;
dTtot  =                        p.ksyn     - p.keT*T  - p.keDT*DT; %total concentration of target

dY     = [dActot; dAp; dTtot];