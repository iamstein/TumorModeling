function[dY] = ode_iv_2cmtc_full(~,Y,p)

Ac     = Y(1);
Ap     = Y(2);
T      = Y(3);
DT     = Y(4);

D      = Ac/p.Vc;

         %TRANSIT       		%DOSE/SYN  %CLEARANCE TERMS
dAc    =  -p.k12*Ac + p.k21*Ap              - p.keD *Ac  + (- p.kon*D*T + p.koff*DT)*p.Vc;
dAp    =   p.k12*Ac - p.k21*Ap;
dT     =                 		p.ksyn      - p.keT * T     - p.kon*D*T + p.koff*DT;
dDT    =                            		- p.keDT*DT     + p.kon*D*T - p.koff*DT;   
		
dY     = [dAc; dAp; dT; dDT];