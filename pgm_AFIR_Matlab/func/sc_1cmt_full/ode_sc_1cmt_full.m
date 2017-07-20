function[dY] = ode_sc_1cmt_full(~,Y,p)

Adepot = Y(1);
A      = Y(2);
T      = Y(3);
DT     = Y(4);

D      = A/p.V;

         %TRANSIT       %DOSE/SYN  %CLEARANCE TERMS
dAdepot=     -p.ka*Adepot;
dA     =  p.F*p.ka*Adepot               - p.keD * A  + (- p.kon*D*T + p.koff*DT)*p.V;
dT     =                     p.ksyn     - p.keT * T     - p.kon*D*T + p.koff*DT;
dDT    =                                - p.keDT*DT     + p.kon*D*T - p.koff*DT;   

dY     = [dAdepot; dA; dT; dDT];