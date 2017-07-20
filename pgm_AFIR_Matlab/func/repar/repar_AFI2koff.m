function[p] = repar_AFI2koff(p)

if isfield(p,'Vc')
    CL     = p.keD*p.Vc;
elseif isfield(p,'V')
    CL     = p.keD*p.V;
end
if isfield(p,'F')
    F      = p.F;
else
    F      = 1;
end

Acc    = p.keT/p.keDT;
Davg   = (F*p.dose)/(CL*p.tau);
Kd     = p.AFI*Davg/Acc; 

p.koff = Kd*p.kon;
