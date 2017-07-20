function[p] = repar_kappa2koff(p)

Acc    = p.keT/p.keDT;
Kd     = p.kappa/Acc;
p.koff  = p.kon*Kd;
