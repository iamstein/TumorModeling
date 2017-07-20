function[p] = repar_kappa2kon(p)

Acc    = p.keT/p.keDT;
Kd     = p.kappa/Acc;
p.kon  = p.koff/Kd;
