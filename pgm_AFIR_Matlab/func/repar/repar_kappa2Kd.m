function[p] = repar_kappa2Kd(p)

Acc    = p.keT/p.keDT;
p.Kd   = p.kappa/Acc;
