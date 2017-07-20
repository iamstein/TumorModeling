function[p] = repar_Txx2kxx(p)

p.ksyn  = p.keDT*p.Tss;
p.keT   = p.ksyn/p.T0;