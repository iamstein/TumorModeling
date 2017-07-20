function[p] = repar_CLV2kxx(p)

p.keD = p.CL/p.Vc;
p.k12 = p.Q/p.Vc;
p.k21 = p.Q/p.Vp;