function[out] = out_iv_2cmtc_full(Y,p)
%convert solve_tmdd2 output matrix (Y) into structure of output vars

Ac  	 = Y(:,1);
Ap    	 = Y(:,2);
T        = Y(:,3);
DT       = Y(:,4);

C    	 = Ac/p.Vc;
D        = C;

out.Ctot = D + DT;
out.Ap   = Ap;
out.Dtot = D + DT;
out.Ttot = T + DT;
out.D    = D;
out.T    = T;
out.DT   = DT;

try
    out.Tpchg= out.T/(p.ksyn/p.keT);
end