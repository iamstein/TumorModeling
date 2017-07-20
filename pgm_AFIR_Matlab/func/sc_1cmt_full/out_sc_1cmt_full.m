function[out] = out_sc_1cmt_full(Y,p)
%convert solve_tmdd2 output matrix (Y) into structure of output vars

Ad  	 = Y(:,1);
A    	 = Y(:,2);
T        = Y(:,3);
DT       = Y(:,4);

C    	 = A/p.V;
D        = C;
Ctot     = D + DT;

out.Ctot = D + DT;
out.Dtot = D + DT;
out.Ttot = T + DT;
out.D    = D;
out.T    = T;
out.DT   = DT;
out.Tpchg= out.T/(p.ksyn/p.keT);