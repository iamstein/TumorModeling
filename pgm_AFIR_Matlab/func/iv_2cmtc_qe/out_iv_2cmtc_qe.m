function[out] = out_iv_2cmtc_qe(Y,p)
%convert solve_tmdd2 output matrix (Y) into structure of output vars

Actot 	 = Y(:,1);
Ap   	 = Y(:,2);
Rtot     = Y(:,3);

Ctot 	 = Actot/p.Vc;
B    	 = Ctot - Rtot - p.Kd;
C    	 = 0.5*( B + sqrt(B.^2 + 4*p.Kd*Ctot) ); % free drug concentration
RC   	 = Rtot.*C./(p.Kd+C);

out.Ctot = Ctot;
out.Rtot = Rtot;
out.Ap   = Ap;
out.C    = C;                    
out.R    = Rtot - RC;
out.RC   = RC; 

out.Dtot = out.Ctot;
out.Ttot = out.Rtot;
out.D    = out.C;
out.T    = out.R;
out.DT   = out.RC;

try
    out.Tpchg= out.T/(p.ksyn/p.keT);
    out.Tocc = 1-out.Tpchg;
end

    