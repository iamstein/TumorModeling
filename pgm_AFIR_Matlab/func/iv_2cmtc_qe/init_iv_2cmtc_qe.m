function[init] = init_iv_2cmtc_qe(p)
%default initial condition

       %Ac, Ap, T   : AmtCentral, AmtPeriph, TotalTarget
init = [0   0   p.ksyn/p.keT];