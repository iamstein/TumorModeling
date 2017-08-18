% This checks the indistinguishability of the 3 component model of
% AFIR_tissue_01.pdf. We want keD3 = 0 so that we can derive AFIRT
% similar to the 2 component model. We derived the new model that is
% indistinguishable from the orignal model by the use of the Similarity
% Transform from godfrey_1989_problem_indisting_pharm.

% Reference
% Godfrey, K.R., Chapman, M.J., 1989. The problem of model 
% indistinguishability in pharmacokinetics. Journal of Pharmacokinetics 
% and Pharmacodynamics 17, 229–267.


function test_indistinguishability_3comp

close all;

% Parameters
VC = 5; VP = 2; VT = 1; 
Dose = random('unif',0,1); 
ksyn1 = random('unif',0,1); ksyn3 = random('unif',0,1); 

k12D = random('unif',0,1); k21D = random('unif',0,1); 
k13D = random('unif',0,1); k31D = random('unif',0,1); 

k13T = random('unif',0,1); k31T = random('unif',0,1); 
k13DT = random('unif',0,1); k31DT = random('unif',0,1); 

kon1 = random('unif',0,1); koff1 = random('unif',0,1); 
kon3 = random('unif',0,1); koff3 = random('unif',0,1); 

keD1 = random('unif',0,1); keD3 = random('unif',0,1); 
keT1 = random('unif',0,1); keT3 = random('unif',0,1); 
keDT1 = random('unif',0,1); keDT3 = random('unif',0,1);

pars = [VC; VP; VT; Dose; ksyn1; ksyn3; k12D; k21D; ...
        k13D; k31D; k13T; k31T; k13DT; k31DT; ...
        kon1; koff1; kon3; koff3; ...
        keD1; keD3; keT1; keT3; keDT1; keDT3];

% Initial value;
D0 = [1; 0; 0];
% Initial time; Final time;
t0 = 0; tf = 50;
% Points per minute; Number of points; Time vector;
ppd = 5; N = (tf-t0)*ppd; tspan = linspace(t0,tf,N);

% Solve ode
[t,Dmod1] = ode45(@RHSmod1,tspan,D0,[],pars);
[~,Dmod2] = ode45(@RHSmod2,tspan,D0,[],pars);
% D rows are the time points. D columns are the variables.
% Ex.: D(a*ppd,:) = [D1(a), D2(a), D3(a)]

% Plot
plot(t,Dmod1,'o',t,Dmod2)
legend('D1','D2','D3','D1''','D2''','D3''')

end

% -------------------------------------------------------------------------
% ODE Right hand side 
% -------------------------------------------------------------------------

function dD = RHSmod1(t,D,pars)

% Retrieve parameters by name
VC = pars(1);     VP = pars(2);     VT = pars(3); 
Dose = pars(4); 
ksyn1 = pars(5);  ksyn3 = pars(6); 

k12D = pars(7);   k21D = pars(8); 
k13D = pars(9);   k31D = pars(10); 

k13T = pars(11);  k31T = pars(12); 
k13DT = pars(13); k31DT = pars(14); 

kon1 = pars(15);  koff1 = pars(16); 
kon3 = pars(17);  koff3 = pars(18); 

keD1 = pars(19);  keD3 = pars(20); 
keT1 = pars(21);  keT3 = pars(22); 
keDT1 = pars(23); keDT3 = pars(24);

% Coefficient matrix A (linear model)
Amod1 = [-(k13D + k12D + keD1), VP/VC*k21D, VT/VC*k31D    ;
         VC/VP*k12D           , -k21D     , 0             ;
         VC/VT*k13D           , 0         , -(k31D + keD3)];
B = [1/VC; 0; 0];

% D(1) = D1; D(2) = D2; D(3) = D3;
dD = Amod1*D + B*Dose;

end

function dD = RHSmod2(t,D,pars)

% Retrieve parameters by name
VC = pars(1);     VP = pars(2);     VT = pars(3); 
Dose = pars(4); 
ksyn1 = pars(5);  ksyn3 = pars(6); 

k12D = pars(7);   k21D = pars(8); 
k13D = pars(9);   k31D = pars(10); 

k13T = pars(11);  k31T = pars(12); 
k13DT = pars(13); k31DT = pars(14); 

kon1 = pars(15);  koff1 = pars(16); 
kon3 = pars(17);  koff3 = pars(18); 

keD1 = pars(19);  keD3 = pars(20); 
keT1 = pars(21);  keT3 = pars(22); 
keDT1 = pars(23); keDT3 = pars(24);

% Coefficient matrix A (linear model)
Amod2 = [-(k13D + k12D + keD1)        , VP/VC*k21D, VT/VC*(keD3 + k31D)                        ;
         VC/VP*(k12D - k13D*keD3/k21D), -k21D     , VT/VP*(k31D*keD3/k21D + keD3^2/k21D - keD3);
         VC/VT*k13D                   , 0         , -(k31D + keD3)                             ];
B     = [1/VC; 0; 0];

% D(1) = D1; D(2) = D2; D(3) = D3;
dD = Amod2*D + B*Dose;

end















