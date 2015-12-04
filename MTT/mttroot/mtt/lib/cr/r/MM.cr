%SUMMARY MM: Double-sided Michaelis-Menten kinetics for R
%DESCRIPTION Parameter 1: kappa
%DESCRIPTION Parameter 2: k_m

OPERATOR MM;

%%% Four port version (stoichiometric) with integral causality
%% Flows on ports 1 & 2 = flow on port 3.
FOR ALL A_f, A_r, v_1, v_2,kappa,k_m
LET MM(R, flow,kappa,k_m,flow,1,	
   A_f,effort,1,	
   A_r,effort,2,
   v_1,flow,3,
   v_2,effort,4) = v_1;

FOR ALL A_f, A_r, v_1, v_2,kappa,k_m
LET MM(R, flow,kappa,k_m,flow,2,	
   A_f,effort,1,	
   A_r,effort,2,
   v_1,flow,3,
   v_2,effort,4) = v_1;

%% Flow on port 4 is induced flow
FOR ALL A_f, A_r, v_1, v_2,kappa,k_m
LET MM(R, flow,kappa,k_m,flow,4,	
   A_f,effort,1,	
   A_r,effort,2,
   v_1,flow,3,
   v_2,effort,4) = kappa*(exp(A_f/RT) - exp(A_r/RT))/(k_m + (exp(A_f/RT) + exp(A_r/RT))/2 );
