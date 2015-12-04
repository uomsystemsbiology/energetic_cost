%SUMMARY DM: Double-sided Michaelis-Menten kinetics for R
%DESCRIPTION Parameter 1: v_max
%DESCRIPTION Parameter 2: k_v 
%DESCRIPTION Parameter 3: rho_v

OPERATOR DM;

%%% Four port version (stoichiometric) with integral causality
%% Flows on ports 1 & 2 = flow on port 3.
FOR ALL A_f, A_r, v_1, v_2,v_max, k_v, rho_v
LET DM(R, v_max, k_v, rho_v,flow,1,	
   A_f,effort,1,	
   A_r,effort,2,
   v_1,flow,3,
   v_2,effort,4) = v_1;

FOR ALL A_f, A_r, v_1, v_2,v_max, k_v, rho_v
LET DM(R, v_max, k_v, rho_v,flow,2,	
   A_f,effort,1,	
   A_r,effort,2,
   v_1,flow,3,
   v_2,effort,4) = v_1;

%% Flow on port 4 is induced flow
FOR ALL A_f, A_r, v_1, v_2,v_max, k_v, rho_v
LET DM(R, v_max, k_v, rho_v,flow,4,	
   A_f,effort,1,	
   A_r,effort,2,
   v_1,flow,3,
   v_2,effort,4) = v_max*( exp(A_f/RT) - exp(A_r/RT) )/
                   ( k_v +  (1-rho_v)*exp(A_f/RT) +  rho_v*exp(A_r/RT) );


%%% Five port version (stoichiometric) with integral causality
%% Flow on port 1  = flow on port 3.
FOR ALL A_f, A_r, v_f, v_r, v, v_max,k_v,rho_v
LET DM(R, v_max,k_v,rho_v,flow,1,	
   A_f,effort,1,	
   A_r,effort,2,
   v_f,flow,3,
   v_r,flow,4,
   v,effort,5) = v_f;

%% Flow on port 2  = flow on port 4.
FOR ALL A_f, A_r, v_f, v_r, v, v_max,k_v,rho_v
LET DM(R, v_max,k_v,rho_v,flow,2,	
   A_f,effort,1,	
   A_r,effort,2,
   v_f,flow,3,
   v_r,flow,4,
   v,effort,5) = v_r;

%% Flow on port 5 is induced flow
FOR ALL A_f, A_r, v_f, v_r, v, v_max,k_v,rho_v
LET DM(R, v_max,k_v,rho_v,flow,5,	
   A_f,effort,1,	
   A_r,effort,2,
   v_f,flow,3,
   v_r,flow,4,
   v,effort,5) = v_max*( exp(A_f/RT) - exp(A_r/RT) )/
                   ( k_v +  (1-rho_v)*exp(A_f/RT) +  rho_v*exp(A_r/RT) );
