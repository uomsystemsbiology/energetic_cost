%SUMMARY lMA: Linearised Mass-action kinetics for Re
%% INCOMPLETE!! EXPERIMENTAL!!
%DESCRIPTION Parameter 1: kappa

OPERATOR lMA;

%%% Two port version (standard)
FOR ALL n_out, A_f, A_r, v_p, v_m, RT
LET lMA(R, v_p,v_m, RT,flow,n_out,	
   A_f,effort,1,	
   A_r,effort,2) = (v_p/RT)*A_f - (v_m/RT)*A_r;


%%% Four port version (stoichiometric) with integral causality
%% Flows on ports 1 & 2 = flow on port 3.
FOR ALL A_f, A_r, v_1, v_2, v_p, v_m, RT
LET lMA(R, v_p,v_m, RT,flow,1,	
   A_f,effort,1,	
   A_r,effort,2,
   v_1,flow,3,
   v_2,effort,4) = v_1;

FOR ALL A_f, A_r, v_1, v_2, v_p, v_m, RT
LET lMA(R, v_p,v_m, RT,flow,2,	
   A_f,effort,1,	
   A_r,effort,2,
   v_1,flow,3,
   v_2,effort,4) = v_1;

%% Flow on port 4 is induced flow
FOR ALL A_f, A_r, v_1, v_2, v_p, v_m, RT
LET lMA(R, v_p, v_m, RT ,flow,4,	
   A_f,effort,1,	
   A_r,effort,2,
   v_1,flow,3,
   v_2,effort,4) = (v_p/RT)*A_f - (v_m/RT)*A_r;

%% Need to include and modify rest of MA.cr here.

%%% Five port version (stoichiometric) with integral causality
%% Flow on port 1  = flow on port 3.
FOR ALL A_f, A_r, v_f, v_r, v, v_p,v_m, RT
LET lMA(R, v_p,v_m, RT,flow,1,	
   A_f,effort,1,	
   A_r,effort,2,
   v_f,flow,3,
   v_r,flow,4,
   v,effort,5) = v_f;

%% Flow on port 2  = flow on port 4.
FOR ALL A_f, A_r, v_f, v_r, v, v_p,v_m, RT
LET lMA(R, v_p,v_m, RT,flow,2,	
   A_f,effort,1,	
   A_r,effort,2,
   v_f,flow,3,
   v_r,flow,4,
   v,effort,5) = v_r;

%% Flow on port 5 is induced flow
FOR ALL A_f, A_r, v_f, v_r, v, v_p,v_m, RT
LET lMA(R, v_p,v_m, RT,flow,5,	
   A_f,effort,1,	
   A_r,effort,2,
   v_f,flow,3,
   v_r,flow,4,
   v,effort,5) = (v_p/RT)*A_f - (v_m/RT)*A_r; 


