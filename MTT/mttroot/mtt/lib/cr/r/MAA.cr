%SUMMARY MAA: Asymmetrical Mass-action kinetics for Re
%DESCRIPTION Parameter 1: kappa_f Parameter 2: kappa_r

OPERATOR MAA;

%%% Two port version (standard)
FOR ALL n_out, A_f, A_r, kappa_f,kappa_r
LET MAA(R, flow,kappa_f,kappa_r,flow,n_out,	
   A_f,effort,1,	
   A_r,effort,2) = kappa_f*exp(A_f/RT) - kappa_r*exp(A_r/RT);


%%% Four port version (stoichiometric) with integral causality
%% Flows on ports 1 & 2 = flow on port 3.
FOR ALL A_f, A_r, v_1, v_2,kappa_f,kappa_r
LET MAA(R, flow,kappa_f,kappa_r,flow,1,	
   A_f,effort,1,	
   A_r,effort,2,
   v_1,flow,3,
   v_2,effort,4) = v_1;

FOR ALL A_f, A_r, v_1, v_2,kappa_f,kappa_r
LET MAA(R, flow,kappa_f,kappa_r,flow,2,	
   A_f,effort,1,	
   A_r,effort,2,
   v_1,flow,3,
   v_2,effort,4) = v_1;

%% Flow on port 4 is induced flow
FOR ALL A_f, A_r, v_1, v_2,kappa_f,kappa_r
LET MAA(R, flow,kappa_f,kappa_r,flow,4,	
   A_f,effort,1,	
   A_r,effort,2,
   v_1,flow,3,
   v_2,effort,4) = kappa_f*exp(A_f/RT) - kappa_r*exp(A_r/RT);

%%% Four port version (stoichiometric) with flow imposed on port 1
%% Flow on port 2 = flow on port 1
FOR ALL A_1,A_2,A_3,A_4, v_1,v_2,v_3,v_4, kappa_f,kappa_r
LET MAA(R, flow,kappa_f,kappa_r,flow,2,	
   v_1,flow,1,	
   A_2,effort,2,
   v_3,flow,3,
   A_4,effort,4) = v_1;

%% Flow on port 4 = flow on port 1
FOR ALL A_1,A_2,A_3,A_4, v_1,v_2,v_3,v_4, kappa_f,kappa_r
LET MAA(R, flow,kappa_f,kappa_r,flow,4,	
   v_1,flow,1,	
   A_2,effort,2,
   v_3,flow,3,
   A_4,effort,4) = v_1;

%% Effort on port 1
FOR ALL A_1,A_2,A_3,A_4,v_1,v_2,v_3,v_4,kappa_f,kappa_r
LET MAA(R,flow,kappa_f,kappa_r,effort,1,
   v_1,flow,1,	
   A_2,effort,2,
   v_3,flow,3,
   A_4,effort,4) = RT*log( (v_1/kappa_f) + (kappa_r/kappa_f)*exp(A_2/RT) );


%%% Four port version (stoichiometric) with flow imposed on port 2

%% Flow on port 1 = flow on port 2
FOR ALL A_1,A_2,A_3,A_4, v_1,v_2,v_3,v_4, kappa_f,kappa_r
LET MAA(R, flow,kappa_f,kappa_r,flow,1,	
   A_1,effort,1,	
   v_2,flow,2,
   v_3,flow,3,
   A_4,effort,4) = v_2;

%% Flow on port 4 = flow on port 2
FOR ALL A_1,A_2,A_3,A_4, v_1,v_2,v_3,v_4, kappa_f,kappa_r
LET MAA(R, flow,kappa_f,kappa_r,flow,4,	
   A_1,effort,1,	
   v_2,flow,2,
   v_3,flow,3,
   A_4,effort,4) = v_2;

%% Effort on port 2
FOR ALL A_1,A_2,A_3,A_4,v_1,v_2,v_3,v_4,kappa_f,kappa_r
LET MAA(R,flow,kappa_f,kappa_r,effort,2,
   A_1,effort,1,	
   v_2,flow,2,
   v_3,flow,3,
   A_4,effort,4) = RT*log( (v_2/kappa_r) + (kappa_f/kappa_r)*exp(A_1/RT) );

