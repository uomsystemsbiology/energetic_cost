%SUMMARY sMA: Flow source Re component

OPERATOR MA;

%% sRe component - externally imposed flow
%%% Four port version (stoichiometric) with integral causality
%% Flows on ports 1 & 2 = flow on port 3.
FOR ALL A_f, A_r, v_1, v_2,kappa
LET MA(R, flow,1,	
   A_f,effort,1,	
   A_r,effort,2,
   v_1,flow,3,
   v_2,effort,4) = v_1;

FOR ALL A_f, A_r, v_1, v_2,kappa
LET MA(R, flow,2,	
   A_f,effort,1,	
   A_r,effort,2,
   v_1,flow,3,
   v_2,effort,4) = v_1;

%% Flow on port 4 is induced flow
FOR ALL A_f, A_r, v_1, v_2,kappa
LET MA(R, flow,4,	
   A_f,effort,1,	
   A_r,effort,2,
   v_1,flow,3,
   v_2,effort,4) = v_1;
