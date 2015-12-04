%SUMMARY GHK: Goldman-Hodgkin-Katz kinetics 
%DESCRIPTION Parameter 1: kappa

OPERATOR GHK_fun;

%%% AE version
FOR ALL V, V_0
LET GHK(AE, V_0,effort,2,	
   V,effort,1) =  RT*log(GHK_fun(V,V_0));

%%% Three port version (standard)
%% Ports 1 & 2

FOR ALL A_f, A_r, kappa, V, V_0
LET GHK(R, flow,kappa,V_0,flow,1,	
   A_f,effort,1,	
   A_r,effort,2,
   V  ,effort,3) = kappa*GHK_fun(V,V_0)*(exp(A_f/RT)-exp(A_r/RT));

FOR ALL A_f, A_r, kappa, V, V_0
LET GHK(R, flow,kappa,V_0,flow,2,	
   A_f,effort,1,	
   A_r,effort,2,
   V  ,effort,3) = kappa*GHK_fun(V,V_0)*(exp(A_f/RT)-exp(A_r/RT));


%% Port 3 - zero flow
FOR ALL A_f, A_r, kappa, V, V_0
LET GHK(R, flow,kappa,V_0,flow,3,	
   A_f,effort,1,	
   A_r,effort,2,
   V  ,effort,3) = 0;


%%% Five port version (stoichiometric) with integral causality
%% Flows on ports 1 & 2 = flow on port 3.
FOR ALL A_f, A_r, v_1, v_2, V, kappa
LET GHK(R, flow,kappa,flow,1,	
   A_f,effort,1,	
   A_r,effort,2,
   v_1,flow,3,
   v_2,effort,4,
   V,effort,5) = v_1;

FOR ALL A_f, A_r, v_1, v_2, V, kappa
LET GHK(R, flow,kappa,flow,2,	
   A_f,effort,1,	
   A_r,effort,2,
   v_1,flow,3,
   v_2,effort,4,
   V,effort,5) = v_1;

%% Flow on port 4 is induced flow
FOR ALL A_f, A_r, v_1, v_2, V, kappa
LET GHK(R, flow,kappa,flow,4,	
   A_f,effort,1,	
   A_r,effort,2,
   v_1,flow,3,
   v_2,effort,4,
   V,effort,5) = kappa*GHK_fun(V)*(exp(A_f/RT)-exp(A_r/RT));
%%kappa*GHK_fun(V,exp(A_f/RT),exp(A_r/RT));

%% Flow on port 5 is zero
FOR ALL A_f, A_r, v_1, v_2, V, kappa
LET GHK(R, flow,kappa,flow,5,	
   A_f,effort,1,	
   A_r,effort,2,
   v_1,flow,3,
   v_2,effort,4,
   V,effort,5) = 0;


