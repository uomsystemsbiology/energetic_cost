%SUMMARY mMA: Moodulated Mass-Action kinetics for Re
%DESCRIPTION Parameter 1: i_cr - 

OPERATOR mMA;

%%% Three port version (standard)
FOR ALL n_out, A_f, A_r, m, i_cr
LET mMA(R,i_cr,flow,n_out,	
   A_f,effort,1,	
   A_r,effort,2,
   m,effort,3) = mMA_cr(i_cr,m)*(exp(A_f/RT) - exp(A_r/RT));


%%% Five port version (stoichiometric) with integral causality
%% Flows on ports 1 & 2 = flow on port 3.
FOR ALL A_f, A_r, v_1, v_2, m, i_cr
LET mMA(R, i_cr,flow,1,	
   A_f,effort,1,	
   A_r,effort,2,
   v_1,flow,3,
   v_2,effort,4,
   m,effort,5) = v_1;

FOR ALL A_f, A_r, v_1, v_2, m, i_cr
LET mMA(R, i_cr,flow,2,	
   A_f,effort,1,	
   A_r,effort,2,
   v_1,flow,3,
   v_2,effort,4,
   m,effort,5) = v_1;

%% Flow on port 4 is induced flow
FOR ALL A_f, A_r, v_1, v_2, m, i_cr
LET mMA(R, i_cr,flow,4,	
   A_f,effort,1,	
   A_r,effort,2,
   v_1,flow,3,
   v_2,effort,4,
   m,effort,5) = mMA_cr(i_cr,m)*(exp(A_f/RT) - exp(A_r/RT));

%% Flow on port 5 is modulation flow (0)
FOR ALL A_f, A_r, v_1, v_2, m, i_cr
LET mMA(R, i_cr,flow,5,	
   A_f,effort,1,	
   A_r,effort,2,
   v_1,flow,3,
   v_2,effort,4,
   m,effort,5) =  0;

%% ALL of the following needs to be revised.

%%% Four port version (stoichiometric) with flow imposed on port 1
%% Flow on port 2 = flow on port 1
FOR ALL A_1,A_2,A_3,A_4, v_1,v_2,v_3,v_4, kappa
LET mMA(R, flow,kappa,flow,2,	
   v_1,flow,1,	
   A_2,effort,2,
   v_3,flow,3,
   A_4,effort,4) = v_1;

%% Flow on port 4 = flow on port 1
FOR ALL A_1,A_2,A_3,A_4, v_1,v_2,v_3,v_4, kappa
LET mMA(R, flow,kappa,flow,4,	
   v_1,flow,1,	
   A_2,effort,2,
   v_3,flow,3,
   A_4,effort,4) = v_1;

%% Effort on port 1
FOR ALL A_1,A_2,A_3,A_4,v_1,v_2,v_3,v_4,kappa
LET mMA(R,flow,kappa,effort,1,
   v_1,flow,1,	
   A_2,effort,2,
   v_3,flow,3,
   A_4,effort,4) = RT*log( (v_1/kappa) + exp(A_2/RT) );


%%% Four port version (stoichiometric) with flow imposed on port 2

%% Flow on port 1 = flow on port 2
FOR ALL A_1,A_2,A_3,A_4, v_1,v_2,v_3,v_4, kappa
LET mMA(R, flow,kappa,flow,1,	
   A_1,effort,1,	
   v_2,flow,2,
   v_3,flow,3,
   A_4,effort,4) = v_2;

%% Flow on port 4 = flow on port 2
FOR ALL A_1,A_2,A_3,A_4, v_1,v_2,v_3,v_4, kappa
LET mMA(R, flow,kappa,flow,4,	
   A_1,effort,1,	
   v_2,flow,2,
   v_3,flow,3,
   A_4,effort,4) = v_2;

%% Effort on port 2
FOR ALL A_1,A_2,A_3,A_4,v_1,v_2,v_3,v_4,kappa
LET mMA(R,flow,kappa,effort,2,
   A_1,effort,1,	
   v_2,flow,2,
   v_3,flow,3,
   A_4,effort,4) = RT*log( (v_2/kappa) + exp(A_1/RT) );


%% AE version
FOR ALL mu
LET mMA(AE, effort,2, mu, effort,1) = exp(mu/RT);

;END;
