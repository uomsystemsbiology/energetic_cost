%SUMMARY icp    inverse chemical potential constitutive relationship

%% mu is the chemical potential
%% Standard version: flow (f) given chem pot (mu)
FOR ALL mu, k, RT
LET icp(R, k, RT, flow, 1, 
        mu, effort, 1)
         = k*exp(mu/RT);

%% Inverse version:  chem pot (mu_ given flow (f)
FOR ALL f, k, RT
LET icp(R, k, RT, effort, 1, 
        f, flow, 1)
         = RT*log(f/k);

