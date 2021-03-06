%SUMMARY ecp    Chemical potential constitutive relationship - electrogenic capacitor

% This formula gives the Gibb's free energy for an electrogenic capacitor

% C version integral causality
FOR ALL  concentration, k_e, RT
LET ecp(C,k_e,RT, effort, 1, 
        concentration, state, 1)
         = RT*k_e*concentration;

% C version derivative causality
FOR ALL  potential, k_e, RT
LET ecp(C,k_e,RT, state, 1, 
        potential, effort, 1)
         = (potential/RT)/k_e;

