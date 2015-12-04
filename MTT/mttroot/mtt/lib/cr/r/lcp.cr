%SUMMARY lcp    Chemical potential constitutive relationship - linear case

% This formula gives the Gibb's free energy for a reaction component
% with given concentration.
% k_e is the corresponding free-energy constant 

% C version integral causality
FOR ALL  concentration, x_ss,RT
LET lcp(C,x_ss,RT, effort, 1, 
        concentration, state, 1)
         = (RT/x_ss)*concentration;

% C version derivative causality
FOR ALL  potential, x_ss,RT
LET lcp(C,x_ss,RT, state, 1, 
        potential, effort, 1)
         = (potential)/(RT/x_ss);

