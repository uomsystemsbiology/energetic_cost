%SUMMARY gate_fun: Gate modulation function
%DESCRIPTION Parameter 1: kappa

OPERATOR gate_fun_fun;

%%% AE version
FOR ALL FV, i_cr
LET gate_fun(AE, i_cr ,effort,2,	
   FV,effort,1) =  RT*log(mma_cr(i_cr,FV));

