function [N_f,N_r] = N2Nfr (N)
	 
## usage: [N_f,N_r] = N2Nfr (N)
## 
## Get forward and backward stoichiometric matrices from
## stoichiometric matrix

  aN = abs(N);
  N_f = aN.*(N<0.5);
  N_r = aN.*(N>0.5);
	 
endfunction
