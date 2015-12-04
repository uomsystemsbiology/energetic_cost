function V = stoich2flow (X,stoich)
	 
  ## usage: V = stoich2flow (X,stoich)
  ## 
  ## Computes flow V from state X using stoichimetry

  ## Nonlinear version
  kappa = [diag(stoich.kappa) -diag(stoich.kappa)];
  K = diag(stoich.K_c);
  N_fr = stoich.N_fr;
  V = kappa*exp(N_fr'*log(K*X));
	 
endfunction
