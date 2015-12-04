function [dx, v] = stoich2ss_fun (x,t)
  ## usage: dx = stoich2ss_fun (x,t)
  ## 
  ## 
	 
  global STOICH2SS

  Kappa = diag(STOICH2SS.kappa);
  NNf = STOICH2SS.NNf;
  NNr = STOICH2SS.NNr;
  KK = STOICH2SS.KK;
  N = STOICH2SS.N;
  x_e = STOICH2SS.x_e;
  L_Xx = STOICH2SS.L_Xx;
  L_xX = STOICH2SS.L_xX;
  G_X = STOICH2SS.G_X;
  X_0 =  STOICH2SS.X_0;


  ## Create full state vector
  X_i = L_Xx*x + G_X*X_0;		# Internal state
  X = [X_i;x_e];			# Composite state

  ## Flows
  v_plus  = Kappa*exp(NNf'*log(KK.*X));
  v_minus = Kappa*exp(NNr'*log(KK.*X));
  v = v_plus - v_minus;

  ## State derivative
  dx = L_xX*N*v;

endfunction
