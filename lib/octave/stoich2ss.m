function [X_ss V_ss x_ss] = stoich2ss (Nf,Nr,kappa,K,Nef,Ner,Ke,X_0,L_xX,L_Xx,G_X)
  ## usage: [X_ss V_ss] = stoich2ss (Nf,Nr,kappa,K,Nef,Ner,Ke,X_0,L_xX,L_Xx,G_X)
  ## 
  ## 

  warning("stoich2ss is deprecated - please use stoich_sim with t=0 instead");
  x_e = ones(size(Ke));
 
  global STOICH2SS	 

  STOICH2SS.kappa = kappa;
  STOICH2SS.NNf = [Nf;Nef];
  STOICH2SS.NNr = [Nr;Ner];
  STOICH2SS.KK = [K;Ke];
  STOICH2SS.N = Nr-Nf;
  STOICH2SS.x_e = x_e;
  STOICH2SS.X_0 = X_0;		# State
  STOICH2SS.L_xX = L_xX;
  STOICH2SS.L_Xx = L_Xx;
  STOICH2SS.G_X = G_X;

  ## Start vector - reduced order
  x_0 = L_xX*X_0
  x_ss = fsolve("stoich2ss_fun",x_0);

  ## Regenerate state vector
  X_ss = L_Xx*x_ss + G_X*X_0;

  # Compute ss flows.
  [dx,V_ss] =  stoich2ss_fun (x_ss);

endfunction
