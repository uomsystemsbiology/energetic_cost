function [dx V] = stoich_sim_fun (x,t)
  
  ## usage:  dx = stoich_sim_fun (x,t)
  ##
  ## Function for simulation of reduced-order stoichiometric equations
  ## global  STOICH contains simulation information.
  ## Simulate using stoich_sim
  
  global STOICH

  if nargin<2
    t = STOICH.t;
  endif


  name = STOICH.name;
  N = STOICH.N;
  L = STOICH.L;
  G_X = STOICH.G_X;
  L_xX = STOICH.L_xX;
  X_0 = STOICH.X_0;
  par = STOICH.par;
  ##I_free = STOICH.I_free;
  J_free = STOICH.J_free;

  ## External inputs via fixed states
  X0_fun = STOICH.X0_fun;
  if (length(X0_fun)>0)
     i_fix = STOICH.i_fix;
     X_fix = eval(sprintf("%s(X_0,t)",X0_fun));
     if (length(X_fix)==length(i_fix))
       X_0(i_fix) = X_fix;
     else
	 error(sprintf("function X_fix = %s: X_fix must have %i elements", ...
				X0_fun, length(i_fix)));
     endif
  endif


  X = L*x + G_X*X_0;		# Convert x to X

  ## Avoid numerical issues !!??
  ## small = 1e-10;
  ## i_neg = find(X<small);
  ## X(i_neg) = small;

##  u = eval(sprintf("%s_input(X,[],t,par);",name));
  u = zeros(100,1);		# Dummy
  V = eval(sprintf("%s_odeo(X,u,t,par);",name));
  
  V = J_free*V;			# Set some flows to zero.

  ## External inputs via fixed flows
  V0_fun = STOICH.V0_fun;
  if (length(V0_fun)>0)
    j_fix = STOICH.j_fix;
    V(j_fix) = eval(sprintf("%s(X,t)",V0_fun));
  endif

  ## Enzyme modulation
  enzyme_fun = STOICH.enzyme_fun;
  if length(enzyme_fun)>0
     U = eval(sprintf("%s(X,t)",enzyme_fun));
     V = U.*V;
  endif

  ##dx = L_xX*I_free*N*V;
  dx = L_xX*N*V;

endfunction
