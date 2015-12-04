function [x X V V_fix P] = stoich_sim (name,iname,t,X_0,i_fix,J_fix,enzyme_fun,X0_fun,V0_fun);
  
  ## usage:  [x X V] = stoich_sim (name,iname,t,X_0,i_fix,j_fix,enzyme_fun)
  ## Simulation via stoichiometric bond graph formulation
  ## Uses reduced order system.
  ## x reduced-order state
  ## X full-state
  ## name, iname name of system and inverse system
  ## t time vector
  ## X_0 initial state.
  ## i_fix, j_fix, index vectors of fixed species and reaction
  ## OR j_fix can be a structure with fields j_fix and j_zero.
  ## enzyme_fun : name of enzyme-modulation function of form fun(X,t)

  # Copyright (C) 2014 by Peter J. Gawthrop
  
  global STOICH

  ## Default to no fixed states.
  if nargin<5
     i_fix = [];
  endif

  if nargin<6
     J_fix.j_fix = [];
     J_fix.j_zero = [];
  endif

  if nargin<7
    enzyme_fun = "";
  endif

 if nargin<8
    X0_fun = "";
  endif

if nargin<9
    V0_fun = "";
  endif

 ## Extract the J_fix parameters
 if isstruct(J_fix)
   j_fix = J_fix.j_fix;
   j_zero = J_fix.j_zero;
 else
   j_fix = J_fix; 
   j_zero = [];
 endif

 ## Preset parameters
  if !isfield (STOICH, "par")
    par = eval(sprintf("%s_numpar;",name));
  else
    par = STOICH.par;
    ## warning("using preset parameters");
  endif

  ## Preset thermodynamic parameters
  presetK_c = isfield (STOICH, "K_c");
  if presetK_c
    K_c = STOICH.K_c;
  endif

  ## Set up stoichimetric matrices
  STOICH = dm2stoich (name,iname);
  N = STOICH.N;
  [n_X,n_V] = size(N);
  STOICH.name = name;
  
  STOICH.enzyme_fun = enzyme_fun;
  STOICH.X0_fun = X0_fun;
  STOICH.V0_fun = V0_fun;
  STOICH.par = par;

  ## Put in the fixed states (indexed by i_fix)
  ##STOICH = stoich2stoich(STOICH,i_fix,j_fix,j_zero);
  if (length(j_zero)==0)
    # Use N^cd -- assume that flowstats non-zero
    STOICH = stoich2stoich(STOICH,i_fix); 
  else
    # Use N^d 
    STOICH = stoich2stoich(STOICH,i_fix,j_zero);
  endif

  ## n_fix = length(i_fix);
  ## if n_fix>0
  ##   Z = zeros(n_fix,n_V);
  ##   N_fix = STOICH.N;
  ##   N_fix(i_fix,:) = Z;
  ##   STOICH.N = N_fix;
  ##   STOICH = stoich2stoich(STOICH);
  ## endif

  ## Sanity check
  if (length(X_0)!=n_X)
    error(sprintf("X_0 should have %i elements, not %i", n_X, length(X_0)));
  endif
  
  STOICH.X_0 = X_0;		# State

  n_x = STOICH.n_x;
  i_free = setdiff(1:n_X,i_fix);
  j_free = setdiff(1:n_V,j_fix);

  ##Set up I_fix - diagonal with zeros for fixed elements.
  STOICH.i_fix = i_fix;
  I_fix = zeros(n_X,n_X);
  for i=1:n_X
    if ismember(i,i_fix);
       I_fix(i,i) = 1;
    endif
  endfor
  ##STOICH.I_free = eye(n_X)-I_fix;
  STOICH.i_fix = i_fix;
 
  ## Set up J_fix - diagonal with zeros for fixed elements.
  STOICH.j_fix = j_fix;
  J_fix = zeros(n_V,n_V);
  for j=1:n_V
    if ismember(j,j_fix);
       J_fix(j,j) = 1;
    endif
  endfor
  STOICH.J_free = eye(n_V)-J_fix;

  ## Reduced-order initial state
  x_0 = STOICH.L_xX*X_0;

  if (length(t)>1)
    ## Reduced-order sim
    x = lsode("stoich_sim_fun",x_0,t)';
  elseif (length(t)==1)
   ##disp(sprintf("Find steady-state evaluated with time = %g",t));
   STOICH.t = t;
   OPT.AutoScaling = "off";
   OPT.TolX = 1e-14;
   [x fvec info] = fsolve("stoich_sim_fun",x_0,OPT);
   if (info!=1)
     warning(sprintf("fsolve has terminated with info=%i",info));
   endif
  else
    error("t must not be empty");
  endif

  ## x_i = x_0;
  ## x = [];
  ## dt = mean(diff(t))
  ## nt = length(t)
  ## for i = 1:nt
  ##     x = [x x_i];
  ##     x_i = x_i + dt*stoich_sim_fun (x_i,t(i));
  ##     ## Avoid numerical issues !!??
  ##     small = 1e-8;
  ##     i_neg = find(x_i<small);
  ##     x_i(i_neg) = small;
  ## endfor

  if nargout>1
    ## Recreate full state vector and flows
    X = []; V = []; V_fix = [];
    for i=1:length(t);
	
      ## State vector
      x_i = x(:,i);
      X_i =  STOICH.L*x_i + STOICH.G_X*X_0;
      X = [X X_i];
      if nargin>2
	## Flows
	t_i = t(i);
	[dx_i,V_i] = stoich_sim_fun(x_i,t_i);
	V = [V V_i];
      endif
    endfor
  endif

  if nargout>3
    ## Flows needed to fix the state
    dX = N*V;
    V_fix = -dX(i_fix,:);
  endif

  if nargout>4			# Do energy stuff

 ## Preset thermodynamic parameters
    if !presetK_c
      ## Old method
      ## Evaluate free-energy constants.
      par = STOICH.par;
      eval(sprintf("[K_c i_lin] = %s_K(par);",name));
    endif

    P.K_c = K_c;
    
    ## Get RT value
    ## sym = eval(sprintf("%s_sympar;",name));
    ## RT = par(sym.RT);
    RT = 1;			# Put this in later
    
    ## Chemical potentials
    mu = RT*log(diag(K_c)*X);
    
    ## ## Correct for the linear versions
    if exist("i_lin")
      mu(i_lin,:) = RT*diag(K_c(i_lin))*X(i_lin,:);
    endif
    P.mu = mu;
       
    ## Affinity
    A = -N'*mu;
    P.A = A;

    ## Dissipated power (in Re)
    P_r = A.*V;
    P.P_r = P_r;

    ## Dissipated power (in Re): internal
    P_ri = P_r(j_free,:);
    P.P_ri = P_ri;

    ## Dissipated power (in Re): external
    P_re = P_r(j_fix,:);
    P.P_re = P_re;


    ## Capacitive power
    P_c = mu.*dX;
    P.P_c = P_c;

    ## Capacitive power - internal
    P_ci = P_c(i_free,:);	
    P.P_ci = P_ci;

    ## Capacitive power - external
    P_ce = P_c(i_fix,:); 	# External flows are -dX
    P.P_ce = P_ce;

  endif
  
endfunction
