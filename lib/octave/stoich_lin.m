function [SYS sys comp m_0] = stoich_lin (stoich,normalise)
	 
  ## usage: [SYS sys comp m_0] = stoich_lin (stoich)
  ## 
  ## SYS: full-order system
  ## sys: reduced-order system
  ## comp: linearised component values
  ## m_0: stoichiometric matrices with chemostats

  if nargin<2
     normalise = 0;
  endif
	 
  ## Extract from stoich data structure
  n_X = stoich.n_X;
  n_x = stoich.n_x;
  n_V = stoich.n_V;
  N = stoich.N;			# Stoichiometric matrix
  N_f = stoich.N_f;		# Forward stoichiometric matrix
  N_r = stoich.N_r;		# Backward stoichiometric matrix
  K = stoich.K_c;		# Thermodynamic constants
  kappa = stoich.kappa;		# Reaction gains (=v_max/k_v)
  
  if isfield(stoich,"k_v")
    k_v = stoich.k_v;
    rho_v = stoich.rho_v;
    MassAction = 0;
  else
    MassAction = 1;
  endif

  X_ss = stoich.X_ss;

  if isfield(stoich,"i_fix")
    i_fix = stoich.i_fix;		# Fixed (ie external) state indices
  else
    i_fix = [];
  endif

  if isfield(stoich,"j_fix")
    j_fix = stoich.j_fix;		# Fixed (ie external) flows
  else
    j_fix = [];
  endif

  if isfield(stoich,"module")
    module = stoich.module;
  else
    module = 0;
  endif

  ## Sanity checks
  if (length(X_ss)!=n_X)
     error(sprintf("X_ss should have %i elements, not %i",n_X,length(X_ss)));
  endif
  if (length(K)!=n_X)
     error(sprintf("K should have %i elements, not %i",n_X,length(K)));
  endif
  if (length(kappa)!=n_V)
     error(sprintf("kappa should have %i elements, not %i",n_V,length(kappa)));
  endif

  ## Non-fixed state indices
  ##[n_X, n_V] = size(N);		# dimensions
  i_free = setdiff([1:n_X],i_fix);
  j_free = setdiff([1:n_V],j_fix);
  
  ## Steady-state flows
  mu_ss = log(K.*X_ss);		# Chemical potential
  Vplus0  = exp(N_f'*mu_ss);	# Forward nominal reaction flow
  Vminus0 = exp(N_r'*mu_ss);	# Backward nominal reaction flow

  if MassAction
    Den = ones(n_V,1);
  else
    Den = (1 + ((1-rho_v).*Vplus0 + rho_v.*Vminus0)./k_v );
  endif

  v_ss = kappa.*(Vplus0 - Vminus0)./Den; # Mass-action reaction flow
    
  ## Linearise about steady state

  ## Derivatives of affinities wrt X
  dA_f = N_f'/diag(X_ss);
  dA_r = N_r'/diag(X_ss);

  ## Derivatives of Vplus0,  Vminus0 wrt X
  dVplus0  =  diag(Vplus0)*dA_f;
  dVminus0  = diag(Vminus0)*dA_r;

  ## ## ## Derivatives of Vplus0,  Vminus0 wrt X
  ## check_norm(dVplus0-diag(Vplus0)*N_f'/diag(X_ss), "dVplus0",0)
  ## check_norm(dVminus0-diag(Vminus0)*N_r'/diag(X_ss), "dVminus0",0)

  ## Full-order system
  C = diag(kappa./Den)*(dVplus0 - dVminus0);
  if !MassAction
     C = C - diag(v_ss./Den)*(diag((1-rho_v)./k_v)*dVplus0 + diag(rho_v./k_v)*dVminus0);
  endif
  if normalise
    D = eye(n_V);
  else
    D = diag(v_ss./kappa);
  endif
  
  ## Create N_d: N with zero rows for fixed states
  N_d = N;
  N_d(i_fix,:) = zeros(length(i_fix),n_V);
  ##N_d(:,j_fix) = zeros(n_X,length(j_fix));
  A = N_d*C;
  B = N_d*D;
  SYS = ss(A,B,C,D);
  
  ## Reduced-order system
  if nargout>1

    ## Find the relevant transformation matrixes
    m_0.N = N_d;
    m_0 = stoich2stoich(m_0,i_fix);
    L_xX = m_0.L_xX;
    L_Xx = m_0.L_Xx;
    G_X_cs = m_0.G_X_cs;
    N_d = m_0.N_d;
    
    if (length(i_fix)==0)		# Old version
      a = L_xX*A*L_Xx;
      b = L_xX*B;
      c = c_out = C*L_Xx;
      d = D;
    else			# New version
      n_fix = length(j_fix);
      c = c_0 = C*L_Xx;
      c_0(j_fix,:) = zeros(n_fix,m_0.n_x); # Zap j_fix rows
      d_x = C*G_X_cs;		# Chemostat
      d_v = zeros(n_V,n_fix);	# Flowstat
      for jj = 1:n_fix
	d_v(j_fix(jj),jj) = 1;
      endfor

      d = [d_x d_v];
      a = L_xX*N_d*c_0;
      b = L_xX*N_d*d;

      if module
	## Species flows as chemostat outputs
	N_s = N - N_d;
	c_cs = N_s(i_fix,:)*c;
	d_cs = N_s(i_fix,:)*d;
	
	## Reaction affinities as flowstat outputs
	dAdX = dA_r-dA_f;		# dA/dX
	dAdx = dAdX*L_Xx;		# dA/dx
	c_fs = dAdx(j_fix,:);
	c_out = [c_cs;c_fs];

        [n_stat junk] = size(c_out);
	d_fs = zeros(length(j_fix),n_stat);
	d_out = [d_cs;d_fs];

	io_name = {stoich.species{i_fix} stoich.reaction{j_fix}};
      else
	c_out = c;
	d_out = d;
      endif
    endif
    if exist("io_name")
      sys = ss(a,b,c_out,d_out,"INNAME",io_name,"OUTNAME",io_name);
    else
      sys = ss(a,b,c_out,d_out);
    endif
  endif

  ## Component linearisation.
  if nargout>2
     if normalise
       comp.gamma = ones(length(kappa),1);
     else
       comp.gamma = v_ss./kappa;
     endif
    comp.K = 1./X_ss;

    if MassAction
      Den = ones(n_V,1);
    else
      Den = (1 + ((1-rho_v).*Vplus0 + rho_v.*Vminus0)./k_v );
    endif

    kappa_f = kappa.*Vplus0./Den;
    kappa_r = kappa.*Vminus0./Den;

    if !MassAction
       kappa_f = kappa_f - k_v.\(1-rho_v).*(v_ss./Den).*Vplus0;
       kappa_r = kappa_r + k_v.\rho_v    .*(v_ss./Den).*Vminus0;
   endif
    comp.kappa_f = kappa_f;
    comp.kappa_r = kappa_r;
    comp.kappa_0 = (kappa_f + kappa_r)/2;
    comp.kappa_1 = (kappa_f - kappa_r)/2;
  endif
endfunction
