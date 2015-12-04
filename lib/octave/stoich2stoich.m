function M = stoich2stoich (m, i_fix, j_fix)
	 
  ## usage: M = stoich2stoich (m)
  ## 
  ## 

  if nargin<2
     i_fix = [];
  endif

  if nargin<3
     j_fix = [];
  endif

  ## if nargin<4
  ##    j_zero = [];
  ## endif


  ## Copy
  M = m;
	 
  ## Extract info
  N = m.N;			# Stoichiometric matrix
  r = rank(N);			# Rank
  [n_X,n_V] = size(N);		# Sizes
  M.N = N;

  ## Extract I matrices for chemostat
  I_X = I_cs = I_cd = eye(n_X);
  for i=1:n_X
    if ismember(i,i_fix)
      I_cd(i,i) = 0;
    else
      I_cs(i,i) = 0;
    endif
  endfor
  M.I_cd = I_cd;
  M.I_cs = I_cs;

  ## Extract I matrices for flowstat
  I_V = I_fs = I_fd = eye(n_V);
  for i=1:n_V
    if ismember(i,j_fix)
      I_fd(i,i) = 0;
    else
      I_fs(i,i) = 0;
    endif
  endfor
  M.I_fd = I_fd;
  M.I_fs = I_fs;

  ## Set i_fix rows to zero
  N_cd = I_cd*N;

  ## Set j_fix cols to zero
  N_d = N_cd*I_fd;

  n_fix = length(i_fix);
  ## [n_X,n_V] = size(N);
  ## N(i_fix,:) = zeros(n_fix,n_V);

  ## Set up matrices: chemostat
  N_s = N-N_d;			# Static part
  M.N_d = N_d;
  M.N_s = N_s;

  ## Set up matrices: flowstat
  

  ## Null space stuff
  [G,K] = stoich2null(N_d);
  M.G = G;
  M.K = K;

  ## Range space stuff
  [RR i_x]  = rref (N_d');	# Reduced row echelon form
  n_x = length(i_x);		# Number of reduced-order states

  ## Transformations 
  
  ## Transform x to X
  M.L = M.L_Xx = L_Xx = RR'(:,1:n_x);
  
  ## Transform X to x
  T = eye(n_X);
  M.L_xX = L_xX = T(i_x,:);
  
  ## Transform X to x_d
  I_d = setdiff(1:(n_X),i_x);
  M.L_dX = L_dX = T(I_d,:);
  
  ## For reconstructing X from x
  M.G_X = G_X = L_dX'*G;
  G_X_alt = eye(n_X)-L_Xx*L_xX;	# Alternative computation
  check_norm(G_X-G_X_alt,"G_X");

  ## Chemostat versions: chemostat
   M.G_X_cs = G_X(:,i_fix);
   M.G_X_d = G_X;
   M.G_X_d(i_fix,:) = zeros(n_fix,n_X);

  ## Some dimensions
  [M.n_X M.n_V] = size(N);
  M.n_x = n_x;

  ## Indices of x
  M.i_x = i_x;

  ## Indices of conserved moieties 
  [n_G m_G] = size(G);
  for i=1:n_G
      i_cm{i} = find(G(i,:)!=0);
  endfor
  M.i_cm = i_cm;

endfunction
