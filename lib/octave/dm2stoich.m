function [m tex var_names] = dm2stoich (sys_name,dsys_name,X_0)

  ## usage:  [stoich tex var_names] = dm2stoich (sys_name,dsys_name)
  ## stoich.N Stoichiometric matrix
  ## 
	 
  ## Generate stoichiometric matrix etc.

  ##Defaults

  ## Don't care version
  X_00 = ones(1000,1);

  if nargin<3
    X_0 = X_00;
  endif
  
  ## Function names
  name_dm = sprintf("%s_dm",sys_name);
  name_struc = sprintf("%s_struc",sys_name);

  eval(sprintf("par = %s_numpar;",sys_name));
  ## Descriptor matrices for system in integral causality
   if exist(name_dm)
    eval(sprintf("[A,B,C,D] = %s(par,X_0);", name_dm));
  else
    error(sprintf("%s.m does not exist. Please execute ""mtt %s dm m""", name_dm,sys_name))
  endif

  ## Names for system in integral causality
  if exist(name_struc)
    eval(sprintf("[input_name, output_name,state_name] = %s;", name_struc));
    i_name = cellstr(long2short(mtt_strip_name(input_name,sys_name)));
    o_name = cellstr(long2short(mtt_strip_name(output_name,sys_name)));
    s_name = cellstr(long2short(mtt_strip_name(state_name,sys_name)));
  else
    error(sprintf("%s.m does not exist. Please execute ""mtt %s struc m""", name_struc,sys_name))
  endif

  ## Are we using the new form of Re with separate v_f v_r?
  [n_V,n_X] = size(C);	# Induced flows (number of Re)
  [n_XX,n_VV] = size(B);	# Generated flows (number of Re OR 2x number of Re
  single_v = (n_VV==n_V);
  double_v = (n_VV==(2*n_V));
  if (single_v==double_v)
     error("Wrong number of flows: double_v==single_v")
  endif

  ## If no inverse system; generate N from BG and rest numerically
  if (length(dsys_name)==0)
    ##warning("dm2stoich: No inverse system, using stoich2stoich");
    if double_v
      m.N_f = -B(:,1:2:n_VV);	# Odd cols
      m.N_r =  B(:,2:2:n_VV);	# Even cols
      m.N = -m.N_f + m.N_r;
      m.reaction = i_name(1:2:n_VV,:);	# Save reaction names
    else
      m.N = B;			# N Stoichiometric matrix
      m.N_f = (m.N<0).*abs(m.N);
      m.N_r = (m.N>0).*abs(m.N);
      m.reaction = i_name;	# Save reaction names
    endif
    m.N_fr = [m.N_f m.N_r];
    m = stoich2stoich (m);	# Compute null space matrices etc.
    m.species = s_name;		# Save species names
    return
  endif

  m.species = s_name;		# Save species names
  m.reaction = i_name;	# Save reaction names
  
  ## Function names
  dname_dm = sprintf("%s_dm",dsys_name);
  dname_struc = sprintf("%s_struc",dsys_name);

  ## Descriptor matrices for system in derivative causality
  if exist(dname_dm)
    eval(sprintf("[AA,BB,CC,DD,EE] = %s(par,X_00);", dname_dm));
  else
    error(sprintf("%s.m does not exist. Please execute ""mtt %s dm m""", dname_dm,dsys_name))
  endif
 ## Names for system in integral causality
  if exist(name_struc)
    eval(sprintf("[input_name, output_name,state_name] = %s;", name_struc));
    i_name = cellstr(mtt_strip_name(input_name,sys_name));
    o_name = cellstr(mtt_strip_name(output_name,sys_name));
    s_name = cellstr(mtt_strip_name(state_name,sys_name));
  else
    error(sprintf("%s.m does not exist. Please execute ""mtt %s struc m""", name_struc,sys_name))
  endif


  ## Names for system in derivative causality
  if exist(dname_struc)
    eval(sprintf("[input_name, output_name,state_name,nonstate_name] = %s;", dname_struc));
    di_name = cellstr(mtt_strip_name(input_name,dsys_name));
    do_name = cellstr(mtt_strip_name(output_name,dsys_name));
    ds_name = cellstr(mtt_strip_name(state_name,dsys_name));
    dn_name = cellstr(mtt_strip_name(nonstate_name,dsys_name));
  else
    error(sprintf("%s.m does not exist. Please execute ""mtt %s struc m""", dname_struc,dsys_name))
  endif

  ## Dimensions from inverse system
  [NN,MM] = size(EE);
  n_z = NN-rank(EE); 		# Non-states
  n_x = NN - 2*n_z;		# States
  n_X = n_x + n_z;		# Nonstates + states

  [n_V,n_u] = size(DD);		# Rates and independent rates.
  n_v = n_V-n_u;		# Dependent rates

  ## Other dimensions from forward system.
  I_ext = find(sum(abs(D))>0);	# Indices of external flows.
  I_int = setdiff([1:n_V],I_ext);
  n_v_e = length(I_ext);
  n_v_i = length(I_int);

  ## Save some dimensions for Actual system
  m.n_x = n_z;			# Actual
  m.n_X = n_X;
  m.n_z = n_x;			# Actual
  m.n_V = n_V;
  m.n_v = n_v;
  m.n_u = n_u;

 ## Create matrices relating dependent and independent flows to all
  ## flows.
  ## Assumes that do_list contains all flows and di_list the independent
  ## flows.
  ## Transformation from v to v_i (T_vi) and v to v_d (T_vd)
  I_i = [];		
  for i = 1:n_u
    for j = 1:n_V
      if strcmp(di_name{i},o_name{j})
	I_i = [I_i j];
      endif
    endfor
  endfor
  I_d = setdiff(1:n_V,I_i);
  T = eye(n_V);
  T_vi = T(I_i,:);		
  T_vd = T(I_d,:);

  ## Matrix to convert between original and inverse output list.
  I_o = [];		
  for i = 1:n_V
    for j = 1:n_V
      if strcmp(do_name{i},o_name{j})
	I_o = [I_o j];
      endif
    endfor
  endfor
  T = eye(n_V);
  T_o = T(:,I_o);

  ## Sanity check
  if (rank(T_o)!=n_V)
     disp("system reaction names");
     o_name
     disp("inverse system reaction names");
     do_name
     error("Reaction names of system and inverse do not match: see above");
  endif

  ## Create matrices relating all states to independent and dependent states
  ## Assumes that s_list contains all states and dn_list the dependent states
  I_x = [];		
  for i = 1:n_x
    for j = 1:(n_x+n_z)
      if strcmp(ds_name{i},s_name{j})
	I_x = [I_x j];
      endif
    endfor
  endfor
  I_z = setdiff(1:(n_x+n_z),I_x);
  T = eye(n_x+n_z);
  L_xX = T(I_z,:);
  L_dX = T(I_x,:);

  ## Convert C and D matrices from derivative system to
  ## original coordinate system
  CCC = T_o*CC;
  DDD = T_o*DD;

  ## Create matrix relating dependent to independent states
  L_dx = AA(1:n_x,n_X+1:n_X+n_z);
 
  ## The matrices
  N = B;			# N Stoichiometric matrix

  ## Decompose N.
  N_i = N(:,I_int);		# Reaction flows
  N_e = N(:,I_ext);		# Input flows
  
  ## Forward and reverse components of N
  N_f = (N<0).*abs(N);
  N_r = (N>0).*abs(N);

  ## Forward and reverse components of N_i
  N_if = (N_i<0).*abs(N_i);
  N_ir = (N_i>0).*abs(N_i);

  ## Forward and reverse components of N_e
  N_ef = (N_e<0).*abs(N_e);
  N_er = (N_e>0).*abs(N_e);

  K_Vu = K = DDD;			# K (from v_i to v)
  K_du = K_d = T_vd*K;			# K_d (from v_i to v_d)
  K_Vx = M = CCC(:,n_X+1:n_X+n_z);	# M (from zdot to v)
  K_dx = M_d = T_vd*M;    		# M_d (from zdot to v_d) 
  L_Xx = L = N*M;			# L (from z to entire state)
  G = L_dX - L_dx*L_xX;			# Left null space matrix
  G_X = L_dX'*G;			# For reconstructing X from x

  ##N_C_R
  iN_C_R = CCC(n_u+1:n_V,n_X+1:n_X+n_z);


  ## Reduced-order system
  #N_d = round(inv(M_d));
  #N_i = N_d*K_d;

  ## Sanity checks
  tol = 1e-14;
  if norm(N*K)>tol
    warning("NK is not zero");
  endif

  if norm(G*N)>tol
    warning("GN is not zero");
  endif

  ## Linearisation
  J_VX = C;			# dV/dX
  J_XX = N*J_VX;		# dXdot/dX
  J_Vx = J_VX*L_Xx;		# dV/dx
  J_Xx = N*J_Vx;		# dXdot/dx
  J_xx = L_xX*J_Xx;		# dxdot/dx

  J_vX = J_VX(1:n_v,:);		# dv/dX
  J_vx = J_Vx(1:n_v,:);		# dv/dx

  names = {"N","K","L","G","M","K_d","M_d","G_X", ...
	   "N_f","N_r", ...
	   "N_i", "N_if","N_ir", ...
	   "N_e", "N_ef","N_er", ...
	   "T_vi","T_vd","T_o", ...
	   "L_xX","L_dX","L_dx","L_Xx", ...
	   "K_Vx","K_dx","K_du","K_Vu", ...
	   "I_i" "I_d", ...
	   "J_VX" "J_XX" "J_Vx" "J_Xx" "J_xx", ...
	   "J_vX" "J_vx"};

  for i = 1:length(names)
    name = names{i};
    cmd = sprintf("mat = %s;", name);
    eval(cmd);
    ## Replace zero-dimension matrices by scalar zero.
    if length(mat(:))==0
      mat = 0;
    end
    cmd = sprintf("m.%s = mat;", name);
    eval(cmd);
     cmd = sprintf("tex.%s = mat2tex(%s,""%s"",sys_name);",name,name,name);
    eval(cmd);
  endfor

  ## Save variable names
  var_names.X = s_name;
  var_names.x = var_names.X(I_z,1);
  var_names.X_d = var_names.X(I_x,1);

  var_names.V = i_name;
  var_names.v = var_names.V(I_i,1);
  var_names.V_d = var_names.V(I_d,1);

endfunction
