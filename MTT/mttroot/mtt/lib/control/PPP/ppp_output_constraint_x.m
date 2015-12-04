function [Gamma,Gamma_x,gamma] = ppp_output_constraint_x (A,B,C,D,A_u,Tau,Min,Max,Order,i_y)

  ## usage:  [Gamma,gamma] = ppp_output_constraint (A,B,C,D,A_u,Tau,Min,Max,Order)
  ##
  ## Derives the output constraint matrices Gamma and gamma_o and
  ## gamma_x
  ## Avoids the on-line usage of ppp_output_constraint
  ## For Constraints Min and Max at times Tau
  ## Order=0 - output constraints
  ## Order=1 - output derivative constraints
  ## etc
  ## NOTE You can stack up Gamma and gamma matrices for create multi-output constraints.

  ## Copyright (C) 2007 by Peter J. Gawthrop
  ## 	$Id: ppp_output_constraint_x.m,v 1.1 2007/10/09 08:58:15 peterg Exp peterg $	

  ## Sizes
  [n_x,n_u,n_y] = abcddim(A,B,C,D); # System dimensions
  [n_U,m_U] = size(A_u);	# Number of basis functions
  [n,n_tau] = size(Tau);		# Number of constraint times
  
  if n_tau==0			# Nothing to be done
    Gamma = [];
    Gamma_x = [];
    gamma_0 = [];
    return
  endif

  ## Defaults
  if nargin<10
    Order = zeros(1,n_tau);
  endif

  if nargin<11
    i_y = 1;			# First output
  endif

  if n != 1
    error("Tau must be a row vector");
  endif
  
  n = length(Min);
  m = length(Max);
  o = length(Order);

  if (n != n_tau)||(m != n_tau)||(o != n_tau)
    error("Tau, Min, Max and Order must be the same length");
  endif
  

  ## Compute Gamma
  Gamma = [];
  zero_x = zeros(n_x,1);
  for i=1:n_U
    one = [1:n_U]'==i*ones(n_U,1)
    U = one;			# Only ith element non-zero
    y_i = ppp_ystar (A,B,C,D,zero_x,A_u,U,Tau); # Compute y* for ith input for each tau
    y_i = y_i(i_y,:); # Pluck out output i_y
    Gamma_i = [];
    if (Min>-inf)
      Gamma_i = [Gamma_i; -y_i']; # Min part of column
    endif
    if (Max<inf)
      Gamma_i = [Gamma_i;  y_i'];	# Add column for Max
    endif
    Gamma = [Gamma Gamma_i]; # Put in parts for Min and max
  endfor
  
  ## Compute Gamma_x
  Gamma_x = [];
  zero_U = zeros(n_U,1);
  for i=1:n_x
    one = [1:n_x]'==i*ones(n_x,1)
    x_0 = one;			# Only ith element non-zero
    y_i = ppp_ystar (A,B,C,D,x_0,A_u,zero_U,Tau); # Compute y* for ith state for each tau
    y_i = y_i(i_y,:); # Pluck out output i_y
    Gamma_i = [];
    if (Min>-inf)
      Gamma_i = [Gamma_i; -y_i']; # Min part of column
    endif
    if (Max<inf)
      Gamma_i = [Gamma_i;  y_i'];	# Add column for Max
    endif
    Gamma_x = [Gamma_x Gamma_i]; # Put in parts for Min and max
  endfor


  ## Compute gamma
#   zero_U = zeros(size(U));
#   y_x = ppp_ystar (A,B,C,D,zero_x,A_u,zero_U,Tau);	# Output if U is zero
  gamma = [];
  if (Min>-inf)
    gamma = [gamma; -(Min)'];
  endif
  if (Max<inf)
    gamma = [gamma; (Max)'];
  endif
  
endfunction


