function K_c = Keq2K (K_eq,N,K,G, k_c_1,tol)
	 
## usage: K_c = Keq2K (K_eq,N,K,G [,k_c_1])
## 
## 
## Compute BG C parameters K_c from equilibrium constants K_eq.
## NB K_eq must be thremodynamically consistent.

  if nargin<6
    tol = 1e-3;			# Error tolorance.
  endif

  logK_eq = log(K_eq);

  ##First check that Keq is thermodynamically consistent.
  check = norm(K'*logK_eq)/norm(logK_eq);
  if (check>tol)
     warning(sprintf("Thermodynamically inconsistent K_eq: check = %g",check));
  endif
  NN = -N';			# Transformation from mu to A
  pNN = pinv(NN);		# Pseudo inverse

  logK_c = pNN*logK_eq;		# Use pseudo inverse to solve for log K_c
  K_c = exp(logK_c);		# Solution

  ## Sanity check - reverse the solution
  err = norm(NN*logK_c - logK_eq)/norm(logK_eq);
  if (err>tol)
    warning(sprintf("Incorrect solution: err = %g",err));
  endif

  if nargin>4
    [n_g n_x] = size(G);
    n_k = length(k_c_1);
    if n_k!=n_g
       error(sprintf("k_c_1 should be length %i, not %i",n_g,n_k));
    else
	K_c_1 = exp(G'*log(k_c_1));
	K_c = K_c.*K_c_1;
    endif
    
  endif
  

endfunction
