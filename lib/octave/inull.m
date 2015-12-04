function Z = inull (A)
	 
  ## usage: Z = inull (A)
  ## Z is a null matrix with integer entries
  ## Used for stoichiometric analysis
  ## See stoich2null

  tol = 1000*eps;

  [m,n] = size(A);	 
  [R,pivcol] = rref(A);

  r = length(pivcol);
  nopiv = 1:n;
  nopiv(pivcol) = [];
  Z = zeros(n,n-r,class(A));
  if n > r
    Z(nopiv,:) = eye(n-r,n-r,class(A));
    if r > 0
      Z(pivcol,:) = -R(1:r,nopiv);
    endif
  endif

  ## Make it integer (??)
  minval = min(min(abs(Z)+(abs(Z)<tol)));
  ## if (abs(minval-1)>tol)
  ##    warning(sprintf("inull: min value = %g, multiplying by %g",minval,1/minval));
  ##    Z = Z/minval;
  ## endif
  
  if (abs(minval-1)>tol)
    warning(sprintf("inull: min value = %g",minval));
   endif

  discrepancy = norm(Z-round(Z));
  if discrepancy>tol
     warning("inull: non-integer matrix returned - discrepancy %g", discrepancy);
  endif

endfunction
