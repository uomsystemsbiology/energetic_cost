function [G,K] = stoich2null (N)

  ## usage: [G,K] = stoich2null (N)
  ## Compute G: left null matrix (conserved moieties)
  ## Compute K: right null matrix (pathways)

  tol = 100*eps;

  G = inull(N')';
  K = inull(N);

  if norm(N*K)>tol
     warning("stoich2null: NK not zero")
  endif

  if norm(G*N)>tol
     warning("stoich2null: GN not zero")
  endif

endfunction
