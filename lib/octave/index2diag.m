function D = index2diag (ind,n)
  ## usage: D = index2diag (ind)
  ## D is diagonal with unit entries corresp. to ind.
  ## 
  for i=1:n
    d(i) = ismember(i,ind);
  endfor
  D = diag(d);
endfunction

