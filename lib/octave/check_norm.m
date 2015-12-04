function check_norm (val,name,small)
## usage: check_norm (val[,name,small])
## Give a warning if norm(val)>small
## (small defaults to 1e-6)

  if nargin<2
    name = '.';
  endif

  if nargin<3
     small = 1e-6;
  endif

 if norm(val)>small
    warning(sprintf("norm(%s_error) = %g",name,norm(val)));
 endif
endfunction
