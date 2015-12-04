function dx = HH_fun (x,t)

  ## usage:  dx = HH_fun (x,t)
  ##
  ## 

  ## TT is used for fsolve
  global U_HH
    
  par = HH_numpar;		# Parameters
  xx = zeros(8,1);			# Not used
  yy = 0;			# Not used

  u = U_HH;
  dx = HH_ode(x,u,t,par);
endfunction
