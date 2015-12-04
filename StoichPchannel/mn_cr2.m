function [dx x_0 alpha beta] = mn_cr2 (x,v)

  ## usage:  dx = mn_cr2 (x,v)
  ##
  ## 

  V_eq = -65;		# mV Resting potential
  V = -(v - V_eq);	# HH equations are relative to V_eq and negative


  ## HodHux52 (23)
  alpha = 0.07*exp(V/20);

  ## HodHux52 (24)
  vv = V + 30;
  if (abs(vv)<eps)
    beta = 0.5;
  else
    beta = 1./(exp(vv/10)+1);
  end
  
  dx = alpha.*(1-x) - beta.*x;

  x_0 = alpha./(alpha+beta);
endfunction
