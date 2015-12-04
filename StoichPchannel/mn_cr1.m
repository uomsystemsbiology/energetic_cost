function [dx x_0 alpha beta] = mn_cr1 (x,v)

  ## usage:  dx = mn_cr1 (x,v)
  ##
  ## 

  V_eq = -65;		# mV Resting potential
  V = -(v - V_eq);	# HH equations are relative to V_eq and negative

   ## HodHux52 (20)
  vv = V + 25;
  if (abs(vv)<eps)
    alpha = 1;
  else
    alpha = 0.1*vv./(exp(vv/10) -1);
  end
  
  ## HodHux52 (21)
  beta = 4*exp(V/18);

  ## HodHux52 (15)
  dx = alpha.*(1-x) - beta.*x;

  ## Steady-state
  x_0 = alpha./(alpha+beta);

endfunction
