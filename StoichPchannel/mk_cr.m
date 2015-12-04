function [dx x_0 alpha beta] = mk_cr (x,v)
  ## usage:  [dx x_0 alpha beta] = mk_cr (x,v)
  ##
  ## Constitutive relation for Mk - potassium memristor
  ## v is voltage across memristor.
  
  V_eq = -65;		# mV Resting potential
  V = -(v - V_eq);	# HH equations are relative to V_eq and negative

## HodHux52 (12)
  vv = V + 10;		
  if (abs(vv)<1e-6)
    alpha = 0.1;
  else
    alpha = 0.01*vv./(exp(vv/10) -1);
  end

## HodHux52 (13)
  beta = 0.125*exp(V/80);

## HodHux52 (7)
  dx = alpha.*(1-x) - beta.*x;

## Steady-state
  x_0 = alpha./(alpha+beta);
endfunction
