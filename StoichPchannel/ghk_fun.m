function val = ghk_fun (FV,V_0)

  ## usage:  val = ghk_fun (V,V_0)
  ##
  ## V membrane voltage
  ## V_0: if V_0 == 0 then use GHK equation
  ##      if V_0 != 0 then use HH equation

  if nargin<2
    V_0 = 0;
  endif

  ## Get some constants
  const = ThermoConstants;
  F = const.F;
  V_N = const.V_N;

  ##Convert FV to volts
  V = FV/F;
  V_bar = V/V_N;		# Normalised voltage
  V_tilde = V_bar-V_0/V_N;	# Differential voltage (V_0=0 for GHK)

  ## GHK computation (avoiding 0/0)
  tol = 1e-6;
  den = exp(V_tilde)-1;
  if abs(den)<tol
    val = 1;
  else
    val = V_tilde./den;
  endif
   
endfunction
