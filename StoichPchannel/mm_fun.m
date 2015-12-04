function val = mm_fun (FV,V_0,k)

  ## usage:  val = mm_fun (V,V_0,k)
  ##
  ## V membrane voltage
  ## V_0: if V_0 == 0 then use MM equation
  ##      if V_0 != 0 then use HH equation

  if nargin<2
    V_0 = 0;
  endif

  if nargin<3
    k = 1
  endif


  ## Get some constants
  const = ThermoConstants;
  F = const.F;
  V_N = const.V_N;

  ##Convert FV to volts
  V = FV/F;
  V_bar = V/V_N;		# Normalised voltage
  V_tilde = V_bar-V_0/V_N;	# Differential voltage (V_0=0 for MM)

  ## MM computation 
  val = 1./(exp(V_tilde) + k);

endfunction
