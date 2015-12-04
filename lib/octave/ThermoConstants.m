function const = ThermoConstants (T)

  ## Computes some useful constants

  if nargin<1
    T = 300;			# K
  end

  R = 8.3144621;		# Gas constant J K^-1 mol^-1
  F = 9.64853399*10^4;		# Faraday constant C mol^-1
  RT = R*T;

  V_N = (RT)/F;	

  A = 	6.02214129e23; 		# Avogadro's constant mol^-1

  const.R = R;
  const.RT = RT;
  const.F = F;
  const.V_N = V_N;
  const.A = A;

endfunction
