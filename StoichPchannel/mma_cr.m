function kappa = mma_cr (i_cr,FV)
	 
  ## usage: kappa = mma_cr (i_cr,V)
  ## Gating function for voltage modulated gates
  ## Uses Keener Sneyd formulation

## Get some constants
const = ThermoConstants;
F = const.F;
V_N = const.V_N;

##Convert FV to volts
V = FV/F;

  
  one = ones(size(V));

  V_eq = -65; 			# mV (KeeSne 5.1  p 206)
  v = V*1000 - V_eq;		# They use voltage relative to V_eq in
				# mV
  
  if (i_cr==1)		# n HH gate

    ## Hodgkin-Huxley
    HH.alpha = 0.01*(-v+10)./(exp( (-v+10)/10 ) - 1);
    HH.beta  = 0.125*exp(-v/80);

    ## Physical parameters (see CompareModels.m)
    ## V_g =  25.8519960916881;
    ## k_c =  0.0500000000000000;
    ## k_o =  0.00869005086328556;

    V_g = V_N;
    k_c =  5.7537;
    k_o =  1;

  elseif (i_cr==2)	# m HH gate

    ## Hodgkin-Huxley
    HH.alpha = 0.1*(-v+25)./(exp( (-v+25)/10 ) - 1);
    HH.beta  = 4*exp(-v/18);
    
    ## Physical parameters
    ## V_g = 9;
    ## k_c = 0.2;
    ## k_o =  0.00261286753285733;
    ## V_g =  8.61733203056271;
    ## k_c =  0.200000000000000;
    ## k_o =  0.00189597393241766;
    V_g = V_N/3;
    k_c =  105.49;
    k_o =  1;
  elseif (i_cr==3)	# h HH gate
	 
    ## Hodgkin-Huxley
    HH.beta = 0.07*exp(-v/20);
    HH.alpha  = 1./(exp((30-v)/10) + 1);

    ## Physical parameters
    ## V_g = 7;
    ## k_c = 0.05;
    ## k_o =    6.84412974503901e-06;
    ## V_g =  6.46299902292203;
    ## k_c =  0.0500000000000000;
    ## k_o =    3.16405718378869e-06;
    V_g = V_N/4;
    k_c =  1;
    k_o =    6.3281e-05;

  else
    error(sprintf("Gate function i_cr %i does not exist",i_cr));
  endif

  ## HH stuff in msec - change to sec.
  HH.alpha = HH.alpha*1e3;
  HH.beta = HH.beta*1e3;

  ## HH time constant.
  HH.tau = 1./(HH.alpha + HH.beta);
  
  ## Physical
  Ph.alpha = k_c*exp(V/V_g);
  Ph.beta = k_o*one;

  ## Physical time-constant
  Ph.tau = 1./(Ph.alpha + Ph.beta);

  ## Choose kappa to make tau same as HH case.
  kappa = Ph.tau./HH.tau;

endfunction
