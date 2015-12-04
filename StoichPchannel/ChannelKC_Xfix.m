function X_fix = ChannelKC_Xfix (X,t)
	 
	 ## usage: X_fix = ChannelKC_Xfix (X,t)
	 ## 
	 ## 

  global STOICH

  const = ThermoConstants;
  F = const.F;

  i_fix = STOICH.i_fix;
  X_0 = STOICH.X_0;

  X_fix = X_0(i_fix);		# Leave the same

  X_fix(1) = F*0.1*sin(2*pi*1e3*t/10);
	 
endfunction
