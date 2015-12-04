function [V_N R F] = RTF (T,scale)

  ## usage:  ratio = RTF (T)
  ## Computes RT/F
  ## 

  if nargin<1
    T = 300;			# K
  end

  if nargin<2
    scale = 1;		
  end

  R = 8.3144621;		# Gas constant J K^-1 mol^-1
  F = 9.64853399*10^4/scale;	# Faraday constant C mol^-1

  V_N = (R*T)/F;

endfunction
