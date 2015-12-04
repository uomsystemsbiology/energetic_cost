function info = struc_info (name,par)
	 
  ## usage: info = struc_info (name,par)
  ## 
  ## Finds out various structural info
  ## Use <struc_info name input>  (shell script) to generate name_input_info.m
	 
	 
  if nargin<2
    par.signal="input";
  endif
  
  
  ## Causality
  info = eval(sprintf("%s_%s_info;", name, par.signal));
  info.is_effort = strcmp(info.cause,"effort");
  info.is_flow = strcmp(info.cause,"flow");
  info.i_effort = find(info.is_effort);
  info.i_flow = find(info.is_flow);
endfunction
