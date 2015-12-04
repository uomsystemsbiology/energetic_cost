function ind = name2index (names)
	 
  ## usage: ind = name2index (names)
  ## 
  ## names is a cell array of strings
	 
  for i = 1:length(names)
      eval(sprintf("ind.%s = %i;",names{i},i));
  endfor
  
endfunction

