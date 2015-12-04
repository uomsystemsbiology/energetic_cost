function str = mat2latex (A,name)
	 
  ## usage: str = mat2latex (A,name)
  ## 
  ## Converts matrix A   into a LaTeX string with optional name etc.
	 
  str = strrep(strrep(mat2str(A)," ","&"),";","\\\\\n")(2:end-1);
  
  if nargin>1
    str = sprintf("%s =\n\\begin{pmatrix}\n%s\n\\end{pmatrix}\n", name, str);
  endif

endfunction
