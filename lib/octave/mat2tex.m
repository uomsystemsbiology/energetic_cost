function tex = mat2tex (mat)

  ## usage:  tex = mat2tex (mat)
  ##
  ## Convert matrix mat to latex form.

  if isscalar(mat)
    tex = mat2str(mat);
  else
    tex = strrep(strrep(mat2str(mat)," "," & "),";"," \\\\\n")(2:end-1);
  endif
  tex = strrep(tex, "-0 ", "0 ");
endfunction
