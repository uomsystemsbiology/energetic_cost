function tex = mat2texmat (mat, matname,label)

  ## usage:  tex = mat2texmat (mat,name)
  ##
  ## Produce matrix suitable for use in align.

  if nargin<3
    label = matname;
  else
    label = sprintf("%s_%s",label,matname);
  end

  texmat = mat2tex(mat);

  tex = sprintf("%s &= \n\t\\begin{pmatrix}\n%s\n\t\\end{pmatrix}",matname,texmat);
  tex = sprintf("%s \\label{eq:%s}\\\\\n",tex, label);
endfunction